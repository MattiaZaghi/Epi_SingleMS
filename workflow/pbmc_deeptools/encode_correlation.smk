# Snakefile -- ENCODE bulk vs scPBMC pseudobulk correlation (corrplot), on Condor.
# Max parallelization: one bigwig job per ENCODE BAM (5 parallel), then a single
# multiBigwigSummary, then the R corrplot. Reuses the exact pipeline bamCoverage
# command for the ENCODE side; scPBMC pseudobulk bigwigs are now taken from the
# per-cluster, QC-passed-cells tracks built by
# cluster_tracks/Snakefile_cluster_tracks.smk (a.k.a. deeptools_tracks.smk),
# instead of the old whole-possorted_bam bam_to_bw output.
#
# Launch (from this analysis folder so its .snakemake lock is separate):
#   cd /date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/encode_correlation
#   snakemake -s scripts/Snakefile --profile htcondor --jobs 200 --cores 200
import glob, os

BIN   = "/date/gcb/gcb_MZ/nanoCTAR_YD/MS/.snakemake/conda/67977e66f74713228faaa1470a4bdd30_/bin"
BASE  = "/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS"
REFS  = BASE + "/Refs"
OUT   = BASE + "/encode_correlation"
BW    = OUT + "/bigwigs"
RSCRIPT = "/home/mattia/Multi_nanoCTRNA/workflow/pbmc_deeptools/scripts/03_corrplot.R"
CONDA_SH = "/home/mattia/miniconda3_n/etc/profile.d/conda.sh"

# --- ENCODE BAMs -> track name (ENCODE_<mark>_<exp>) ---
ENCODE_BAMS = {}
for bam in sorted(glob.glob(REFS + "/PBMC_*.bam")):
    base = os.path.basename(bam)[:-4]                 # PBMC_H3K27ac_ENCSR105EMQ
    ENCODE_BAMS["ENCODE_" + base[len("PBMC_"):]] = bam
ENCODE_NAMES = sorted(ENCODE_BAMS)

# --- scPBMC pseudobulk bigwigs ------------------------------------------
# These used to be reused straight from the pipeline's bam_to_bw output
# (`{BASE}/{s}/{m}/bigwig/all_reads.bw`), which is built from the raw,
# unfiltered `possorted_bam.bam` -- i.e. it includes barcodes that never
# passed cell calling. That is no longer appropriate for a rigorous ENCODE
# comparison (see Xie et al. / Zhang et al., both of whom pseudobulk only
# QC-passed cells before track generation).
#
# Bigwigs are now built per (experiment, cluster) from QC-passed cells only,
# by cluster_tracks/Snakefile_cluster_tracks.smk (a.k.a. deeptools_tracks.smk),
# which sinto-splits each experiment's BAM by cell-type cluster and writes:
#   {results_dir}/{experiment}/bigwig/{cluster}.bw
# and, if that run set `merge_across_experiments: true` in its config,
# additionally:
#   {results_dir}/merged/bigwig/{cluster}.bw
#
# Confirmed layout for the PBMC runs: one cluster_tracks/deeptools_tracks.smk
# run per modality, each with its own `outdir:` under IGV_clusters/, e.g.
#   {BASE}/IGV_clusters/results_ATAC/PBMC_CT_ATAC/bigwig/{cluster}.bw
#   {BASE}/IGV_clusters/results_ATAC/PBMC_DynaTag_ATAC/bigwig/{cluster}.bw
#   {BASE}/IGV_clusters/results_H3K27ac/PBMC_CT_H3K27ac/bigwig/{cluster}.bw
#   ... etc. (one results_{modality}/ folder per modality's Snakemake run)
# so the glob has to walk that extra results_{modality}/ level; the
# {experiment}/bigwig/{cluster}.bw part underneath is unchanged.
#
# The results_{modality}/ folders are shared with other (non-PBMC) runs of
# the same cluster_tracks pipeline, so the raw glob below also picks up
# experiment folders that have nothing to do with this PBMC comparison.
# Restrict to experiments that actually belong to a known PBMC sample --
# extend PBMC_SAMPLES if you add more PBMC replicates/conditions later.
IGV_CLUSTERS = BASE + "/IGV_clusters"
INCLUDE_MERGED = False  # True -> include the merged bigwigs from cluster_tracks/deeptools_tracks.smk
PBMC_SAMPLES = ["PBMC_CT", "PBMC_DynaTag"]

SC_TRACKS = {}
skipped_non_pbmc = set()
for bw in sorted(glob.glob(IGV_CLUSTERS + "/results_*/*/bigwig/*.bw")):
    experiment = os.path.basename(os.path.dirname(os.path.dirname(bw)))  # "{experiment}" or "merged"
    cluster    = os.path.basename(bw)[:-3]                               # strip ".bw"
    if experiment == "merged":
        if not INCLUDE_MERGED:
            continue
    elif not any(experiment == s or experiment.startswith(s + "_") for s in PBMC_SAMPLES):
        skipped_non_pbmc.add(experiment)
        continue
    SC_TRACKS[f"{experiment}_{cluster}"] = bw
SC_NAMES = sorted(SC_TRACKS)

if skipped_non_pbmc:
    print("Ignored non-PBMC experiment folders under IGV_clusters:", sorted(skipped_non_pbmc))

# ordered track set for the summary (labels aligned to bigwigs)
ALL_LABELS = ENCODE_NAMES + SC_NAMES
ALL_BWS    = [f"{BW}/{n}.bw" for n in ENCODE_NAMES] + [SC_TRACKS[n] for n in SC_NAMES]

print("ENCODE bigwig jobs:", ENCODE_NAMES)
print("scPBMC pseudobulk tracks:", SC_NAMES)

rule all:
    input: OUT + "/corrplot.png"

# ---- one job per ENCODE BAM: sort/index if needed + bamCoverage (pipeline cmd) ----
rule encode_bigwig:
    output: BW + "/{name}.bw"
    params:
        bam = lambda wc: ENCODE_BAMS[wc.name],
        bin = BIN, out = BW
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 300,
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    shell:
        r"""
        mkdir -p {params.out}
        src="{params.bam}"
        # Work on a BAM inside the output dir so we never write into Refs (which
        # may be read-only): symlink if already coordinate-sorted, else sort here.
        work={params.out}/{wildcards.name}.input.bam
        if {params.bin}/samtools view -H "$src" | grep -q "SO:coordinate"; then
            ln -sf "$src" "$work"
        else
            {params.bin}/samtools sort -@ {threads} -o "$work" "$src"
        fi
        {params.bin}/samtools index -@ {threads} "$work"
        # paired-end? count paired reads on chr1 via the index (strict-mode-safe:
        # no `| head` SIGPIPE that would trip pipefail).
        paired=$({params.bin}/samtools view -c -f 1 "$work" chr1 2>/dev/null || echo 0)
        if [ "$paired" -eq 0 ]; then ext="--extendReads 200"; else ext="--extendReads"; fi
        echo "{wildcards.name}: paired reads on chr1 = $paired -> '$ext'"
        {params.bin}/bamCoverage -b "$work" -o {output} -p {threads} \
          --binSize 10 --effectiveGenomeSize 2913022398 --centerReads --smoothLength 300 --normalizeUsing RPGC --exactScaling  $ext
        """

# ---- signal matrix over genome-wide 10kb bins across all tracks ----
rule summary:
    input:
        bws = ALL_BWS
    output:
        npz = OUT + "/summary_1kb.npz",
        tab = OUT + "/counts_1kb.tab",
    params:
        bin = BIN,
        labels = " ".join(ALL_LABELS),
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    threads: 20
    resources:
        mem_mb = 32000,
        runtime = 300,
    shell:
        r"""
        {params.bin}/multiBigwigSummary bins -b {input.bws} --labels {params.labels} \
          --binSize 1000 -p {threads} -o {output.npz} --outRawCounts {output.tab}
        """

# ---- Spearman correlation -> R corrplot (png/pdf) + matrix csv ----
rule corrplot:
    input: OUT + "/counts_1kb.tab"
    output: OUT + "/corrplot.png"
    params:
        rscript = RSCRIPT,
        conda_sh = CONDA_SH,
    threads: 2
    resources:
        mem_mb = 16000,
        runtime = 120,
    shell:
        r"""
        source {params.conda_sh}
        conda run -n nanoctarna-analysis Rscript {params.rscript}
        """
