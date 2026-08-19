import os
import re

# ---------------------------------------------------------------------------
#  Per-cluster pseudobulk tracks from a merged nanoCTAR object.
#
#  For each experiment (one cellranger BAM) sinto splits the BAM into one BAM
#  per predicted.id_Human_MS cluster, then deepTools makes a bigwig per cluster.
#  Splitting is per-experiment because the same barcode recurs across experiments.
#
#  Run:
#    snakemake -s cluster_tracks/Snakefile_cluster_tracks.smk \
#              --configfile cluster_tracks/config.yaml --cores 16 --use-conda
# ---------------------------------------------------------------------------

OUTDIR      = config["outdir"]
BCDIR       = config.get("barcodes_dir", "barcodes")
BAMS        = config["bams"]                       # {experiment: /path/possorted_bam.bam}
EXPERIMENTS = list(BAMS)
BARCODETAG  = config.get("barcodetag", "CB")
MERGE       = config.get("merge_across_experiments", False)


def read_clusters(tsv):
    cl = set()
    with open(tsv) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1]:
                cl.add(p[1])
    return sorted(cl)


CLUSTERS = {e: read_clusters(os.path.join(BCDIR, f"{e}.tsv")) for e in EXPERIMENTS}

# cluster -> experiments that contain it (for the optional merge stage)
CLUSTER2EXP = {}
for e in EXPERIMENTS:
    for c in CLUSTERS[e]:
        CLUSTER2EXP.setdefault(c, []).append(e)

# constrain {experiment} to the real names so it can never match the literal
# "merged" path component (which would make the two bamcoverage rules ambiguous)
wildcard_constraints:
    experiment="|".join(re.escape(e) for e in EXPERIMENTS),
    cluster="[^/]+",


def targets():
    t = []
    # NB: do NOT list the per-cluster <experiment>/bam/<cluster>.bam here. They
    # are a SIDE EFFECT of the single sinto_split pass, whose only declared
    # output is the .split.done flag -- no rule has <cluster>.bam as an output,
    # so requesting them directly raises MissingInputException. The per-cluster
    # .bw targets below already force sinto_split (bamcoverage depends on the
    # flag), which is what creates every cluster BAM.
    t += [os.path.join(OUTDIR, e, "bigwig", f"{c}.bw")
          for e in EXPERIMENTS for c in CLUSTERS[e]]
    t += [os.path.join(OUTDIR, e, "fragments", "fragments.tsv.gz")
          for e in EXPERIMENTS]
    # one peakset per experiment (MACS runs on the whole experiment's pooled
    # fragments). NB: iterate only over experiments -- do NOT reuse `c` here, it
    # would leak the last cluster from the CLUSTER2EXP loop above and misname the
    # files (e.g. a PBMC_DynaTag peak file landing in the PBMC_CT directory).
    #t += [os.path.join(OUTDIR, e, "peaks", "macs_frag", f"{e}_peaks.bed")
          #for e in EXPERIMENTS]
    if MERGE:
        t += [os.path.join(OUTDIR, "merged", "bigwig", f"{c}.bw") for c in CLUSTER2EXP]
    return t
# ============================================================================
#  Per-modality MACS2 settings
#  narrow  -> ATAC, H3K4me3, H3K27ac   (punctate)
#  broad   -> H3K27me3, H3K9me3, ...    (domains)
# ============================================================================
MACS_FLAGS = {
    "ATAC":     "-q 0.01 --call-summits",
    "H3K4me3":  "-q 0.01 --call-summits",
    "H3K27ac":  "-q 0.01 --call-summits",
    "H3K27me3": "--broad --broad-cutoff 0.05 --max-gap 1000 --llocal 100000 -q 0.01",
}
MACS_FLAGS_DEFAULT = "--broad --broad-cutoff 0.05 --max-gap 1000 --llocal 100000 -q 0.01"
BROAD_MODALITIES   = {"H3K27me3", "H3K9me3", "H3K36me3"}

# The {experiment} wildcard is the BAM/experiment NAME (e.g. "PBMC_CT_H3K27me3"),
# not the modality -- so keying MACS_FLAGS on it never matched and everything fell
# through to the broad default. Derive the mark from the name (it is embedded as a
# token, e.g. ..._H3K27me3 / ..._ATAC) so each modality gets the right flags.
KNOWN_MARKS = set(MACS_FLAGS) | BROAD_MODALITIES | {"ATAC", "H3K4me1", "H3K9ac"}

def modality_of(experiment):
    for tok in experiment.split("_"):           # exact token first
        if tok in KNOWN_MARKS:
            return tok
    for m in sorted(KNOWN_MARKS, key=len, reverse=True):   # substring fallback
        if m in experiment:
            return m
    return None

def macs_flags(wildcards):
    return MACS_FLAGS.get(modality_of(wildcards.experiment), MACS_FLAGS_DEFAULT)

def macs_peak_src(wildcards):                 # extension MACS2 actually writes
    return "broadPeak" if modality_of(wildcards.experiment) in BROAD_MODALITIES else "narrowPeak"

rule all:
    input:
        targets()


# ---- 1) split one experiment's BAM into per-cluster BAMs (single sinto pass)
rule sinto_split:
    input:
        bam=lambda wc: BAMS[wc.experiment],
        cells=lambda wc: os.path.join(BCDIR, f"{wc.experiment}.tsv"),
    output:
        flag=os.path.join(OUTDIR, "{experiment}", "bam", ".split.done"),
    params:
        bamdir=os.path.join(OUTDIR, "{experiment}", "bam"),
        tag=BARCODETAG,
    threads: config.get("threads", 8)
    conda: "/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml"
    shell:
        r"""
        mkdir -p {params.bamdir}
        sinto filterbarcodes -b {input.bam} -c {input.cells} \
            --barcodetag {params.tag} --outdir {params.bamdir} -p {threads}
        for b in {params.bamdir}/*.bam; do
            samtools sort -@ {threads} -o "$b".tmp "$b" && mv "$b".tmp "$b"
            samtools index "$b"
        done
        touch {output.flag}
        """
rule bam_to_fragments:
    input:
        flag=os.path.join(OUTDIR, "{experiment}", "bam", ".split.done"),
        bam=lambda wc: BAMS[wc.experiment],
    output:
        raw = os.path.join(OUTDIR, "{experiment}", "fragments", "fragments_raw.tsv"),
    conda: "/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml"
    threads: 20
    resources:
        mem_mb=16000
    shell:
        'sinto fragments -b {input.bam} -f {output.raw} -p {threads}'

rule sort_sinto_output:
    input:
        raw = os.path.join(OUTDIR, "{experiment}", "fragments", "fragments_raw.tsv"),
    output:
        fragments = os.path.join(OUTDIR, "{experiment}", "fragments", "fragments.tsv.gz"),
        index      = os.path.join(OUTDIR, "{experiment}", "fragments", "fragments.tsv.gz.tbi"),
    conda: "/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml"
    resources:
        mem_mb = 16000
    shell:
        'sort -k1,1 -k2,2n {input.raw} | bgzip > {output.fragments}; '
        'tabix -p bed {output.fragments} '

rule run_macs_fragments:
    input:
        fragments = os.path.join(OUTDIR, "{experiment}", "fragments", "fragments.tsv.gz"),
    output:
        peaks=os.path.join(OUTDIR, "{experiment}", "peaks", "macs_frag", "{experiment}_peaks.bed"),
        bedpe=temp(os.path.join(OUTDIR, "{experiment}", "peaks", "macs_frag", "{experiment}_fragments.bedpe")),
    params:
        macs_outdir=os.path.join(OUTDIR, "{experiment}", "peaks", "macs_frag"),
        macs_genome=config.get("macs_genome", "hs"),
        flags=macs_flags,                      # narrow/broad flags picked per modality
    conda: "/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml"
    resources:
        mem_mb = 16000
    shell:
        r'''
        mkdir -p {params.macs_outdir}

        zcat {input.fragments} | grep -v '^#' | cut -f1-3 > {output.bedpe}

        macs2 callpeak -t {output.bedpe} -g {params.macs_genome} -f BEDPE \
            -n {wildcards.experiment} --outdir {params.macs_outdir} --keep-dup all {params.flags}

        # MACS2 writes *_peaks.narrowPeak OR *_peaks.broadPeak depending on --broad;
        # normalise to one canonical .bed without guessing which mode ran.
        if [ -f {params.macs_outdir}/{wildcards.experiment}_peaks.narrowPeak ]; then
            cp {params.macs_outdir}/{wildcards.experiment}_peaks.narrowPeak {output.peaks}
        elif [ -f {params.macs_outdir}/{wildcards.experiment}_peaks.broadPeak ]; then
            cp {params.macs_outdir}/{wildcards.experiment}_peaks.broadPeak {output.peaks}
        else
            echo "ERROR: no MACS2 peak file for {wildcards.experiment}" >&2
            exit 1
        fi
        '''
# ---- 2) one bigwig per (experiment, cluster) -------------------------------
rule bamcoverage:
    input:
        flag=os.path.join(OUTDIR, "{experiment}", "bam", ".split.done"),
    output:
        bw=os.path.join(OUTDIR, "{experiment}", "bigwig", "{cluster}.bw"),
    params:
        bam=os.path.join(OUTDIR, "{experiment}", "bam", "{cluster}.bam"),
        binsize=config.get("binsize", 10),
        norm=config.get("normalize", "RPKM"),
        gsize=config.get("effective_genome_size", 2913022398),
        blacklist=(("--blackListFileName " + config["blacklist"])
                   if config.get("blacklist") else ""),
    threads: config.get("threads", 20)
    conda: 
        "/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml"
    shell:
        r"""
        bamCoverage -b {params.bam} -o {output.bw} -p {threads} --binSize {params.binsize} --effectiveGenomeSize {params.gsize} --normalizeUsing {params.norm} --extendReads --centerReads --smoothLength 300 --exactScaling  
        """


# ---- 3) OPTIONAL: pool the same cluster across experiments -----------------
#  Enable with `merge_across_experiments: true` in the config. To pool by
#  condition/patient instead of everything, edit CLUSTER2EXP to group only the
#  experiments you want (e.g. filter by a prefix of the experiment name).
if MERGE:

    rule merge_cluster:
        input:
            flags=lambda wc: [os.path.join(OUTDIR, e, "bam", ".split.done")
                              for e in CLUSTER2EXP[wc.cluster]],
        output:
            bam=os.path.join(OUTDIR, "merged", "bam", "{cluster}.bam"),
        params:
            bams=lambda wc: " ".join(os.path.join(OUTDIR, e, "bam", f"{wc.cluster}.bam")
                                     for e in CLUSTER2EXP[wc.cluster]),
        threads: config.get("threads", 20)
        conda: "/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output.bam})
            samtools merge -f -@ {threads} {output.bam} {params.bams}
            samtools index {output.bam}
            """

    rule bamcoverage_merged:
        input:
            bam=os.path.join(OUTDIR, "merged", "bam", "{cluster}.bam"),
        output:
            bw=os.path.join(OUTDIR, "merged", "bigwig", "{cluster}.bw"),
        params:
            binsize=config.get("binsize", 50),
            norm=config.get("normalize", "RPKM"),
            gsize=config.get("effective_genome_size", 2913022398),
            blacklist=(("--blackListFileName " + config["blacklist"])
                       if config.get("blacklist") else ""),
        threads: config.get("threads", 8)
        conda: "/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml"
        shell:
            r"""
            bamCoverage -b {input.bam} -o {output.bw} -p {threads} \
                --binSize {params.binsize} --normalizeUsing {params.norm} --effectiveGenomeSize {params.gsize}  \
                --extendReads --centerReads --smoothLength 300 
                
            """