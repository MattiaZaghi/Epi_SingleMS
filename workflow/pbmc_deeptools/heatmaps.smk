# heatmaps.smk -- deeptools computeMatrix + plotHeatmap/plotProfile for the PBMC
# single-cell pseudobulk. For each PEAK SET, plot ALL THREE pseudobulk signals
# (ATAC, H3K27ac, H3K27me3) reference-point-centered over those peaks.
#
# Inputs are the per-cluster tracks produced by deeptools_tracks.smk, under
# IGV_clusters/results_<mark>/<sample>_<mark>/ :
#   signal = bigwig/<sample>.bw                             (per-cluster pseudobulk; cluster==sample for PBMC)
#   macs   = peaks/macs_frag/<sample>_<mark>_peaks.bed      (per-modality MACS2, correct narrow/broad)
# Peak sets = (mark x caller): each mark uses BOTH its MACS2 peaks (above) and its
# cellranger peaks. Applied to both samples (PBMC_CT = CUT&Tag, PBMC_DynaTag = DynaTag).
#
# Run deeptools_tracks.smk for ALL three modality configs FIRST, then:
#   cd /date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/deeptools_heatmaps
#   conda activate nanoscope_base
#   snakemake -s /home/mattia/Multi_nanoCTRNA/workflow/pbmc_deeptools/heatmaps.smk \
#       --profile htcondor --jobs 200 --cores 200
import os

BIN  = "/date/gcb/gcb_MZ/nanoCTAR_YD/MS/.snakemake/conda/67977e66f74713228faaa1470a4bdd30_/bin"  # deeptools
BASE = "/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS"
OUT  = BASE + "/deeptools_heatmaps"
IGV  = BASE + "/IGV_clusters"                    # deeptools_tracks.smk output root

SAMPLES = ["PBMC_CT", "PBMC_DynaTag"]            # both have all 3 modalities
MOD_DIR = {"ATAC": "ATAC_TATAGCCT", "H3K27ac": "H3K27ac_CCTATCCT", "H3K27me3": "H3K27me3_ATAGAGGC"}
SIGNALS = ["ATAC", "H3K27ac", "H3K27me3"]        # bigwig order in EVERY matrix

# peak sets: label -> (mark, caller). Every mark uses BOTH its MACS2 peaks
# (deeptools_tracks.smk) and its cellranger peaks.
PEAKSETS = {
    "ATAC_macs":           ("ATAC",     "macs"),
    "ATAC_cellranger":     ("ATAC",     "cellranger"),
    "H3K27ac_macs":        ("H3K27ac",  "macs"),
    "H3K27ac_cellranger":  ("H3K27ac",  "cellranger"),
    "H3K27me3_macs":       ("H3K27me3", "macs"),
    "H3K27me3_cellranger": ("H3K27me3", "cellranger"),
}

# per-signal colours (same order as SIGNALS): ATAC blue, H3K27ac red, H3K27me3 green
HEATMAP_COLORLIST = "'white,#08519c' 'white,#a50f15' 'white,#238b45'"
PROFILE_COLORS    = "'#1f78b4' '#e31a1c' '#33a02c'"


def bigwig(sample, mod):
    # per-cluster pseudobulk bigwig from deeptools_tracks.smk (cluster == sample for PBMC)
    return f"{IGV}/results_{mod}/{sample}_{mod}/bigwig/{sample}.bw"


def peak_path(sample, peakset):
    mark, caller = PEAKSETS[peakset]
    if caller == "cellranger":
        return f"{BASE}/{sample}/{MOD_DIR[mark]}/cellranger/outs/peaks.bed"
    # per-modality MACS2 peaks from deeptools_tracks.smk (correct narrow/broad)
    return f"{IGV}/results_{mark}/{sample}_{mark}/peaks/macs_frag/{sample}_{mark}_peaks.bed"


def region_label(peakset):
    mark, caller = PEAKSETS[peakset]
    return f"{mark} ({caller})"


wildcard_constraints:
    sample  = "|".join(SAMPLES),
    peakset = "|".join(PEAKSETS),


rule all:
    input:
        [f"{OUT}/heatmap_{s}__{ps}.png" for s in SAMPLES for ps in PEAKSETS]
        + [f"{OUT}/profile_{s}__{ps}.png" for s in SAMPLES for ps in PEAKSETS],


# clean peaks -> 3-col BED on main chromosomes (drops cellranger '#' header + GL scaffolds)
rule prep_bed:
    input:
        lambda wc: peak_path(wc.sample, wc.peakset),
    output:
        OUT + "/beds/{sample}__{peakset}.bed",
    shell:
        r"""
        mkdir -p $(dirname {output})
        grep '^chr' {input} | cut -f1-3 | sort -k1,1 -k2,2n > {output}
        """


# one matrix per (sample, peak set): all 3 signals, centered on peaks, +/-3kb
rule compute_matrix:
    input:
        bws=lambda wc: [bigwig(wc.sample, m) for m in SIGNALS],
        bed=OUT + "/beds/{sample}__{peakset}.bed",
    output:
        OUT + "/matrix_{sample}__{peakset}.gz",
    params:
        bin=BIN,
    threads: 16
    resources:
        mem_mb=24000,
        runtime=240,
    shell:
        r"""
        {params.bin}/computeMatrix reference-point --referencePoint center \
          -S {input.bws} -R {input.bed} -a 3000 -b 3000 --binSize 50 \
          --missingDataAsZero --averageTypeBins mean -p {threads} -o {output}
        """


rule plot_heatmap:
    input:
        OUT + "/matrix_{sample}__{peakset}.gz",
    output:
        OUT + "/heatmap_{sample}__{peakset}.png",
    params:
        bin=BIN,
        colors=HEATMAP_COLORLIST,
        labels=" ".join(SIGNALS),
        region=lambda wc: region_label(wc.peakset),
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        r"""
        {params.bin}/plotHeatmap -m {input} -o {output} --dpi 300 \
          --colorList {params.colors} --missingDataColor white \
          --sortRegions descend --sortUsing mean \
          --refPointLabel center --samplesLabel {params.labels} \
          --regionsLabel "{params.region} peaks" \
          --whatToShow 'plot, heatmap and colorbar' \
          --heatmapHeight 12 --heatmapWidth 4
        """


rule plot_profile:
    input:
        OUT + "/matrix_{sample}__{peakset}.gz",
    output:
        OUT + "/profile_{sample}__{peakset}.png",
    params:
        bin=BIN,
        colors=PROFILE_COLORS,
        labels=" ".join(SIGNALS),
        region=lambda wc: region_label(wc.peakset),
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        r"""
        {params.bin}/plotProfile -m {input} -o {output} --dpi 300 \
          --perGroup --colors {params.colors} --plotWidth 8 \
          --refPointLabel center --samplesLabel {params.labels} \
          --regionsLabel "{params.region} peaks"
        """
