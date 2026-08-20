# Epi_SingleMS

Bioinformatics pipeline and downstream analysis for **single-nucleus/single-cell nano-CUT&Tag + RNA data from human Multiple Sclerosis (MS) brain and spinal cord samples**.

This repo is the MS-specific branch of the more general `Multi_nanoCTRNA` pipeline: it keeps the shared preprocessing/analysis infrastructure but only the code and configs relevant to the human MS brain project (no mouse/embryo development data, no reanalysis of unrelated published datasets, no other-tissue comparisons).

## Pipeline run location

The pipeline is run from `/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/` (a full checkout of this repo alongside its data). Raw per-sample MGI fastqs live under its `raw_data/{sample}/` subfolder, and preprocessing outputs are written as `{sample}/{modality}_{barcode}/` directly under that same root. Several scripts/configs in this repo default to that path — repoint them if you run elsewhere.

## What this repo does

Two tiers, same as the parent pipeline:

1. **Snakemake preprocessing** (`workflow/`, `config/`) — raw FASTQ (MGI) → demultiplexed reads per modality → alignment/counting via Cell Ranger → per-cell barcode QC metrics → GMM-based cell calling → filtered fragment/count matrices.
2. **Downstream R analysis** (`analysis/`, `qc/`, `scripts/`) — builds Seurat/Signac objects per modality and sample, merges/integrates (WNN or Harmony), clusters, annotates against reference brain atlases, and produces QC plots — specifically for the MS brain cohort.

Data originates from the `nanoscope`-derived nano-CUT&Tag protocol, sequenced on the MGI platform (R1+R2 only, modality barcode packed into R2).

## Samples covered

Per `config/config_human_MGI.yaml` (epigenomic fraction only — H3K27ac + H3K27me3):

| Library | Description |
|---|---|
| `MS3BL1_CT`, `MS3BL2_CT` | Pooled MS brain samples (patients MS381, MS058, MS549) |
| `MS549_CT` | MS549 brain sample |

## Repository layout

```
config/            Snakemake config for the human/MGI preprocessing run
workflow/           Snakemake rules and helper scripts
  Snakefile_preprocess.smk    Single-cell nano-CUT&Tag (epigenomic) preprocessing
  Snakefile_preprocess_RNA.smk  RNA-only preprocessing
  Snakefile_prep*.smk         FASTQ-parsing includes shared by the preprocess Snakefiles
  first_eval.smk / .sub       Renders the parametrized "first evaluation" QC report
  scripts/                    debarcode.py, pick_cells.R, matrix/GTF filtering, UMAP plotting helpers
analysis/
  main/                       MS brain Rmd/notebook analyses (WNN + bins integration, first-evaluation report)
qc/                            MS QC comparison plotting (qc/QC_nanoCTAR_hMS.R)
scripts/                        Shared R helper libraries (functions_scCT.R, GetCells.R, functions_firstEval.R)
standalone/                     One-off utility scripts (genome binning, library-complexity curves, h5->Seurat conversion)
envs/                            Conda environment specs for the pipeline and the R analysis layer
10x_barcodes/                   10x Genomics ATAC/RNA barcode whitelists
```

## Running the preprocessing pipeline

The `configfile:` directive is commented out in every Snakefile — pass `--configfile` explicitly.

```bash
# Single-cell nano-CUT&Tag (epigenomic fraction), human/MGI
snakemake -s workflow/Snakefile_preprocess.smk --configfile config/config_human_MGI.yaml --cores N --use-conda

# RNA only
snakemake -s workflow/Snakefile_preprocess_RNA.smk --configfile config/config_human_MGI.yaml --cores N --use-conda

# First-evaluation QC report (after cell calling)
snakemake -s workflow/first_eval.smk --configfile config/config_human_MGI.yaml --profile htcondor --jobs 1
```

Preprocess DAG (single-cell): `demultiplex` → `run_cellranger_ATAC`/`run_cellranger_RNA` → `bam_to_bw` (QC bigwig) → `run_macs_broad` (peaks) → `barcode_metrics_peaks`/`barcode_metrics_all` → `cell_selection` (`pick_cells.R`) → `clean_cellranger_output` → `bam_to_fragments` (sinto) → `sort_sinto_output` → `preseq` → `get_cells` → `create_matrix_{peaks,bins,genes}` (fragtk).

Each `Snakefile_preprocess*.smk` relies on a `ruleorder:` chain rather than a single `all` target — many outputs in `rule all_preprocess` are commented out, so uncomment the target you want to build.

## Downstream analysis

- `analysis/main/human_MS.rmd` / `human_MS.ipynb` — bins-based downstream analysis of the MS brain nano-CUT&Tag+RNA replicates (load → QC → merge → normalize/DR → cluster → annotate).
- `analysis/main/MS_WNN.rmd` — the multimodal **WNN** integration (RNA + H3K27ac + H3K27me3 via `FindMultiModalNeighbors`).
- `analysis/main/first_eval_MGI.Rmd` (+ `render_first_eval.R`) — parametrized, sample-agnostic "first evaluation" QC report; discovers experiments from a working directory rather than hardcoding sample names.
- `analysis/main/first_eval_interactive.rmd` — interactive companion to the report above, for chunk-by-chunk exploration in RStudio (doublet removal, Harmony batch integration, label transfer against reference brain atlases).
- `qc/QC_nanoCTAR_hMS.R` — cross-sample QC comparison plotting for the MS brain dataset.
- `scripts/functions_scCT.R` — shared plotting/QC helper library `source()`d by the Rmds/notebooks.
- `scripts/functions_firstEval.R` — pure-function library backing the first-evaluation report (sample/modality-agnostic, grouping defaults to the `MS` sample prefix).
- `scripts/GetCells.R` — manual/interactive cell gating gadget (Shiny + `gatepoints`).

## Environments

- **Analysis (R)**: `envs/nanoctarna_analysis.yaml` → conda env `nanoctarna-analysis` (R ≥4.1 + Seurat/Signac/rmarkdown stack).
- **Pipeline**: `envs/nanoscope_general.yaml` (samtools, bedtools, deeptools, macs2, sinto, pysam, R+Rsamtools), `envs/nanoscope_debarcode.yaml` (for `debarcode.py`), `envs/nanoscope_base*.yaml`, `envs/nanoscope_amulet.yaml` (AMULET doublet detection).

### External prerequisites not in repo envs

- **cellranger-atac**, **cellranger** (GEX) — via `config['general']['cellranger_software*']`.
- **fragtk** — builds the count matrices.
- **preseq** (`lc_extrap`) — library-complexity estimation.
- **AMULET** — doublet detection, used in `first_eval_interactive.rmd`.
