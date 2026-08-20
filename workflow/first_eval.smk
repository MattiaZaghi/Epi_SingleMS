# first_eval.smk
# ---------------------------------------------------------------------------
# Render the streamlined first-evaluation HTML report as a single job.
# Runs after cell-calling (depends on every experiment's cell_picking/metadata.csv)
# and, under `--profile htcondor`, is submitted to Condor automatically.
#
# Run from the SAME working directory as the preprocessing run (the folder that
# holds the {sample}/{modality}_{barcode}/ outputs), e.g.:
#   cd /date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction
#   snakemake -s /date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/workflow/first_eval.smk \
#       --configfile /date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/config/config_human_MGI.yaml \
#       --profile htcondor --jobs 1
#
# The report itself DISCOVERS experiments from the working directory, so it will
# also pick up samples present on disk but absent from the config; the config is
# only used here to gate the job on the expected metadata files.
# ---------------------------------------------------------------------------
include: 'Snakefile_prep.smk'

REPODIR   = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction"
ANALYSIS_ENV = "nanoctarna-analysis"        # conda env built from envs/nanoctarna_analysis.yaml
CONDA_SH  = "/home/mattia/miniconda3_n/etc/profile.d/conda.sh"

rule all_first_eval:
    input:
        "first_eval/first_eval_report.html"

rule render_first_eval:
    input:
        # Gate on the cell-calling metadata that ALREADY EXISTS on disk (this
        # standalone Snakefile cannot build it -- that happens in the preprocess
        # DAG). Listing only existing files means: re-render when they change, run
        # even while preprocessing is still finishing (the report skips whatever
        # is not yet present). Run this after cell-calling for a complete report.
        metadata = lambda wildcards: [
            f for f in (
                '{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(
                    sample=sample, modality=modality, barcode=barcodes_dict[sample][modality])
                for sample in samples_list for modality in barcodes_dict[sample].keys())
            if os.path.exists(f)],
        # the analysis code: editing any of these re-triggers the render
        code = [REPODIR + "/analysis/main/first_eval_MGI.Rmd",
                REPODIR + "/analysis/main/render_first_eval.R",
                REPODIR + "/scripts/functions_firstEval.R",
                REPODIR + "/scripts/functions_scCT.R"],
    output:
        html = "first_eval/first_eval_report.html",
    params:
        render  = REPODIR + "/analysis/main/render_first_eval.R",
        wd      = os.getcwd(),
        fig_dir = os.getcwd() + "/first_eval/figures",
        env     = ANALYSIS_ENV,
        conda_sh = CONDA_SH,
        include_rna = config['general'].get('include_rna', False),
        group_pattern = config['general'].get('group_pattern', '^(MS)'),
        reference_yaml = config['general'].get('reference_yaml', ''),   # set to enable QC benchmarking against another saved .rds
        dataset_label = config['general'].get('dataset_label', 'Human-MS'),
        compute_tss = config['general'].get('compute_tss', False),
        reuse_counts = config['general'].get('reuse_counts', True),
    threads: 12
    resources:
        mem_mb  = 160000,   # bin-matrix building over many experiments x parallel workers
        runtime = 360,      # minutes
    shell:
        "source {params.conda_sh}; "
        "conda run -n {params.env} Rscript {params.render} "
        "workers={threads} group_pattern='{params.group_pattern}' "
        "wd={params.wd}/ experiment_source=discover "
        "include_rna={params.include_rna} "
        "reference_yaml='{params.reference_yaml}' dataset_label='{params.dataset_label}' "
        "compute_tss={params.compute_tss} reuse_counts={params.reuse_counts} "
        "fig_dir={params.fig_dir} out_html={params.wd}/{output.html} 2>&1"
