include: 'Snakefile_prep.smk'

rule all_preprocess:
    input:
        cellranger=[
            '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        bigwig_all=[
            '{sample}/{modality}_{barcode}/bigwig/all_reads.bw'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        macs_broad=[
            '{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        macs_macs_per_modality=[
            '{sample}/{modality}_{barcode}/peaks/macs/{modality}_peaks.bed'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        macs_fragments=[
            '{sample}/{modality}_{barcode}/peaks/macs_frag/{modality}_peaks.bed'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        peaks_overlap=[
            '{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        barcodes_sum=[
            '{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        cell_pick=[
            '{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        multiqc = 'multiqc/multiqc_report.html',      # single aggregated FastQC report over all post-debarcoding fastqs (all samples)
        #cellranger_cleanup = [
            #'{sample}/{modality}_{barcode}/_clean_cellranger.out'.format(sample=sample,modality=modality,barcode=
            #barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        sort_sinto_output = [
            '{sample}/{modality}_{barcode}/cellranger/outs/fragments_raw.tsv.gz'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        preseq = [
            '{sample}/{modality}_{barcode}/cellranger/outs/yield.txt'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
       # matrix_peaks = [
            #'{sample}/{modality}_{barcode}/matrix/matrix_peaks/'.format(sample=sample,modality=modality,barcode=
            #barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
       # matrix_bins = [
            #'{sample}/{modality}_{barcode}/matrix/matrix_bin_{binsize}/'.format(sample=sample,modality=modality,barcode=
            #barcodes_dict[sample][modality],binsize = binsize) for sample in samples_list for modality in barcodes_dict[sample].keys() for binsize in bins], 
        #matrix_genes = [
            #'{sample}/{modality}_{barcode}/matrix/matrix_genes/'.format(sample=sample,modality=modality,barcode=
            #barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        #bam_RNA=['{sample}/RNA_AAAAGGGG/cellranger/outs/possorted_genome_bam.bam'.format(sample=sample) for sample in samples_list],
        #agg=['{merge}/{modality_n}_{barcode_n}/cellranger/outs/fragments.tsv.gz'.format(merge=merge,modality_n=modality,barcode_n=barcodes_gen[merge][modality]) for merge in params_list for modality in barcodes_gen[merge].keys()],
        #bw_merge_ATAC =['{merge}/ATAC_TATAGCCT/bigwig/all_reads.bw'.format(merge=merge) for merge in params_list],
        #bw_merge_H3K27ac =['{merge}/H3K27ac_CCTATCCT/bigwig/all_reads.bw'.format(merge=merge) for merge in params_list],
        #bw_merge_H3K27me3 =['{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads.bw'.format(merge=merge) for merge in params_list],
        #bam_merge_ATAC =['{merge}/ATAC_TATAGCCT/cellranger/outs/possorted_bam.bam'.format(merge=merge) for merge in params_list],
        #bam_merge_H3K27ac =['{merge}/H3K27ac_CCTATCCT/cellranger/outs/possorted_bam.bam'.format(merge=merge) for merge in params_list],
        #bam_merge_H3K27me3 =['{merge}/H3K27me3_ATAGAGGC/cellranger/outs/possorted_bam.bam'.format(merge=merge) for merge in params_list]


ruleorder: demultiplex > run_cellranger_ATAC >  bam_to_bw > run_macs_broad > barcode_metrics_peaks > barcode_metrics_all > cell_selection > clean_cellranger_output > bam_to_fragments >  sort_sinto_output > preseq  #get_cells >  #> #create_matrix_peaks > create_matrix_bins > create_genebody_and_promoter_matrix 
#merge_bw_ATAC > merge_bw_H3K27ac > merge_bw_H3K27me3 > merge_bam_ATAC > merge_bam_H3K27ac > merge_bam_H3K27me3 run_cellranger_RNA >


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

def macs_flags(wildcards):
    return MACS_FLAGS.get(wildcards.modality, MACS_FLAGS_DEFAULT)

def macs_peak_src(wildcards):                 # extension MACS2 actually writes
    return "broadPeak" if wildcards.modality in BROAD_MODALITIES else "narrowPeak"

# V4 design: debarcode.py reads the 8-bp modality barcode from R3 (sequencing R2),
# flanked by two mosaic ends (ME=AGATGTGTATAAGAGACAG), trims the cassette off R3, and
# keeps the 16-bp cell barcode on R2 (i5) untouched. ME/barcode length/mismatch use
# debarcode.py defaults; pass --me/--barcode_length/--me_mismatch to override.
rule demultiplex:
    input:
        script=workflow.basedir + '/scripts/debarcode.py',
        fastq=lambda wildcards: glob.glob(config['samples'][wildcards.sample]['fastq_path'] + '/**/*{lane}*R[{reads}]*.fastq.gz'.format(
            lane=wildcards.lane,
            reads='12' if config['samples'][wildcards.sample].get('platform', 'illumina') == 'mgi' else '123'),recursive=True)
    output:
        '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R1_{suffix}',
        '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R2_{suffix}',
        '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R3_{suffix}',
    params:
        nbarcodes=lambda wildcards: len(config['samples'][wildcards.sample]['barcodes']),
        out_folder=lambda wildcards: '{sample}/{modality}_{barcode}/fastq/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode),
        platform=lambda wildcards: config['samples'][wildcards.sample].get('platform', 'illumina'),
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 720 # 8 hours should be enough for most reasonable single-cell data, if the command timeouts, increase this value
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    shell:
        "python3 {input.script} -i {input.fastq} -o {params.out_folder} --mismatch 2 --single_cell --platform {params.platform} --barcode {wildcards.barcode} --name {wildcards.sample} 2>&1"

rule run_cellranger_ATAC:
    input:
        lambda wildcards: get_fastq_for_cellranger(config['samples'][wildcards.sample]['fastq_path'],sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    output:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        bai='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam.bai',
        frag='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        peaks='{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    params:
        cellranger_software=config['general']['cellranger_software'],
        cellranger_ref=config['general']['cellranger_ref'],
        fastq_folder=lambda wildcards: os.getcwd() + '/{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    threads: 20
    resources:
        mem_mb = 32000
    shell:
        'rm -rf {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/cellranger/; '
        'cd {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/; '
        '{params.cellranger_software} count --id cellranger --chemistry=ARC-v1 --reference {params.cellranger_ref} --fastqs {params.fastq_folder}'

rule run_cellranger_RNA:
    input:
        lambda wildcards: get_fastq_for_cellranger_rna(config['samples'][wildcards.sample]['fastq_path_RNA']+ '/**/*{lane}*R[12]*.fastq.gz',sample=wildcards.sample)
    output:
        bam_RNA='{sample}/RNA_AAAAGGGG/cellranger/outs/possorted_genome_bam.bam'
    params:
        cellranger_software=config['general']['cellranger_software_RNA'],
        cellranger_ref=config['general']['cellranger_ref_RNA'],
        fastq_folder=lambda wildcards: config['samples'][wildcards.sample]['fastq_path_RNA']#+'/02-FASTQ/20231116_LH00217_0027_A22FMT5LT3/'
    threads: 20
    resources:
        mem_mb = 32000
    shell:
        'rm -rf {wildcards.sample}/RNA_AAAAGGGG/cellranger/; '
        'cd {wildcards.sample}/RNA_AAAAGGGG/; '
        '{params.cellranger_software} count --id cellranger --transcriptome {params.cellranger_ref} --chemistry=ARC-v1 --fastqs {params.fastq_folder} --create-bam=true'

rule bam_to_bw: # For QC reasons
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        bigwig='{sample}/{modality}_{barcode}/bigwig/all_reads.bw'
    threads: 20
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'bamCoverage -b {input.cellranger_bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
        ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

rule run_macs:
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        peaks='{sample}/{modality}_{barcode}/peaks/macs/{modality}_peaks.bed'
    params:
        macs_outdir='{sample}/{modality}_{barcode}/peaks/macs',
        macs_genome=config['general']['macs_genome'],
        flags=macs_flags,                      # narrow/broad flags per modality (as before)
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        r'''
        mkdir -p {params.macs_outdir}
        macs2 callpeak -t {input.cellranger_bam} -g {params.macs_genome} -f BAMPE \
            -n {wildcards.modality} --outdir {params.macs_outdir} --keep-dup 1 {params.flags}

        # MACS2 writes *_peaks.narrowPeak OR *_peaks.broadPeak depending on --broad;
        # normalise to one canonical output without guessing which.
        if [ -f {params.macs_outdir}/{wildcards.modality}_peaks.narrowPeak ]; then
            cp {params.macs_outdir}/{wildcards.modality}_peaks.narrowPeak {output.peaks}
        elif [ -f {params.macs_outdir}/{wildcards.modality}_peaks.broadPeak ]; then
            cp {params.macs_outdir}/{wildcards.modality}_peaks.broadPeak {output.peaks}
        else
            echo "ERROR: no MACS2 peak file for {wildcards.modality}" >&2
            exit 1
        fi
        '''

rule run_macs_broad:
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        broad_peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'
    params:
        macs_outdir='{sample}/{modality}_{barcode}/peaks/macs_broad/',
        macs_genome=config['general']['macs_genome']
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
        '--max-gap 1000 --broad 2>&1 '

rule barcode_metrics_peaks:
    input:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
    output:
        overlap='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'
    params:
        get_cell_barcode=workflow.basedir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list=workflow.basedir + '/scripts/add_sample_to_list.py',
        tmpdir=config['general']['tempdir']
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -f2 | '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" |  '
        'sort -T {params.tmpdir} | uniq -c > {output.overlap} && [[ -s {output.overlap} ]] ; '

# ----------------------------------------------------------------------------
#  Pseudobulk peak calling straight from the cellranger-atac fragments file,
#  instead of the possorted BAM (see run_macs / run_macs_broad above).
#
#  cellranger-atac's fragments.tsv.gz is already Tn5-offset corrected and
#  PCR-deduplicated -- one line per unique fragment, columns 1-3 are exactly
#  chrom/start/end of the whole insert. That is precisely what MACS2's "-f
#  BEDPE" mode expects (its own private format: chrom, fragment start, fragment
#  end -- NOT bedtools' 6-column BEDPE), so no BAM->tagAlign shifting/modeling
#  is needed: just drop the barcode/dup-count columns and hand MACS2 the rest.
#  Because the fragments are already deduplicated, use --keep-dup all (the
#  BAM-based rules above use --keep-dup 1, which matters when raw reads/dups
#  are still in play).
# ----------------------------------------------------------------------------
rule run_macs_fragments:
    input:
        fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'
    output:
        peaks='{sample}/{modality}_{barcode}/peaks/macs_frag/{modality}_peaks.bed',
        bedpe=temp('{sample}/{modality}_{barcode}/peaks/macs_frag/{modality}_fragments.bedpe'),
    params:
        macs_outdir='{sample}/{modality}_{barcode}/peaks/macs_frag',
        macs_genome=config['general']['macs_genome'],
        flags=macs_flags,                      # same narrow/broad-per-modality flags as run_macs
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        r'''
        mkdir -p {params.macs_outdir}

        zcat {input.fragments} | grep -v '^#' | cut -f1-3 > {output.bedpe}

        macs2 callpeak -t {output.bedpe} -g {params.macs_genome} -f BEDPE \
            -n {wildcards.modality} --outdir {params.macs_outdir} --keep-dup all {params.flags}

        # MACS2 writes *_peaks.narrowPeak OR *_peaks.broadPeak depending on --broad;
        # normalise to one canonical output without guessing which (same convention
        # as run_macs above).
        if [ -f {params.macs_outdir}/{wildcards.modality}_peaks.narrowPeak ]; then
            cp {params.macs_outdir}/{wildcards.modality}_peaks.narrowPeak {output.peaks}
        elif [ -f {params.macs_outdir}/{wildcards.modality}_peaks.broadPeak ]; then
            cp {params.macs_outdir}/{wildcards.modality}_peaks.broadPeak {output.peaks}
        else
            echo "ERROR: no MACS2 peak file for {wildcards.modality}" >&2
            exit 1
        fi
        '''

rule barcode_metrics_all:
    input:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
    output:
        all_bcd='{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'
    params:
        get_cell_barcode=workflow.basedir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list=workflow.basedir + '/scripts/add_sample_to_list.py',
        tmpdir=config['general']['tempdir']
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'mkdir -p {params.tmpdir}; '
        ' samtools view -f2 {input.bam}| '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" |  '
        'sort -T {params.tmpdir} | uniq -c > {output.all_bcd} && [[ -s {output.all_bcd} ]] ; '

rule cell_selection:
    input:
        bcd_all='{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt',
        bcd_peak='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
        metadata='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
    output:
        '{sample}/{modality}_{barcode}/cell_picking/cells_10x.png',
        '{sample}/{modality}_{barcode}/cell_picking/cells_picked.png',
        '{sample}/{modality}_{barcode}/cell_picking/metadata.csv',
    params:
        script=workflow.basedir + '/scripts/pick_cells.R',
        out_prefix='{sample}/{modality}_{barcode}/cell_picking/',
    resources:
        mem_mb = 25000
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    shell:
        "Rscript {params.script} --metadata {input.metadata} --fragments {input.fragments} --bcd_all {input.bcd_all} --bcd_peak {input.bcd_peak} --modality {wildcards.modality} --sample {wildcards.sample} --out_prefix {params.out_prefix}"

rule clean_cellranger_output:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta = '{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        peaks = '{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    output:
        '{sample}/{modality}_{barcode}/_clean_cellranger.out'
    params:
        cellranger_folder = '{sample}/{modality}_{barcode}/cellranger/'
    shell:
        'touch {params.cellranger_folder}/tmp.txt;'                     # Create temp empty file to avoid error if the directory is empty
        'ls -d  {params.cellranger_folder}/* | grep -v outs | xargs rm -r; '
        'touch {output}'

rule bam_to_fragments:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        index= '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam.bai',
        frag = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'
    output:
        fragmments = temp('{sample}/{modality}_{barcode}/cellranger/outs/fragments_raw.tsv'),
    conda: '../envs/nanoscope_general.yaml'
    threads: 20
    resources:
        mem_mb=16000
    shell:
        'sinto fragments -b {input.bam} -f {output.fragmments} -p {threads}'

rule sort_sinto_output:
    input:
        fragments = '{sample}/{modality}_{barcode}/cellranger/outs/fragments_raw.tsv',
    output:
        fragments = '{sample}/{modality}_{barcode}/cellranger/outs/fragments_raw.tsv.gz',
        index      = '{sample}/{modality}_{barcode}/cellranger/outs/fragments_raw.tsv.gz.tbi',
    conda: '../envs/nanoscope_general.yaml'
    resources:
        mem_mb = 16000
    shell:
        'sort -k1,1 -k2,2n {input.fragments} | bgzip > {output.fragments}; '
        'tabix -p bed {output.fragments} '

rule preseq:
    input:
        bam = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'  # Fixed input path
    output:
        preseq = '{sample}/{modality}_{barcode}/cellranger/outs/yield.txt'
    resources:
        mem_mb = 16000
    conda:
        "/home/mattia/miniconda3_n/envs/preseq.yml"
    shell:
        """

        preseq lc_extrap -B -P -e 2000000000 -o {output.preseq} {input.bam} 
        """

# FastQC on the debarcoded fastqs (R1/R2/R3 per lane) for one modality/barcode.
# Depends on the demultiplex outputs (via get_fastq_for_cellranger) so it runs right
# after debarcoding, independent of cellranger; the shell globs the whole barcode
# folder so the synthesized R3 (genomic insert) is QC'd too. fastqc needs java, which
# lives in the 'common' env bin -- prepend it to PATH.
rule fastqc:
    input:
        lambda wildcards: get_fastq_for_cellranger(config['samples'][wildcards.sample]['fastq_path'],sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    output:
        directory('{sample}/{modality}_{barcode}/fastqc')
    params:
        fastq_dir=lambda wildcards: '{sample}/{modality}_{barcode}/fastq/barcode_{barcode}'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode),
        env_bin='/home/mattia/miniconda3_n/envs/common/bin',
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        'mkdir -p {output}; '
        'PATH={params.env_bin}:$PATH fastqc -t {threads} -o {output} {params.fastq_dir}/*.fastq.gz'

# Aggregate every sample's FastQC output into ONE report per snakemake run.
rule multiqc:
    input:
        ['{sample}/{modality}_{barcode}/fastqc'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()]
    output:
        html='multiqc/multiqc_report.html'
    params:
        multiqc='/home/mattia/miniconda3_n/bin/multiqc',
    resources:
        mem_mb = 8000
    shell:
        '{params.multiqc} --force -o multiqc -n multiqc_report.html {input}'

rule get_cells:
    input:
        lambda wildcards: ['{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample = wildcards.sample,
                                                                                            modality = modality, 
                                                                                            barcode = barcodes_dict[wildcards.sample][modality]) 
            for modality in barcodes_dict[wildcards.sample].keys()]
    output:
        cells = '{sample}/all_cells.txt'
    params:
        script = workflow.basedir + '/scripts/get_passed_cells_barcodes.awk'
    shell:
        "cat {input} | awk -f {params.script} | sed 's/\"//g'| sort | uniq > {output.cells}"

rule create_matrix_peaks:
    # Requires /home/mattia/.cargo/bin/fragtk installation (https://github.com/stuart-lab//home/mattia/.cargo/bin/fragtk)
    input:
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        peaks = '{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
        cells = '{sample}/all_cells.txt'
    output:
        features = '{sample}/{modality}_{barcode}/matrix/matrix_peaks/features.tsv.gz',
        matrix   = '{sample}/{modality}_{barcode}/matrix/matrix_peaks/matrix.mtx.gz',
        barcodes = '{sample}/{modality}_{barcode}/matrix/matrix_peaks/barcodes.tsv',
        folder   = directory('{sample}/{modality}_{barcode}/matrix/matrix_peaks/'),
    shell:
        '/home/mattia/.cargo/bin/fragtk matrix -f {input.frag} -b {input.peaks} -c {input.cells} -o {output.folder}'

rule create_matrix_bins:
    input:
        bam   = '{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag  = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        cells = '{sample}/all_cells.txt'
    output:
        features   = '{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/features.tsv.gz',
        matrix     = '{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/matrix.mtx.gz',
        barcodes   = '{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/barcodes.tsv',
        folder     = directory('{sample}/{modality}_{barcode}/matrix/matrix_bin_{bins}/'),
        chromsizes = temp('{sample}/{modality}_{barcode}/chromsizes_{bins}.txt'),
        windows    = temp('{sample}/{modality}_{barcode}/windows_{bins}.txt'),
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    shell:
        """
        samtools idxstats {input.bam} | cut -f1-2 | awk '$2 != 0' > {output.chromsizes}; 
        bedtools makewindows -g {output.chromsizes} -w {wildcards.bins} > {output.windows}; 
        /home/mattia/.cargo/bin/fragtk matrix -f {input.frag} -b {output.windows} -c {input.cells} -o {output.folder}; 
        """
    

rule create_genebody_and_promoter_matrix:
    input:
        cellranger_gtf = config['general']['cellranger_ref'] + 'genes/genes.gtf.gz',
        frag           = '{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        cells          = '{sample}/all_cells.txt',
    output:
        genebody_gtf = '{sample}/{modality}_{barcode}/matrix/matrix_genes/genebody_and_promoter.gtf',
        genebody_bed = '{sample}/{modality}_{barcode}/matrix/matrix_genes/genebody_and_promoter.bed',
        gene_names   = '{sample}/{modality}_{barcode}/matrix/matrix_genes/gene_names.txt',
        features     = '{sample}/{modality}_{barcode}/matrix/matrix_genes/features.tsv.gz',
        matrix       = '{sample}/{modality}_{barcode}/matrix/matrix_genes/matrix.mtx.gz',
        barcodes     = '{sample}/{modality}_{barcode}/matrix/matrix_genes/barcodes.tsv',
        folder       = directory('{sample}/{modality}_{barcode}/matrix/matrix_genes/'),
    conda: '/home/mattia/Multi_nanoCTRNA/envs/nanoscope_general.yaml'
    params:
        script = workflow.basedir + '/scripts/filter_cellranger_gtf_file.py',
    shell:
        'python3 {params.script} -i {input.cellranger_gtf} -o {output.genebody_gtf} -n {output.gene_names};'
        'cut -f 1,4,5 {output.genebody_gtf} > {output.genebody_bed};'
        '/home/mattia/.cargo/bin/fragtk matrix -f {input.frag} -b {output.genebody_bed} -c {input.cells} -o {output.folder}; '


#rule generate_csv:
    #input:
        #script=workflow.basedir + '/scripts/csv_generator.py',
        #metadata='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv'.format(sample=sample,modality=modality,barcode=barcode),
        #fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'.format(sample=sample,modality=modality,barcode=barcode)
    #output:
        #csv="/date/gcb/gcb_MZ/multiNanoCT/merge/{modality}_{barcode}/library.csv"
    #params:
        #sample_data=lambda wildcards:config['samples'][wildcards.sample]
    #conda: '/home/mattia/nanoscope/envs/nanoscope_debarcode.yaml'
    #shell:
        #"python3 {input.script} -o {output.csv} --sample_data {wildcards.sample} --fragments_path {input.fragments} --cells_path {input.metadata} 2>&1"

#rule aggregate_fragments:
    #input:
        #csv="/date/gcb/gcb_MZ/multiNanoCT/samples/FC_Droplet_Paired-Tag/{merge}/{modality_n}_{barcode_n}/library.csv" 
    #output:
        #fragments="{merge}/{modality_n}_{barcode_n}/cellranger/outs/fragments.tsv.gz",
        #peaks="{merge}/{modality_n}_{barcode_n}/cellranger/outs/peaks.bed"
    #params:
       # cellranger_software=config['general']['cellranger_software'],
        #cellranger_ref=config['general']['cellranger_ref'],
       # normalization=config['general']['norm']
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #'rm -rf {wildcards.merge}/{wildcards.modality_n}_{wildcards.barcode_n}/cellranger/; '
        #'cd {wildcards.merge}/{wildcards.modality_n}_{wildcards.barcode_n}/; '
        #'{params.cellranger_software} aggr --id=cellranger --reference={params.cellranger_ref} --normalize={params.normalization} --csv={input.csv} --localcores={threads}  --localmem={resources.mem_gb}'

#rule merge_bw_ATAC:
    #input:
        #bw=['{sample}/ATAC_TATAGCCT/bigwig/all_reads.bw'.format(sample=sample) for sample in samples_list]
    #output:
        #bdg ='{merge}/ATAC_TATAGCCT/bigwig/all_reads.bdg',
        #sorted = '{merge}/ATAC_TATAGCCT/bigwig/all_reads_sorted.bdg',
        #bw_merge = "{merge}/ATAC_TATAGCCT/bigwig/all_reads.bw"
    #params:
        #chrom_sizes = config['general']['chrom_sizes']
    #conda: '/home/mattia/miniconda3/envs/bedtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #/home/mattia/UCSC/bigWigMerge {input.bw} {output.bdg}
        
        #bedtools sort -i {output.bdg} > {output.sorted}
        
        #/home/mattia/UCSC/bedGraphToBigWig {output.sorted} {params.chrom_sizes} {output.bw_merge}
        #"""
#rule merge_bw_H3K27ac:
    #input:
        #bw=['{sample}/H3K27ac_CCTATCCT/bigwig/all_reads.bw'.format(sample=sample) for sample in samples_list]
    #output:
        #bdg ='{merge}/H3K27ac_CCTATCCT/bigwig/all_reads.bdg',
        #sorted = '{merge}/H3K27ac_CCTATCCT/bigwig/all_reads_sorted.bdg',
        #bw_merge = "{merge}/H3K27ac_CCTATCCT/bigwig/all_reads.bw"
    #params:
        #chrom_sizes = config['general']['chrom_sizes']
    #conda: '/home/mattia/miniconda3/envs/bedtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #/home/mattia/UCSC/bigWigMerge {input.bw} {output.bdg}
        
        #bedtools sort -i {output.bdg} > {output.sorted}
        
        #/home/mattia/UCSC/bedGraphToBigWig {output.sorted} {params.chrom_sizes} {output.bw_merge}
        #"""
#rule merge_bw_H3K27me3:
    #input:
        #bw=['{sample}/H3K27me3_ATAGAGGC/bigwig/all_reads.bw'.format(sample=sample) for sample in samples_list]
    #output:
        #bdg ='{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads.bdg',
        #sorted = '{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads_sorted.bdg',
        #bw_merge = "{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads.bw"
    #params:
        #chrom_sizes = config['general']['chrom_sizes']
    #conda: '/home/mattia/miniconda3/envs/bedtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #/home/mattia/UCSC/bigWigMerge {input.bw} {output.bdg}
        
        #bedtools sort -i {output.bdg} > {output.sorted}
        
        #/home/mattia/UCSC/bedGraphToBigWig {output.sorted} {params.chrom_sizes} {output.bw_merge}
        #"""
#rule merge_bam_ATAC:
    #input:
        #bam=['{sample}/ATAC_TATAGCCT/cellranger/outs/possorted_bam.bam'.format(sample=sample) for sample in samples_list]
    #output:
        #merged ='{merge}/ATAC_TATAGCCT/cellranger/outs/possorted_bam.bam'
    #conda: '/home/mattia/miniconda3/envs/samtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #samtools merge {output.merged} {input.bam} -@ {threads}

        #samtools index {output.merged} -@ {threads}

        #"""
#rule merge_bam_H3K27ac:
    #input:
        #bam=['{sample}/H3K27ac_CCTATCCT/cellranger/outs/possorted_bam.bam'.format(sample=sample) for sample in samples_list]
    #output:
        #merged = "{merge}/H3K27ac_CCTATCCT/cellranger/outs/possorted_bam.bam"
    #conda: '/home/mattia/miniconda3/envs/samtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #samtools merge {output.merged} {input.bam}  -@ {threads}

        #samtools index {output.merged} -@ {threads}

        #"""
#rule merge_bam_H3K27me3:
    #input:
        #bam=['{sample}/H3K27me3_ATAGAGGC/cellranger/outs/possorted_bam.bam'.format(sample=sample) for sample in samples_list]
    #output:
        #merged ='{merge}/H3K27me3_ATAGAGGC/cellranger/outs/possorted_bam.bam'
    #conda: '/home/mattia/miniconda3/envs/samtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
       #"""
        
        #samtools merge {output.merged} {input.bam}  -@ {threads}

        #samtools index {output.merged} -@ {threads}
        #"""