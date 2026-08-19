include: 'Snakefile_prep_RNA.smk'

rule all_preprocess:
    input:
        bam_RNA=[
            '{sample}/RNA_AAAAGGGG/cellranger/outs/possorted_genome_bam.bam'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()]
       

rule run_cellranger_RNA:
    input:
        lambda wildcards: get_fastq_for_cellranger_rna(config['samples'][wildcards.sample]['fastq_path_RNA']+ '/**/*{lane}*R[12]*.fastq.gz',sample=wildcards.sample)
    output:
        bam_RNA='{sample}/RNA_AAAAGGGG/cellranger/outs/possorted_genome_bam.bam'
    params:
        cellranger_software=config['general']['cellranger_software_RNA'],
        cellranger_ref=config['general']['cellranger_ref_RNA'],
        fastq_folder=lambda wildcards: config['samples'][wildcards.sample]['fastq_path_RNA'],
        mem=64
    threads: 17
    resources:
        mem_mb = 64000
    shell:
        'rm -rf {wildcards.sample}/RNA_AAAAGGGG/cellranger/; '
        'cd {wildcards.sample}/RNA_AAAAGGGG/; '
        '{params.cellranger_software} count --id cellranger --transcriptome {params.cellranger_ref} --chemistry=ARC-v1 --fastqs {params.fastq_folder} --localcores={threads} --localmem={params.mem} --create-bam=true' 

