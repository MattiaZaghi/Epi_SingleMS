include: 'Snakefile_prep.smk'


rule all_preprocess:
    input:
        cellbender_seurat_h5=[
            '/date/gcb/gcb_MZ/multiNanoCT/samples/Human_Brain_doplet_paired_Tag/{sample}/RNA_AAAAGGGG_H3K27ac/cellbender_output_seurat.h5'.format(sample=sample) for sample in samples_list]
       
rule run_cellbender:
    input:
        h5='/date/gcb/gcb_MZ/multiNanoCT/samples/Human_Brain_doplet_paired_Tag/{sample}/RNA_AAAAGGGG_H3K27ac/cellranger/outs/raw_feature_bc_matrix.h5'
    output:
        cellbender_h5='/date/gcb/gcb_MZ/multiNanoCT/samples/Human_Brain_doplet_paired_Tag/{sample}/RNA_AAAAGGGG_H3K27ac/cellbender_output.h5'
    params:
        fpr=config['general']['cellbender_fpr'],
        cuda_flag=lambda wildcards: '--cuda' if config['general']['cellbender_cuda'] else ''
    threads: 1
    resources:
        mem_mb = 16000,
        gpu=1
    conda:
        '/home/mattia/miniconda3/envs/cellbender.yml'
    shell:
        'cellbender remove-background '
        '--input {input.h5} '
        '--output {output.cellbender_h5} '
        '--fpr {params.fpr} ' 
        '{params.cuda_flag}'


rule ptrepack_seurat:
    input:
        cellbender_h5='/date/gcb/gcb_MZ/multiNanoCT/samples/Human_Brain_doplet_paired_Tag/{sample}/RNA_AAAAGGGG_H3K27ac/cellbender_output.h5'
    output:
        seurat_h5='/date/gcb/gcb_MZ/multiNanoCT/samples/Human_Brain_doplet_paired_Tag/{sample}/RNA_AAAAGGGG_H3K27ac/cellbender_output_seurat.h5'
    threads: 1
    resources:
        mem_mb = 16000
    conda:
        '/home/mattia/miniconda3/envs/cellbender.yml'
    shell:
        'ptrepack --complevel 5 {input.cellbender_h5}:/matrix {output.seurat_h5}:/matrix'
