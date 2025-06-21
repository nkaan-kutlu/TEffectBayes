process BEDTOOLS_INTERSECT {
    tag "$sample_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

    publishDir "${params.outdir}/Data_Integration/gene_repeat_intersect/intersect/", mode: 'copy'

    input:
    tuple val(sample_id), path(promoter_bed), path(repeat_bed)

    output:
    path("Data_Integration/gene_repeat_intersect/intersect/${sample_id}_intersect.txt"), emit: gene_repeat_bedtools

    script:
    """
    mkdir -p Data_Integration/gene_repeat_intersect/intersect

    bedtools intersect \\
        -a ${promoter_bed} \\
        -b ${repeat_bed} \\
        -wa -wb -loj -header \\
        > Data_Integration/gene_repeat_intersect/intersect/${sample_id}_intersect.txt
    """
}