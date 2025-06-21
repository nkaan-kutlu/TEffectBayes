process featureCountsR {
    tag "$meta.id (featureCounts)"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/rtools:latest'

    publishDir "${params.outdir}/gene_quant", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    path "gene_quant/counts.tsv", emit: counts
    path "gene_quant/FPKM.tsv", emit: fpkm
    path "gene_quant/counts_annotated.csv", emit: counts_annotated
    path "gene_quant/FPKM_annotated.csv", emit: fpkm_annotated
    path "gene_quant/genes_counttable.csv", emit: genes_counttable

    script:
    """
    mkdir -p gene_quant

    Rscript bin/featureCounts.R \\
        --bam "${bam}" \\
        --gtf "${gtf}" \\
        --outdir "gene_quant"

    cat <<-END_VERSIONS > gene_quant/versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1 | awk '{print \$3}')
    END_VERSIONS
    """

    
}
