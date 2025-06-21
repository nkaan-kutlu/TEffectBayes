process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/nf-core/htslib_samtools_star_gawk:311d422a50e6d829'

    publishDir "${params.outdir}/align/sorted", mode: 'copy'

    input:
    tuple val(meta), path(bam_unsorted)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"),  emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.csi"),  emit: csi,  optional: true
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools cat \\
        ${bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.${extension} \\
        -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}