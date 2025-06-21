process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/nf-core/htslib_samtools_star_gawk:311d422a50e6d829'

    publishDir "${params.outdir}/align", mode: 'copy'

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    val star_ignore_sjdbgtf

    output:
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    

    tuple val(meta), path('*d.out.bam')                              , optional:true, emit: bam
    tuple val(meta), path("${prefix}.sortedByCoord.out.bam")         , optional:true, emit: bam_sorted
    tuple val(meta), path("${prefix}.Aligned.sortedByCoord.out.bam") , optional:true, emit: bam_sorted_aligned
    tuple val(meta), path('*toTranscriptome.out.bam')                , optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam')                 , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')                               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                                   , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')                            , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')                  , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')                          , optional:true, emit: junction
    tuple val(meta), path('*.out.sam')                               , optional:true, emit: sam
    tuple val(meta), path('*.wig')                                   , optional:true, emit: wig
    tuple val(meta), path('*.bg')                                    , optional:true, emit: bedgraph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }
    def ignore_gtf      = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' 'SM:$prefix'"
    def out_sam_type    = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $ignore_gtf \\
        $attrRG \\
        $args

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}
