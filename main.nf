/*
 * HELP MESSAGE
 */

if (params.help) {
    log.info """
    How to use :

    Required Arguments:
    --input         RNA-seq input samplesheet
    --input_chip    ChIP-seq input samplesheet
    --repeat        Repeat annotation gtf file
    --replicate_count   The replicate sample number for ChIP-seq samples that needs to be specified 

    Optional Arguments: 
    --gtf                Annotation gtf file if genome is not specified
    --fasta              fasta file if genome is not specified  
    --genelist           target genelist can be provided if desired
    --interval           the genomic range for the gene upstream region dafault: 5000 kb
   """
    exit 1
}

println """


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { validateInputParameters  } from './subworkflows/local/utils_teffectbayes_pipeline'

include { STAR_GENOMEGENERATE } from './modules/star/genomegenerate/main.nf'
include { STAR_ALIGN } from './modules/star/align/main.nf'
include { SAMTOOLS_SORT } from './modules/samtools/sort/main.nf'
include { featureCountsR } from './modules/rsubread/featureconts/main.nf'

include { salmonTE_quant } from './modules/salmonTE/quant/main.nf'
include { locus_quant } from './modules/salmonTE/locus_quant/main.nf'
include { locus_interval_prep } from './modules/salmonTE/locus_interval_prep/main.nf'

include { histone_count_prep } from './modules/chip_process/histone_count_prep/main.nf'
include { histone_interval_prep } from './modules/chip_process/histone_interval_prep/main.nf'

include { gene_histone_intersection } from './modules/data_integration/gene_histone/main.nf'
include { repeat_histone_intersection } from './modules/data_integration/repeat_histone/main.nf'
include { gene_repeat_intervals } from './modules/data_integration/gene_repeat/gene_repeat_intervals/main.nf'
include { BEDTOOLS_INTERSECT } from './modules/data_integration/gene_repeat/bedtools/intersect/main.nf'

include { BNM_INPUT_PREP_STEP1 } from './modules/bnm_input/bnm_input_prep_step1/main.nf'
include { BNM_INPUT_PREP_STEP2 } from './modules/bnm_input/bnm_input_prep_step2/main.nf'
include { BNM_INPUT_PREP_STEP3 } from './modules/bnm_input/bnm_input_prep_step3/main.nf'
include { BNM_INPUT_PREP_STEP4 } from './modules/bnm_input/bnm_input_prep_step4/main.nf'

include { BNM_CALCULATION } from './modules/bnm/bnm_calculation/main.nf'
include { BNM_INFERENCE } from './modules/bnm/bnm_inference/main.nf'
include { BNM_VISUALIZATION } from './modules/bnm/bnm_visualization/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.input            = params.input ?: null
params.input_chip       = params.input_chip ?: null
params.genome           = params.genome ?: null
params.fasta            = params.fasta ?: null
params.gtf              = params.gtf ?: null
params.repeat           = params.repeat ?: null
params.replicate_count  = params.replicate_count ?: null
params.interval         = params.interval ?: null
params.genelist         = params.genelist ?: null
params.outdir           = params.outdir ?: "results"

workflow {
    if (!params.input) error "ERROR: Input samplesheet must be provided with --input"
    if (!params.input_chip) error "ERROR: Input ChIP-Seq samplesheet must be provided with --input_chip"
    if (!params.repeat) error "ERROR: Please provide a repeat GTF file using --repeat"
    if (!params.genome && (!params.fasta || !params.gtf)) {exit 1, "You must specify either --genome or both --fasta and --gtf"}

    // Validate general input parameters
    validateInputParameters()

    // Load genome attribute files
    genome_fasta_ch = Channel.fromPath(getGenomeAttribute('fasta')).map { [null, it] }
    genome_gtf_ch   = Channel.fromPath(getGenomeAttribute('gtf')).map { [null, it] }
    repeat_gtf_ch = Channel.value(params.repeat).map { validateRepeatGtf(it) }

    // Loading the samplesheets to validate and create channels
    samplesheet_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ",")
        .map { row -> validateInputSamplesheet(row) }

    salmonTE_samplesheet_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ",")
        .map { row -> validateInputSalmonTESamplesheet(row) }

    chip_samplesheet_ch = Channel.fromPath(params.input_chip)
        .splitCsv(header: true, sep: ",")
        .map { row -> validateInputChIPSamplesheet(row) }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Gene Expression Analysis
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    star_index_ch = STAR_GENOMEGENERATE(genome_fasta_ch, genome_gtf_ch)

    aligned_bam_ch = STAR_ALIGN(samplesheet_ch, star_index_ch, genome_gtf_ch)

    sorted_bam_ch = SAMTOOLS_SORT(aligned_bam_ch, genome_fasta_ch)

    gene_quant_ch = featureCountsR(sorted_bam_ch, genome_gtf_ch)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Repeat Expression Analysis
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    salmonte_out_ch = salmonTE_quant(salmonTE_samplesheet_ch)

    locus_out_ch = locus_quant(salmonte_out_ch.map { it.quant_results }, repeat_gtf_ch)

    cell_line_map_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ",")
        .map { row -> tuple(row.sample, row.cell_line) }

    te_intervals_ch = locus_interval_prep(
        locus_out_ch.map { it.te_counttable },
        cell_line_map_ch
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Histone Expression Analysis
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
	
	chip_antibodies_ch = chip_samplesheet_ch.map { it[0] }.collect()
    
	chip_counts_ch = histone_count_prep(chip_samplesheet_ch)

    chip_intervals_ch = histone_interval_prep(chip_samplesheet_ch)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Data Integration
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    chip_csvs = Channel.fromPath("${params.outdir}/histone_quant/*_counttable.csv")

    gene_histone_out_ch = gene_histone_intersection(
        genes_counttable: file("${params.outdir}/gene_quant/genes_counttable.csv"),
        chip_counttables: chip_csvs,
        samplesheet: file(params.input),
        genelist: params.genelist ? file(params.genelist) : null
    )

    repeat_histone_ch = repeat_histone_intersection(
        repeat_intervals: te_intervals_ch,
        chip_intervals: chip_intervals_ch,
        samplesheet: file(params.input)
    )

    gene_repeat_beds_ch = gene_repeat_intervals(
        gene_histone_counttable: gene_histone_out_ch,
        te_counttable: repeat_histone_ch
    )

    gene_repeat_beds_ch
        .map { file -> 
            def fname = file.getBaseName()
            def sample_id = fname.tokenize("_")[0]
            def type = fname.contains("promoter") ? "promoter" : "repeat"
            tuple(sample_id, type, file)
        }
        .groupTuple()
        .filter { sample_id, files -> 
            files*.get(1).containsAll(["promoter", "repeat"])
        }
        .map { sample_id, files -> 
            def promoter = files.find { it[1] == "promoter" }[2]
            def repeat = files.find { it[1] == "repeat" }[2]
            def meta = [id: sample_id]
            tuple(meta, promoter, repeat)
        }
        .set { paired_gene_repeat_beds_ch }

    gene_repeat_bedtools_ch = BEDTOOLS_INTERSECT(paired_gene_repeat_beds_ch)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Bayesian Network Modelling
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    gene_repeat_intersect_transformed_ch = BNM_INPUT_PREP_STEP1(gene_repeat_bedtools_ch, chip_antibodies_ch)

    cell_line_csvs_ch = BNM_INPUT_PREP_STEP2(gene_repeat_intersect_transformed_ch)

    bnm_input_pickles_ch = BNM_INPUT_PREP_STEP3(samplesheet: file(params.input),repeat_histone_ch,cell_line_csvs_ch)

    bnm_gene_input_pickles_ch = BNM_INPUT_PREP_STEP4(samplesheet: file(params.input),samplesheet_chip: file(params.input_chip),bnm_input_pickles_ch)

    bnm_calc_ch = BNM_CALCULATION(samplesheet: file(params.input),samplesheet_chip: file(params.input_chip),bnm_gene_input_pickles_ch)

    bnm_inf_ch = BNM_INFERENCE(bnm_calc_ch)

    bnm_visual_ch = BNM_VISUALIZATION(bnm_calc_ch)

}