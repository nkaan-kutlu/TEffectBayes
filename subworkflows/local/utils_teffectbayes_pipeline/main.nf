//
// Subworkflow with functionality specific to the TEffectbayes pipeline
//

def validateInputParameters() {
    getGenomeAttribute ()
    genomeExistsError()
    validateInputSamplesheet()
    validateInputSalmonTESamplesheet()
    validateInputChIPSamplesheet()
    validateRepeatGtf()
}

//
// Get attribute from genome config file
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    } else if (attribute == 'fasta' && params.fasta) {
        return params.fasta
    } else if (attribute == 'gtf' && params.gtf) {
        return params.gtf
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }

}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Validate the samplesheet channle for salmonTE process
//
def validateInputSalmonTESamplesheet(input) {
    def sample_id = row.sample
    def fastq_files = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
    return tuple(sample_id, fastq_files)
}

//
// Validate channels from chip samplesheet
//
def validateInputChIPSamplesheet(input_chip) {
     def (row) = input_chip
    if (!row.antibody || !row.feature_counts || !row.annotation) {
        error "ERROR: Missing required column(s) for sample ${row.antibody}. All columns (antibody, feature_counts, annotation) must be provided."
    }

    // Check for the count and annotation files
    if (!file(row.feature_counts).exists()) {
        error "ERROR: Feature counts file not found for sample ${row.antibody} at path: ${row.feature_counts}"
    }

    if (!file(row.annotation).exists()) {
        error "ERROR: Annotation file not found for sample ${row.antibody} at path: ${row.annotation}"
    }

    return [ row.antibody, file(row.feature_counts), file(row.annotation) ]
}

//
// Validate repeat gtf
//
def validateRepeatGtf(path_str) {
    def f = file(path_str)
    if (!f.exists()) {
        error "ERROR: repeat GTF file not found at path: ${path_str}"
    }
    if (!path_str.toLowerCase().endsWith(".gtf")) {
        error "ERROR: repeat file must be a .gtf file. Given: ${path_str}"
    }
    return f
}