//
// Subworkflow with functionality specific to the TEffectbayes pipeline
//

def validateInputParameters() {
    genomeExistsError()
    validateRepeatGtf(params.repeat)
    
    if (!params.genome && (!params.fasta || !params.gtf)) {
        error "You must specify either --genome or both --fasta and --gtf"
    }
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
def validateInputSamplesheet(row) {
    if (row == null || row.sample == null) {
        error "Invalid RNA-seq samplesheet row: ${row}"
    }
    def fastq_files = []
    if (row.fastq_1) fastq_files << row.fastq_1
    if (row.fastq_2) fastq_files << row.fastq_2
    return tuple(row.sample, fastq_files, row.condition, row.cell_line)
}


//
// Validate the samplesheet channel for salmonTE process
//
def validateInputSalmonTESamplesheet(row) {
    if (row == null || row.sample == null) {
        error "Invalid SalmonTE samplesheet row: ${row}"
    }
    def fastq_files = []
    if (row.fastq_1) fastq_files << file(row.fastq_1)
    if (row.fastq_2) fastq_files << file(row.fastq_2)
    return tuple(row.sample, fastq_files)
}

//
// Validate channels from chip samplesheet
//
def validateInputChIPSamplesheet(row) {
    if (!row.antibody || !row.feature_counts || !row.annotation) {
        error "Missing required columns for ChIP sample: ${row}"
    }
    if (!file(row.feature_counts).exists()) error "Feature counts not found for ${row.antibody}"
    if (!file(row.annotation).exists()) error "Annotation not found for ${row.antibody}"
    return tuple(row.antibody, file(row.feature_counts), file(row.annotation))
}

//
// Validate repeat gtf
//
def validateRepeatGtf(path_str) {
    def f = file("${workflow.projectDir}/${path_str}")
    if (!f.exists()) {
        error "ERROR: repeat GTF file not found at path: ${path_str}"
    }
    if (!path_str.toLowerCase().endsWith(".gtf")) {
        error "ERROR: repeat file must be a .gtf file. Given: ${path_str}"
    }
    return f
}