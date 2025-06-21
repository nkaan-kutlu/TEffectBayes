process salmonTE_quant {
    tag "TE_quant"

    conda "${moduleDir}/environment.yml"    
    container 'nkaankutlu/salmonTE:latest'
    
    publishDir "${params.outdir}/TE_quant", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_files)
    path outdir

    output:
    path "TE_quant/", emit: quant_dirs

    script:
    """
    mkdir -p TE_quant/
    export PATH=\$PATH:SalmonTE
    SalmonTE.py quant --reference=hs --outpath TE_quant/ \\
        ${fastq_files.join(" ")}
    """
}
