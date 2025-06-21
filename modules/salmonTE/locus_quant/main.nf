process locus_quant {
  tag "Locus-specific TE quantification"
  label 'process_high_memory'

  conda "${moduleDir}/environment.yml"
  container 'nkaankutlu/rtools:latest'
    
  publishDir "${params.outdir}/TE_quant", mode: 'copy'
  
  input:
  path quant_dirs
  path repeat_gtf

  output:
  path "TE_quant/ls_TE_counttable.csv", emit: te_counttable

  script:
  """
    Rscript bin/locus_specific_te_quant.R \\
      --repeat_gtf ${repeat_gtf} \\
      --quant_dirs ${quant_dirs.join(" ")} \\
      --output TE_quant/ls_TE_counttable.csv

   echo \\
    \"\"\"\\
    process: locus_quant
    r_version: \$(R --version | grep 'R version' | awk '{print \$3}')
    container: quay.io/biocontainers/r:4.1.0--h5f2f3f3_0
    \"\"\" > TE_quant/versions.yml
  """
}