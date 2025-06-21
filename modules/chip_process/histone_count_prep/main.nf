process histone_count_prep {
    tag "$antibody"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pybedtools:latest'

    publishDir "${params.outdir}/histone_quant", mode: 'copy'

    input:
    tuple val(antibody), path(feature_counts), path(annotation)
    
    output:
    path "histone_quant/${antibody}_counttable.csv", emit: histone_counttable

    script:
    """
    python3 <<EOF
    import os
    import pandas as pd

    feature_counts_file = "${feature_counts}"
    annotated_peaks_file = "${annotation}"
    antibody = "${antibody}"
    replicate_count = ${params.replicate_count}

    output_dir = "histone_quant"
    os.makedirs(output_dir, exist_ok=True)

    # --- Read and Normalize Counts ---
    df = pd.read_csv(feature_counts_file, sep="\\t", comment="#", low_memory=False)
    count_columns = df.columns[6:]
    df[count_columns] = df[count_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
    scaling_factors = df[count_columns].sum() / 1_000_000
    rpm_table = df.copy()
    rpm_table[count_columns] = rpm_table[count_columns].div(scaling_factors, axis=1)

    # --- Annotated Peaks Merge ---
    df_ann = pd.read_csv(annotated_peaks_file, sep="\\t")
    df_ann.rename(columns={df_ann.columns[0]: "Peak_ID"}, inplace=True)
    df_final = rpm_table.merge(df_ann[['Peak_ID', 'Gene Name', 'Chr', 'Start', 'End', 'Strand']], 
                               left_on="Geneid", right_on="Peak_ID").drop(columns=["Geneid"])

    # --- Summarize Replicates ---
    replicate_columns = df_final.columns[7:]
    samples = replicate_columns.str.extract(r"(.+)_Rep\\d")[0]
    summarized_rpm = df_final.iloc[:, :7].copy()
    for sample in samples.unique():
        summarized_rpm[sample] = df_final[replicate_columns[samples == sample]].sum(axis=1)

    # --- Melt to Long Format ---
    long_chipseq = summarized_rpm.melt(
        id_vars=["Peak_ID", "Gene Name", "Chr", "Start", "End", "Strand"],
        var_name="cell_line",
        value_name="chip_rpm"
    )

    # --- Expand Based on RNA Replicate Count ---
    expanded_chipseq = long_chipseq.loc[long_chipseq.index.repeat(replicate_count)].reset_index(drop=True)

    # --- Output File ---
    output_file = os.path.join(output_dir, f"{antibody}_counttable.csv")
    expanded_chipseq.to_csv(output_file, index=False)
    EOF

    python3 --version | awk '{print $2}' > versions.yml
    """
}