process histone_interval_prep {
    tag "$antibody"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pybedtools:latest'

    publishDir "${params.outdir}/histone_intervals", mode: 'copy'

    input:
    tuple val(antibody), path(feature_counts), path(annotation)

    output:
    path "histone_intervals/*_${antibody}_intervals.csv", emit: histone_intervals

    script:
    """
    python <<EOF
    import os
    import pandas as pd

    antibody = "${antibody}"
    feature_counts_file = "${feature_counts}"
    annotated_peaks_file = "${annotation}"

    output_dir = "histone_intervals"
    os.makedirs(output_dir, exist_ok=True)

    # Read featureCounts
    df = pd.read_csv(feature_counts_file, sep="\\t", comment="#", low_memory=False)
    df = df[['Geneid'] + list(df.columns[6:])]

    # Read annotation
    df_ann = pd.read_csv(annotated_peaks_file, sep="\\t")
    df_ann.rename(columns={df_ann.columns[0]: "Peak_ID"}, inplace=True)
    df_ann = df_ann[['Peak_ID', 'Chr', 'Start', 'End', 'Strand', 'Gene Name']]

    # Merge
    df_final = pd.merge(df, df_ann, left_on="Geneid", right_on="Peak_ID").drop(columns=["Geneid"])
    df_final = df_final[["Peak_ID", "Gene Name", "Chr", "Start", "End", "Strand"] +
                        [col for col in df_final.columns if col not in ["Peak_ID", "Gene Name", "Chr", "Start", "End", "Strand"]]]

    # Extract sample names
    df_counts = df_final.copy()
    sample_names = list(set(col.split(f"_{antibody}_Rep")[0] for col in df_counts.columns[6:]))

    # Process each sample
    for sample in sample_names:
        sample_columns = [col for col in df_counts.columns if sample in col]
        df_sample = df_counts[['Peak_ID', 'Gene Name', 'Chr', 'Start', 'End', 'Strand'] + sample_columns].copy()
        df_sample = df_sample[(df_sample[sample_columns] > 0).all(axis=1)]
        df_sample.drop(columns=sample_columns, inplace=True)

        output_file = os.path.join(output_dir, f"{sample}_${antibody}_intervals.csv")
        df_sample.to_csv(output_file, index=False)
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: $(python --version 2>&1 | awk '{print $2}')
        pandas: $(python -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
    """
}
