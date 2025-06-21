process BNM_INPUT_PREP_STEP3 {
    tag "bnm_input"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'

    publishDir "${params.outdir}/BNM/input/", mode: 'copy'
    
    input:
    path samplesheet
    path repeat_histone_intersect
    path cell_line_csvs

    output:
    path "BNM/input/*.csv", emit: bnm_input_csvs
    path "BNM/input/*.pkl", emit: bnm_input_pickles

    script:
    """
    mkdir -p BNM/input

    python <<EOF
    import os
    import re
    import glob
    import pickle
    import pandas as pd
    
    # Load samplesheet and get unique cell lines
    samplesheet = pd.read_csv("${samplesheet}")
    cell_lines = samplesheet['cell_line'].unique()

    # Define the input and output directories
    overlaps_dir = "${cell_line_csvs}"
    intersect_dir = "${repeat_histone_intersect}"
    output_dir = "BNM/input"
    os.makedirs(output_dir, exist_ok=True)

    # Function to dynamically detect repeat count
    def get_repeat_count(df):
        return max(
            int(col.split("_")[0].replace("repeat", ""))
            for col in df.columns if col.startswith("repeat") and "_start" in col
        )

    # Function to check interval overlap
    def check_overlap(promoter_repeat_df, df_hist, repeat_num, histone_name):
        repeat_start_col = f"repeat{repeat_num}_start"
        repeat_end_col = f"repeat{repeat_num}_end"

        # Merge on 'Chromosome' and check interval overlap
        merged = promoter_repeat_df.merge(df_hist, on="Chromosome", suffixes=("", "_histone"))

        # Condition for overlap: (repeat_start <= histone_start) & (repeat_end >= histone_end)
        overlap_condition = (merged[repeat_start_col] <= merged["End"]) & (merged[repeat_end_col] >= merged["Start"])

        # Get rows where overlap exists
        overlapping_rows = merged.loc[overlap_condition, ["Chromosome", repeat_start_col, repeat_end_col]]

        # Mark overlap in original df
        promoter_repeat_df.loc[
            promoter_repeat_df.set_index(["Chromosome", repeat_start_col, repeat_end_col])
            .index.isin(overlapping_rows.set_index(["Chromosome", repeat_start_col, repeat_end_col]).index),
            f"{histone_name}_repeat{repeat_num}"
        ] = 1

    # Loop over each cell line
    for cell_line in cell_lines:
        print(f"Processing cell line: {cell_line}")

        # Load Promoter-Repeat Overlap Data
        promoter_repeat_path = os.path.join(overlaps_dir, f"intersect_{cell_line}.csv")
        promoter_repeat_df = pd.read_csv(promoter_repeat_path)

        # Detect number of repeats dynamically
        repeat_count = get_repeat_count(promoter_repeat_df)
        print(f"Detected {repeat_count} repeats for {cell_line}")

        # Load Repeat-Histone Intersection Tables
        histone_files = glob.glob(os.path.join(intersect_dir, f"{cell_line}_INTERSECT_*.csv"))

        # Create a dictionary for histone data
        histone_tables = {
            os.path.basename(f).split("_")[-1].replace(".csv", "").lower(): pd.read_csv(f)
            for f in histone_files
        }

        # Initialize new columns for overlap results
        for histone in histone_tables.keys():
            for i in range(1, repeat_count + 1):
                promoter_repeat_df[f"{histone}_repeat{i}"] = 0

        # Apply the overlap check for each histone table and repeat interval
        for histone, df_hist in histone_tables.items():
            print(f"  Processing histone: {histone}")
            for i in range(1, repeat_count + 1):
                print(f"    Checking overlap for Repeat{i}")
                check_overlap(promoter_repeat_df, df_hist, i, histone)

        # Save the final table
        output_path = os.path.join(output_dir, f"{cell_line}_overlaps_final.csv")
        promoter_repeat_df.to_csv(output_path, index=False)
        print(f"Processing complete! Final table saved to {output_path}")

    # Define output directory
    final_csv_path = os.path.join(output_dir, "bnm_input.csv")
    final_pickle_path = os.path.join(output_dir, "bnm_input.pkl")
    EOF
    """
}