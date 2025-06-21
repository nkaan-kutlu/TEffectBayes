process BNM_INPUT_PREP_STEP1 {
    tag "gene_repeat_intersect_transformed"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'

    publishDir "${params.outdir}/Data_Integration/gene_repeat_intersect/transformed/", mode: 'copy'
    
    input:
    path gene_repeat_bedtools
    val chip_antibodies

    output:
    path "Data_Integration/gene_repeat_intersect/transformed/*_intersect_transformed.csv", emit: gene_repeat_intersect_transformed

    script:
    """
    mkdir -p Data_Integration/gene_repeat_intersect/transformed

    python <<EOF
    import os
    import pandas as pd
    import numpy as np
    import glob

    input_dir = "${gene_repeat_bedtools}"
    out_dir = "Data_Integration/gene_repeat_intersect/transformed"
    os.makedirs(out_dir, exist_ok=True)
    
    # Define column names
    antibodies = "${chip_antibodies}".strip("[]").replace("'", "").split(",")
    rpm_columns = [f"{ab.strip()}_rpm" for ab in antibodies]
    
    base_columns = [
    "Chromosome", "Start", "End",
    "gene_id", "gene_name", "gene_type",
    "sample", "cell_line", "FPKM"]
    
    repeat_columns = [
    "Chromosome_repeat", "Start_repeat", "End_repeat",
    "repeat_id", "TPM", "sample_repeat"]
    
    column_names = base_columns + rpm_columns + repeat_columns
    # Get all intersect.txt files
    intersect_files = glob.glob(os.path.join(input_dir, "*_intersect.txt"))

    # Loop over each file
    for file in intersect_files:
        # Load the data
        df = pd.read_csv(file, sep="\t", header=None, names=column_names)

        # Replace invalid values with NaN in repeat columns
        invalid_values = [".", "-1", None]
        df[["Chromosome_repeat", "Start_repeat", "End_repeat", "repeat_id", "TPM", "sample_repeat"]] = df[["Chromosome_repeat", "Start_repeat", "End_repeat", "repeat_id", "TPM", "sample_repeat"]].replace(invalid_values, pd.NA)

        # Identify and print genes with no repeat intersections
        missing_repeat_info = df[df[["Chromosome_repeat", "Start_repeat", "End_repeat", "repeat_id", "TPM", "sample_repeat"]].isnull().any(axis=1)]
        print(f"Genes with missing repeat information in {file}:")
        print(missing_repeat_info["gene_name"].unique())

        # Eliminate rows of genes with no repeats
        df_cleaned = df.dropna(subset=["Chromosome_repeat", "Start_repeat", "End_repeat", "repeat_id", "TPM", "sample_repeat"])

        # Drop redundant repeat columns
        df_cleaned = df_cleaned.drop(columns=["Chromosome_repeat", "sample_repeat"])

        # Replace NaN values in histone_rpm columns with '0'
        df_cleaned[rpm_columns] = df_cleaned[rpm_columns].fillna('0')

        # Group by relevant columns and aggregate the repeat info
        group_cols = [
        "Chromosome", "Start", "End",
        "gene_id", "gene_name", "gene_type",
        "sample", "cell_line", "FPKM"] + rpm_columns

        grouped = df_cleaned.groupby(
            group_cols,
            observed=False).apply(lambda group: group[['Start_repeat', 'End_repeat', 'repeat_id', 'TPM']].to_dict(orient="records"))

        # Reset the index to transform the grouped data into a DataFrame
        grouped_df = grouped.reset_index(name='repeats')

        # Expand the repeats into separate columns
        final_data = []

        for _, row in grouped_df.iterrows():
            base_data = {
                "Chromosome": row["Chromosome"],
                "Start": row["Start"],
                "End": row["End"],
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "gene_type": row["gene_type"],        
                "sample": row["sample"],
                "cell_line": row["cell_line"],
                "gene_FPKM": row["FPKM"]
            }

            # Dynamically add RPM columns
            for rpm_col in rpm_columns:
                base_data[rpm_col] = row[rpm_col]

            # Add repeat information to the same row
            for i, repeat in enumerate(row["repeats"], start=1):
                base_data.update({                
                    f"repeat{i}_start": repeat["Start_repeat"],
                    f"repeat{i}_end": repeat["End_repeat"],
                    f"repeat{i}_id": repeat["repeat_id"],
                    f"repeat{i}_TPM": repeat["TPM"],
                })
        
            final_data.append(base_data)

        final_df = pd.DataFrame(final_data)
    
        # Fill NaN values temporarily for the repeat columns before conversion
        repeat_columns = [col for col in final_df.columns if 'repeat' in col and ('start' in col or 'end' in col)]
        final_df[repeat_columns] = final_df[repeat_columns].fillna(0)
    
        # Now convert the repeat start and end columns to numeric (either float or integer)
        for col in repeat_columns:
            final_df[col] = pd.to_numeric(final_df[col], errors='coerce', downcast='integer')  # Use 'integer' for conversion
        
        final_df[repeat_columns] = final_df[repeat_columns].replace(0, pd.NA)
    
        # Identify repeat columns dynamically
        repeat_tpm_cols = [col for col in final_df.columns if "repeat" in col and "_TPM" in col]
        repeat_id_cols = [col for col in final_df.columns if "repeat" in col and "_id" in col]
    
        # Convert TPM columns to nullable integer type
        final_df[repeat_tpm_cols] = final_df[repeat_tpm_cols].astype("Float64")  # Pandas nullable integer
    
        # Convert repeat_id columns to Pandas nullable string type
        final_df[repeat_id_cols] = final_df[repeat_id_cols].astype("string")
    
        # Save the final DataFrame to a CSV file
        output_filename = os.path.join(out_dir, f"{os.path.basename(file).replace('.txt', '_intersect_transformed.csv')}")
        final_df.to_csv(output_filename, sep=",", index=False)
        print(f"Saved processed data to {output_filename}")
    
    EOF
    """

}
