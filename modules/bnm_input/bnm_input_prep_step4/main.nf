process BNM_INPUT_PREP_STEP4 {
    tag "bnm_gene_input"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'

    publishDir "${params.outdir}/BNM/gene_input/", mode: 'copy'
    
    input:
    path samplesheet
    path samplesheet_chip
    path bnm_input_pickles
    
    output:
    path "BNM/gene_input/*.csv", emit: bnm_gene_input_csvs
    path "BNM/gene_input/*.pkl", emit: bnm_gene_input_pickles

    script:
    """
    mkdir -p BNM/gene_input

    python <<EOF
    import os
    import re
    import pickle
    import pandas as pd
    
    # Read samplesheet
    samplesheet_df = pd.read_csv("${samplesheet}")

    # Create dictionary: condition -> set of cell lines
    cell_line_groups = samplesheet_df.groupby("condition")["cell_line"].apply(set).to_dict()

    # Read chip samplesheet
    chip_df = pd.read_csv("${samplesheet_chip}")

    # Get unique antibody names (case-insensitive matching for safety)
    antibody_names = chip_df['antibody'].dropna().unique().tolist()

    input_file = ("${bnm_input_pickles}")
    output_dir = "BNM/gene_input"
    os.makedirs(output_dir, exist_ok=True)

    def process_file(file_path):
        try:
            # Load the pickle file
            with open(file_path, "rb") as f:
                final_df = pickle.load(f)

            # Assign biological_condition based on cell_line
            def assign_biological_condition(cell_line):
                for condition, lines in cell_line_groups.items():
                    if cell_line in lines:
                        return condition
                return "Unknown"

            final_df['biological_condition'] = final_df['cell_line'].apply(assign_biological_condition)
            
            # Extract unique gene names
            unique_gene_names = final_df['gene_name'].unique()

            # Loop over each unique gene name and process it
            for gene in unique_gene_names:
                # Filter rows for the current gene name
                gene_df = final_df[final_df['gene_name'] == gene].copy()
                
                # Rename the `gene_FPKM` column to `{gene_name}_FPKM`
                if 'gene_FPKM' in gene_df.columns:
                    new_column_name = f"{gene}_FPKM"
                    gene_df.rename(columns={'gene_FPKM': new_column_name}, inplace=True)
                
                # Identify repeat-related columns
                repeat_id_cols = [col for col in gene_df.columns if col.startswith("repeat") and col.endswith("_id")]
                repeat_tpm_cols = [col for col in gene_df.columns if col.startswith("repeat") and col.endswith("_TPM")]
                histone_repeat_cols = [col for col in gene_df.columns if any(h.lower() in col.lower() for h in antibody_names) and "repeat" in col.lower()]
        
                print("Histone Columns:", histone_repeat_cols)
                
                # Replace NaN values in histone columns with 0
                gene_df[histone_repeat_cols] = gene_df[histone_repeat_cols].fillna(0)
                
                # Drop repeat columns that have NaN values in both `repeat_id` and `repeat_TPM`
                valid_repeats = []
                for repeat_id_col, repeat_tpm_col in zip(repeat_id_cols, repeat_tpm_cols):
                    if not gene_df[[repeat_id_col, repeat_tpm_col]].isnull().all().all():
                        valid_repeats.append((repeat_id_col, repeat_tpm_col))
                
                # Create a mapping of repeat_id to TPM column
                repeat_mapping = {}
                for repeat_id_col, repeat_tpm_col in valid_repeats:
                    repeat_id = gene_df[repeat_id_col].iloc[0]  # Get repeat name
                    if repeat_id not in repeat_mapping:
                        repeat_mapping[repeat_id] = {
                            "TPM": repeat_tpm_col,
                            "histones": histone_repeat_cols[:]
                        }
                
                    m = re.search(r'repeat(\d+)_TPM', repeat_tpm_col)
                    if m: 
                        repeat_num = m.group(1)
                        # Filter histone_repeat_cols to only include those with the matching repeat number.
                        filtered_histones = [h for h in histone_repeat_cols if h.endswith(f"repeat{repeat_num}")]
                        repeat_mapping[repeat_id]["histones"] = filtered_histones
                        
                        print(f"Repeat: {repeat_id}, Repeat number: {repeat_num}, Filtered Histones: {filtered_histones}")
                
                print("Repeat Mapping:", repeat_mapping)
                
                # Generate dynamic antibody_rpm columns
                rpm_columns = [f"{antibody}_rpm" for antibody in antibody_names if f"{antibody}_rpm" in gene_df.columns]

                # Extract columns to keep
                columns_to_keep = ['sample', 'gene_name', f"{gene}_FPKM", 'biological_condition'] + rpm_columns
                
                # Rename repeat TPM and histone columns, checking for column existence
                for repeat_name, data in repeat_mapping.items():
                    tpm_col = data["TPM"]
                    renamed_tpm_col = f"{repeat_name}_TPM"
                    if tpm_col in gene_df.columns:  # Check if TPM column exists
                        gene_df.rename(columns={tpm_col: renamed_tpm_col}, inplace=True)
                        columns_to_keep.append(renamed_tpm_col)
                    
                    for histone_col in data["histones"]:
                        if histone_col in gene_df.columns:  # Check if histone column exists
                            histone_type = histone_col.split("_")[0]  # Extract histone mark
                            renamed_histone_col = f"{histone_type}_{repeat_name}"
                            gene_df.rename(columns={histone_col: renamed_histone_col}, inplace=True)
                            columns_to_keep.append(renamed_histone_col)
                
                print("Columns to Keep:", columns_to_keep)
                
                # Drop all repeat_id columns
                gene_df.drop(columns=repeat_id_cols, inplace=True)
                
                # Keep only relevant columns
                gene_df = gene_df[columns_to_keep]

                # Define the output file path with the gene name
                output_path_pickle = os.path.join(output_dir, f"{gene}_bnm_input.pkl")
                output_path_csv = os.path.join(output_dir, f"{gene}_bnm_input.csv")
                
                # Save the DataFrame to a pickle file
                with open(output_path_pickle, "wb") as f:
                    pickle.dump(gene_df, f)
                
                gene_df.to_csv(output_path_csv, index=False)
                
        except Exception as e:
            print(f"Error processing file {file_path}: {e}")

    if __name__ == "__main__":
        process_file(input_file)
    EOF
    """
}