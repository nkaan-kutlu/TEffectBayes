process BNM_INPUT_PREP_STEP2 {
    tag "merge_by_cellline"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'

    publishDir "${params.outdir}/Data_Integration/gene_repeat_intersect/cell_line/", mode: 'copy'
    
    input:
    path gene_repeat_intersect_transformed

    output:
    path "Data_Integration/gene_repeat_intersect/cell_line/*.csv", emit: cell_line_csvs
    path "Data_Integration/gene_repeat_intersect/cell_line/*.pkl", emit: cell_line_pickles

    script:
    """
    mkdir -p Data_Integration/gene_repeat_intersect/cell_line

    python <<EOF
    import os
    import pickle
    import pandas as pd
    from concurrent.futures import ThreadPoolExecutor

    # Define the input and output directories
    input_dir = "${gene_repeat_intersect_transformed}"
    output_dir = "Data_Integration/gene_repeat_intersect/cell_line"
    os.makedirs(output_dir, exist_ok=True)

    # Get list of all transformed files
    transformed_files = [f for f in os.listdir(input_dir) if f.endswith("intersect.csv")]

    # Get unique cell lines and their associated sample files
    cell_line_files = {}

    # Group files by cell_line based on the `cell_line` column in each CSV file
    for file in transformed_files:
        df = pd.read_csv(os.path.join(input_dir, file))
        cell_lines = df['cell_line'].unique()
        for cell_line in cell_lines:
            if cell_line not in cell_line_files:
                cell_line_files[cell_line] = []
            cell_line_files[cell_line].append(os.path.join(input_dir, file))

    def process_cell_line(cell_line, files):
        print(f"Processing {cell_line}...")
        cell_line_dfs = []
        for file in files:
            df = pd.read_csv(file)
            df_cell_line = df[df['cell_line'] == cell_line]
            cell_line_dfs.append(df_cell_line)
        final_df = pd.concat(cell_line_dfs, axis=0)
        repeat_columns = [col for col in final_df.columns if 'repeat' in col and ('start' in col or 'end' in col)]
        final_df[repeat_columns] = final_df[repeat_columns].fillna(0)
        for col in repeat_columns:
            final_df[col] = pd.to_numeric(final_df[col], errors='coerce', downcast='integer')
        final_df[repeat_columns] = final_df[repeat_columns].replace(0, pd.NA)
        return final_df

    if __name__ == "__main__":
        with ThreadPoolExecutor(max_workers=4) as executor:
            results = list(executor.map(process_cell_line, cell_line_files.keys(), cell_line_files.values()))

        for cell_line, final_df in zip(cell_line_files.keys(), results):
            pickle_path = os.path.join(output_dir, f"intersect_{cell_line}.pkl")
            csv_path = os.path.join(output_dir, f"intersect_{cell_line}.csv")
            final_df.to_pickle(pickle_path)
            final_df.to_csv(csv_path, index=False)
            print(f"Finished processing {cell_line}. Results saved to {pickle_path} and {csv_path}")
    EOF
    """
}