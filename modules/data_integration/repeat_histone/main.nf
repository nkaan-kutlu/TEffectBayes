process repeat_histone_intersection {
    tag "merge_chip_rna"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pybedtools:latest'

    publishDir "${params.outdir}/Data_Integration/repeat_histone_intersect/", mode: 'copy'

    input:
    path samplesheet
    path te_intervals
    path histone_intervals
    
    output:
    path "Data_Integration/repeat_histone_intersect", emit: repeat_histone_intersect

    script:
    """
    mkdir -p Data_Integration/repeat_histone_intersect
    
    python3 -u << EOF
    import os
    import glob
    import pandas as pd
    import pybedtools
    
    # Load samplesheet and get unique cell lines
    samplesheet = pd.read_csv("${samplesheet}")
    cell_lines = samplesheet['cell_line'].unique()
    
    chip_dir = "${histone_intervals}"
    te_dir = "${te_intervals}"
    out_dir = "Data_Integration/repeat_histone_intersect"
    os.makedirs(out_dir, exist_ok=True)
    
    # Helper function to extract chip modification type
    def extract_chip_mod(chip_filename):
        base = os.path.basename(chip_filename)
        parts = base.split("_")
        if len(parts) >= 3:
            return parts[1].upper()
        else:
            return "UNKNOWN"
    
    # Loop over each cell line
    for cl in cell_lines:
        print(f"Processing cell line: {cl}")
        
        # Define TE file for the cell line
        te_file = os.path.join(te_dir, f"{cl}_TE_intervals.csv")
        if not os.path.exists(te_file):
            print(f"  TE file not found for {cl} in {te_dir}. Skipping cell line.")
            continue

        # Load TE file and reformat
        df_te = pd.read_csv(te_file)
        df_te['Start'] = df_te['Start'].astype(int)
        df_te['End'] = df_te['End'].astype(int)
        df_te['Start_corrected'] = df_te[['Start', 'End']].min(axis=1)
        df_te['End_corrected'] = df_te[['Start', 'End']].max(axis=1)
        df_te.drop(columns=['Start', 'End'], inplace=True)
        df_te.rename(columns={'Start_corrected': 'Start', 'End_corrected': 'End'}, inplace=True)
        df_te = df_te[['Chromosome', 'Start', 'End', 'Strand', 'repeat_id']]
        te_bt = pybedtools.BedTool.from_dataframe(df_te)

        # Find all ChIP files for this cell line
        chip_pattern = os.path.join(chip_dir, f"{cl}_*_intervals.csv")
        chip_files = glob.glob(chip_pattern)
        if not chip_files:
            print(f"  No chip files found for {cl}, skipping.")
            continue

        for chip_file in chip_files:
            chip_mod = extract_chip_mod(chip_file)
            print(f"  Intersecting with {os.path.basename(chip_file)} (modification: {chip_mod})")

            df_chip = pd.read_csv(chip_file)
            df_chip = df_chip[['Chr', 'Start', 'End', 'Strand', 'Peak_ID', 'Gene Name']]
            df_chip.rename(columns={'Chr': 'Chromosome'}, inplace=True)
            chip_bt = pybedtools.BedTool.from_dataframe(df_chip)

            intersect_bt = chip_bt.intersect(te_bt, f=0.50, r=True, wa=True, wb=True)

            if len(intersect_bt) == 0:
                print(f"    No intersections found in {os.path.basename(chip_file)}")
                continue

            df_intersect = intersect_bt.to_dataframe(names=[
                'chip_Chromosome', 'chip_Start', 'chip_End', 'chip_Strand', 'Peak_ID', 'Gene Name',
                'te_Chromosome', 'te_Start', 'te_End', 'te_Strand', 'repeat_id'
            ])

            df_intersect = df_intersect.drop(columns=['chip_Chromosome', 'chip_Start', 'chip_End', 'chip_Strand'])
            df_intersect = df_intersect.rename(columns={
            'Gene Name': 'gene_name',
            'te_Chromosome': 'Chromosome',
            'te_Start': 'Start',
            'te_End': 'End',
            'te_Strand': 'Strand'})

            df_out = pd.DataFrame({
            'cell_line': cl,
            'Peak_ID': df_intersect['Peak_ID'],
            'gene_name': df_intersect['gene_name'],
            'Chromosome': df_intersect['Chromosome'],
            'Start': df_intersect['Start'],
            'End': df_intersect['End'],
            'Strand': df_intersect['Strand'],
            'repeat_id': df_intersect['repeat_id']})

            df_out = df_out.drop_duplicates()
            out_filename = os.path.join(out_dir, f"{cl}_INTERSECT_{chip_mod}.csv")
            df_out.to_csv(out_filename, index=False)

    EOF
    """
}
