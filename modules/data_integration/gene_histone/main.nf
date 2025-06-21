process gene_histone_intersection {
    tag "merge_chip_rna"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'
    
    publishDir "${params.outdir}/Data_Integration", mode: 'copy'

    input:
    path genes_counttable         
    path histone_counttable         
    path samplesheet
    path genelist optional true       

    output:
    path "Data_Integration/gene_histone_counttable.csv", emit: gene_histone_counttable

    environment:
    interval = params.interval ?: 5000
    genelist_path = "${params.genelist ?: ''}"

    script:
    """
    python <<EOF
    import os
    import pandas as pd

    output_dir = "Data_Integration"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "gene_histone_counttable.csv")

    # Load RNA-seq gene counts
    df_gene = pd.read_csv("${genes_counttable}", index_col=0)
    df_gene = df_gene.rename(columns={"strand": "Strand"})

    # Promoter region intervals
    df_gene['new_start'] = df_gene.apply(lambda x: max(1, x['start'] - interval) if x['Strand'] == '+' else x['end'] + interval, axis=1)
    df_gene['new_end'] = df_gene.apply(lambda x: x['start'] - 1 if x['Strand'] == '+' else x['end'] + 1, axis=1)
    df_gene = df_gene.drop(columns=['start', 'end'])

    df_gene.rename(columns={'seqnames': 'Chromosome', 'new_start': 'Start', 'new_end': 'End'}, inplace=True)
    df_gene['Start'] = df_gene['Start'].astype(int)
    df_gene['End'] = df_gene['End'].astype(int)
    df_gene['Strand'] = df_gene['Strand'].apply(lambda x: x if x in {"+", "-", "."} else ".").astype("category")

    # Load samplesheet to extract sample-cell_line mapping
    df_samplesheet = pd.read_csv("${samplesheet}")
    sample_map = df_samplesheet.set_index("sample")["cell_line"].to_dict()

    # Map sample to cell line
    df_gene["cell_line"] = df_gene["sample"].map(sample_map)

    # Read ChIP tables
    chip_dir = os.path.dirname("${histone_counttable}")
    chip_files = [os.path.join(chip_dir, f) for f in os.listdir(chip_dir) if f.endswith("_counttable.csv")]

    final_df = df_gene.copy()
    rpm_columns = []

    for chip_file in chip_files:
        df_chip = pd.read_csv(chip_file)
        df_chip = df_chip[['Gene Name', 'cell_line', 'chip_rpm']]
        df_chip['cell_line'] = df_chip['cell_line'].str.split('_').str[0]
        df_chip.dropna(how='all', inplace=True)

        histone_mod = os.path.basename(chip_file).split("_")[0]
        col_name = f"{histone_mod}_rpm"
        rpm_columns.append(col_name)

        df_chip.rename(columns={'Gene Name': 'gene_name', 'chip_rpm': col_name}, inplace=True)

        final_df = pd.merge(final_df, df_chip, how='left', 
                            left_on=['cell_line', 'gene_name'], 
                            right_on=['cell_line', 'gene_name'])

        final_df.drop_duplicates(subset=['sample', 'cell_line', 'gene_name'], inplace=True)

    # Final column order
    base_cols = [
        'gene_id', 'sample', 'FPKM', 'Chromosome',
        'width', 'Strand', 'gene_type', 'gene_name',
        'Start', 'End', 'cell_line'
    ]
    all_cols = base_cols + rpm_columns
    final_df = final_df[all_cols]
       
    # Load genelist if provided
    genelist_file = "${genelist_path}"
    
    if genelist_file and os.path.isfile(genelist_file):
        with open(genelist_file, 'r') as f:
            gene_list = [line.strip() for line in f if line.strip()]
        print(f"Filtering with gene list of {len(gene_list)} genes from: {genelist_file}")
        filtered_df = final_df[final_df['gene_name'].isin(gene_list)]
        filtered_df.to_csv(output_path, index=True)
    else:
        print("No gene list provided or file not found, writing full dataframe.")
        final_df.to_csv(output_path, index=True)
    EOF
    """
}