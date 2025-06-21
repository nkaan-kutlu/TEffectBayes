process locus_interval_prep {
    tag "locus_interval_prep"
    label 'process_high'

  conda "${moduleDir}/environment.yml"
  container 'nkaankutlu/pytools:latest'
      
  publishDir "${params.outdir}/TE_intervals", mode: 'copy'

    input:
    path te_counttable
    path samplesheet
    path outdir

    output:
    path "TE_intervals/*_TE_intervals.csv", emit: te_intervals

    script:
    """
    python3 -u << EOF
    import os
    import pandas as pd
    from tqdm import tqdm

    # Inputs
    count_file = "${te_counttable}"
    samplesheet_file = "${samplesheet}"
    outdir = os.path.join("${outdir}", "TE_intervals")
    os.makedirs(outdir, exist_ok=True)
    
    # Load count table in chunks
    chunk_size = 10000
    total_rows = sum(1 for _ in open(count_file)) - 1
    chunks = []
    with tqdm(total=total_rows, desc="Loading CSV") as pbar:
        for chunk in pd.read_csv(count_file, index_col=0, chunksize=chunk_size):
                chunks.append(chunk)
                        pbar.update(len(chunk))
    df_te = pd.concat(chunks)
    
    # Rename columns for PyRanges
    df_te.rename(columns={'seqnames': 'Chromosome', 'start': 'Start', 'end': 'End', 'strand': 'Strand'}, inplace=True)
    df_te['Start'] = df_te['Start'].astype(int)
    df_te['End'] = df_te['End'].astype(int)
    df_te['Strand'] = df_te['Strand'].apply(lambda x: x if x in ['+', '-', '.'] else '+')
    
    # Reverse Start-End for negative strand
    df_te.loc[df_te['Strand'] == '-', ['Start', 'End']] = df_te.loc[df_te['Strand'] == '-', ['End', 'Start']].values
    
    # Load samplesheet and create sample to cell_line mapping
    samples_df = pd.read_csv(samplesheet_file)
    sample_map = dict(zip(samples_df['sample'], samples_df['cell_line']))
    
    # Map cell_line
    df_te['cell_line'] = df_te['sample'].map(sample_map)
    
    # Output
    columns = ['cell_line', 'repeat_id', 'Chromosome', 'Strand', 'Start', 'End']
    for cell_line, group in df_te.groupby('cell_line'):
        output_df = group[columns].drop_duplicates()
    output_file = os.path.join(outdir, f"{cell_line}_TE_intervals.csv")
    output_df.to_csv(output_file, index=False)
    
    EOF
    """
}