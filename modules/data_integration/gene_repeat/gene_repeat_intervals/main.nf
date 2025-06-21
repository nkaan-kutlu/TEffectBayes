process gene_repeat_intervals {
    tag "merge_gene_repeat"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pybedtools:latest'

    publishDir "${params.outdir}/Data_Integration/gene_repeat_intersect/intervals/", mode: 'copy'

    input:
    path gene_histone_counttable
    path te_counttable

    output:
    path "Data_Integration/gene_repeat_intersect/intervals/*.bed", emit: gene_repeat_intervals

    script:
    """
    python <<EOF
    import os
    import pandas as pd
    import pybedtools
    from concurrent.futures import ThreadPoolExecutor

    output_dir = "Data_Integration/gene_repeat_intersect/intervals"
    os.makedirs(output_dir, exist_ok=True)

    def ensure_start_end(df, start_col="Start", end_col="End", strand_col="Strand"):
        mask = df[start_col] > df[end_col]
        df.loc[mask, [start_col, end_col]] = df.loc[mask, [end_col, start_col]].values
        return df

    # Load gene-histone data
    df_gene = pd.read_csv("${gene_histone_counttable}", index_col=0)
    df_gene = ensure_start_end(df_gene)

    # Load TE count table
    df_te = pd.read_csv("${te_counttable}", index_col=0)
    df_te.rename(columns={
        'seqnames': 'Chromosome',
        'start': 'Start',
        'end': 'End',
        'strand': 'Strand'
    }, inplace=True)
    df_te = ensure_start_end(df_te)

    def convert_to_bed(df, filename, cols):
        df[cols].to_csv(filename, sep="\\t", header=False, index=False)
        return pybedtools.BedTool(filename)

    def process_sample(sample):
        print(f"Processing {sample}...")

        promoters_sample = df_gene[df_gene["sample"] == sample]
        repeats_sample = df_te[df_te["sample"] == sample]

        if promoters_sample.empty and repeats_sample.empty:
            print(f"No data for {sample}, skipping...")
            return None

        if not promoters_sample.empty:
            promoters_sample = ensure_start_end(promoters_sample)
            promoter_bed_file = os.path.join(output_dir, f"{sample}_promoter.bed")
            convert_to_bed(promoters_sample, promoter_bed_file, [
                "Chromosome", "Start", "End",
                "gene_id", "gene_name", "gene_type",
                "sample", "cell_line", "FPKM"] + [col for col in promoters_sample.columns if col.endswith('_rpm')]
            )
            print(f"Saved: {promoter_bed_file}")

        if not repeats_sample.empty:
            repeats_sample = ensure_start_end(repeats_sample)
            repeat_bed_file =  os.path.join(output_dir, f"{sample}_repeat.bed")
            convert_to_bed(repeats_sample, repeat_bed_file, [
                "Chromosome", "Start", "End",
                "repeat_id", "TPM", "sample"
            ])
            print(f"Saved: {repeat_bed_file}")

        return True

    samples = df_gene["sample"].unique()
    print(f"Processing {len(samples)} samples...", flush=True)

    with ThreadPoolExecutor(max_workers=min(4, len(samples))) as executor:
        executor.map(process_sample, samples)

    EOF
    """
}
