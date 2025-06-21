#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Rsubread)
  library(edgeR)
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

option_list <- list(
  make_option(c("--bam"), type = "character", help = "Input BAM file"),
  make_option(c("--gtf"), type = "character", help = "GTF annotation file"),
  make_option(c("--outdir"), type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

bam_file <- opt$bam
gtf_file <- opt$gtf
outdir <- opt$outdir

# FeatureCounts
counts <- featureCounts(
  files = bam_file,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  countMultiMappingReads = FALSE,
  isPairedEnd = FALSE,
  nthreads = 24
)

counts_df <- as.data.frame(counts$counts)
write.table(counts$counts, file = file.path(outdir, "counts.tsv"), quote = FALSE, col.names = TRUE, sep = "\t")

# FPKM
dge <- DGEList(counts = counts$counts, genes = counts$annotation[, c("GeneID", "Length")])
dge <- calcNormFactors(dge)
FPKM <- rpkm(dge)
write.table(FPKM, file = file.path(outdir, "FPKM.tsv"), quote = FALSE, col.names = TRUE, sep = "\t")

# GTF filtering
gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)
gtf_filtered <- gtf_df %>% select(1:5, 10:12)

# Read counts and FPKM again for annotation
counts <- read.table(file.path(outdir, "counts.tsv"))
FPKM <- read.table(file.path(outdir, "FPKM.tsv"))

counts$gene_id <- rownames(counts)
FPKM$gene_id <- rownames(FPKM)

counts <- counts %>% select(gene_id, everything())
FPKM <- FPKM %>% select(gene_id, everything())

counts_annotated <- inner_join(counts, gtf_filtered, by = "gene_id") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  select(gene_id, gene_name, gene_type, seqnames, start, end, width, strand, everything())

FPKM_annotated <- merge(FPKM, gtf_filtered, by = "gene_id") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  select(gene_id, gene_name, gene_type, seqnames, start, end, width, strand, everything())

write.csv(counts_annotated, file.path(outdir, "counts_annotated.csv"), row.names = FALSE)
write.csv(FPKM_annotated, file.path(outdir, "FPKM_annotated.csv"), row.names = FALSE)

# Pivot long format
fpkm_long <- FPKM_annotated %>%
  pivot_longer(cols = starts_with("SRR"), names_to = "sample", values_to = "FPKM") %>%
  mutate(sample = str_remove(sample, "\\.bam$"))

write.csv(fpkm_long, file.path(outdir, "genes_counttable.csv"), row.names = FALSE)
