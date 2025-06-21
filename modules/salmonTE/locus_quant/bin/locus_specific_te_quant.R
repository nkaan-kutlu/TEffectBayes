#!/usr/bin/env Rscript

library(optparse)
library(stringr)
library(tidyr)
library(dplyr)
library(purrr)
library(rtracklayer)

option_list <- list(
  make_option(c("--quant_dirs"), type = "character", help = "Space-separated list of quant.sf parent directories", action = "store"),
  make_option(c("--repeat_gtf"), type = "character", help = "Repeat annotation GTF file"),
  make_option(c("--outpath"), type = "character", help = "Output path for final CSV")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load annotation file
te_gtf <- rtracklayer::import.gff(opt$repeat_gtf)
te_df <- as.data.frame(te_gtf)
te_df$repeat_id <- as.character(te_df$repeat_id)

# List all quant.sf files from directories
quant_dirs <- unlist(strsplit(opt$quant_dirs, " "))
quant_files <- file.path(quant_dirs, "quant.sf")
samples <- basename(quant_dirs)

# Process each file
te_count_data <- list()
for (i in seq_along(quant_files)) {
  sample <- samples[i]
  file_path <- quant_files[i]
  sample_file <- read.table(file_path, header = TRUE)
  colnames(sample_file)[1] <- "repeat_id"
  sample_file <- sample_file %>%
    mutate(sample = sample) %>%
    select(repeat_id, TPM, sample)
  sample_file$repeat_id <- as.character(sample_file$repeat_id)
  sample_file <- merge(sample_file, te_df, by = "repeat_id")
  repeat_count <- sample_file %>%
    group_by(repeat_id) %>%
    summarise(total_length = sum(width), .groups = 'drop')
  sample_file <- sample_file %>%
    left_join(repeat_count, by = "repeat_id") %>%
    mutate(TPM = TPM * (width / total_length)) %>%
    select(-total_length)
  te_count_data[[sample]] <- sample_file
}

final_data <- bind_rows(te_count_data)
dir.create(opt$outpath, showWarnings = FALSE, recursive = TRUE)
write.csv(final_data, file.path(opt$outpath, "ls_TE_counttable.csv"), row.names = FALSE)
