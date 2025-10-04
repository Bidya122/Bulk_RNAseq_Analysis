install.packages("dplyr") # Load required libraries
library(dplyr) 
setwd("~/bulk_RNA_analysis/fastq/quants") # Set working directory to your quants folder
fc_files <- list.files(pattern = "_featurecounts.txt$") # List all featureCounts txt files (excluding .summary files)
# Function to read each file and keep only gene counts
read_fc <- function(file) {
  df <- read.delim(file, comment.char="#")  # featureCounts has comment lines starting with #
  # Usually first 5 columns are annotation, counts start at column 7 (adjust if needed)
  df_counts <- df[, c(1, 7)]  # gene_id and count
  colnames(df_counts)[2] <- gsub("_featurecounts.txt", "", file) # rename count column
  return(df_counts)
}
# Read and merge all files
fc_list <- lapply(fc_files, read_fc)
fc_merged <- Reduce(function(x, y) full_join(x, y, by="Geneid"), fc_list)
# Write merged counts to a file
write.table(fc_merged, "merged_featurecounts.txt", sep="\t", row.names=FALSE, quote=FALSE)
# Check first few rows
head(fc_merged)

