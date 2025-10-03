#!/bin/bash
# fastqc_all_trimmed.sh
# Run FastQC on all trimmed fastq.gz files in the fastq folder

INPUT_DIR="$HOME/bulk_RNA_analysis/fastq"
OUTPUT_DIR="$HOME/bulk_RNA_analysis/fastq/trimmed_fastqc_results"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR_trimmed"

# Loop through all trimmed fastq.gz files
for file in "$INPUT_DIR"/*_trimmed.fastq.gz
do
    echo "Running FastQC on $file..."
    fastqc -o "$OUTPUT_DIR" "$file"
done

echo "✅ FastQC analysis completed. Results are in $OUTPUT_DIR"

