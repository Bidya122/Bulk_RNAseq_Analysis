#!/bin/bash

# Give path to your aligned reads (.bam files)
cd ~/bulk_RNA_analysis/fastq/hisat2_alignments

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 -a /home/mohantybidya39/tools/Homo_sapiens.GRCh38.115.gtf \
        -o ~/bulk_RNA_analysis/fastq/quants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "✅ Completed $bam in $runtime minutes."
    echo "------------------------------------"
done
