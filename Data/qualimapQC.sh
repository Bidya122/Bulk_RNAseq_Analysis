#!/bin/bash
# Activate conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base   # or the environment where Qualimap is installed


mkdir -p ~/bulk_RNA_analysis/fastq/qualimap

for bam in ~/bulk_RNA_analysis/fastq/hisat2_alignments/*.bam
do
    sample=$(basename "$bam" .bam)
    echo "Running Qualimap for $sample ..."

    qualimap rnaseq \
        -bam "$bam" \
        -gtf ~/tools/Homo_sapiens.GRCh38.115.gtf \
        -outdir ~/bulk_RNA_analysis/fastq/qualimap/"$sample" \
        --java-mem-size=12G

    echo "Finished $sample"
done

