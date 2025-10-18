#!/bin/bash

# Nanopolish Methylation Processing Script
# This script processes methylation data from Oxford Nanopore sequencing using nanopolish.
# - Step 1: Merges all individual FASTQ files into a single file.
# - Step 2: Preprocesses the data by indexing FASTQ files.
# - Step 3: Aligns reads to the reference genome.
# - Step 4: Calls methylation using nanopolish.
# - Step 5: Calculates methylation frequency.

# Constants
reference="reference_hg38.fa"  # Reference genome file
summary="sequencing_summary.txt" # Sequencing summary file
cal_freq="calculate_methylation_frequency.py"  # Script for calculating methylation frequency
threads=16  # Number of threads for parallel processing

# Variables
ID="sample_id"  # Sample ID
fast5="fast5_directory"  # Directory containing FAST5 files
fastq="fastq_directory"  # Directory containing FASTQ files

# Step 1: Merge all FASTQ files into a single file
echo "Step 1: Merging individual FASTQ files into a single file"
cat ${fastq}/*.fastq > ${ID}.fastq

# Step 2: Data preprocessing (indexing)
echo "Step 2: Indexing data for nanopolish"
nanopolish index -d ${fast5}/ ${ID}.fastq
# Alternative: Include sequencing summary file if available
# nanopolish index -d ${fast5} -s ${summary} ${ID}.fastq

# Step 3: Align reads to the reference genome
echo "Step 3: Aligning reads to the reference genome"
minimap2 -a -x map-ont ${reference} ${ID}.fastq | samtools sort -T tmp -o ${ID}.sorted.bam
samtools index ${ID}.sorted.bam 

# Step 4: Call methylation
echo "Step 4: Calling methylation"
nanopolish call-methylation -t ${threads} -r ${ID}.fastq -b ${ID}.sorted.bam -g ${reference} > ${ID}_methylation_calls.tsv

# Step 5: Calculate methylation frequency
echo "Step 5: Calculating methylation frequency"
python3 ${cal_freq} ${ID}_methylation_calls.tsv > ${ID}_methylation_frequency.tsv
python3 ${cal_freq} ${ID}_methylation_calls.tsv -s > ${ID}_methylation_split_frequency.tsv
