#!/bin/bash

# This script performs DNA basecalling, alignment, and methylation calling using various tools.
# It includes conversion of FAST5 to POD5 format, followed by basecalling with Dorado,
# alignment to the GRch38 human reference genome, methylation inference with Remora, and final conversion
# Dependencies: pod5, dorado, remora, samtools

# Number of threads
THREADS=32

# Define CUDA devices and Remora device based on the system configuration
CUDA_DEVICES="cuda:0,1,2,3"  # Modify based on available GPU setup
REMORA_DEVICE=0  # Change as per your specific GPU device

# Set up file paths and directories
FAST5_INPUT_DIR="./input"
POD5_OUTPUT_FILE="converted.pod5"
POD5_DIR="./pod5s"
CALLS_BAM="calls.bam"
ALIGN_BAM="calls.align.bam"
REFERENCE_GENOME="hg38.fa"
REMORA_MODEL="train_results/model_best.pt"
REMORA_OUTPUT_BAM="can_infer.bam"
SORTED_BAM="methy.sort.bam"

# Step 1: Convert FAST5 to POD5 format
echo "Converting FAST5 files to POD5 format..."

# Option 1: Bulk conversion
pod5 convert fast5 "$FAST5_INPUT_DIR"/*.fast5 --output "$POD5_OUTPUT_FILE" -t "$THREADS"

# Option 2: One-to-one conversion for managing memory usage
# Uncomment the following line to perform one-to-one conversion
# pod5 convert fast5 "$FAST5_INPUT_DIR"/*.fast5 --output "$POD5_DIR" --one-to-one "$FAST5_INPUT_DIR" -t "$THREADS"

# Step 2: Perform basecalling using Dorado
echo "Running basecalling with Dorado..."
dorado basecaller dna_r9.4.1_e8_hac@v3.3/ "$POD5_DIR" --emit-moves -x "$CUDA_DEVICES" > "$CALLS_BAM"

# Step 3: Align reads to the reference genome
echo "Aligning basecalled reads to the reference genome..."
dorado aligner "$REFERENCE_GENOME" "$CALLS_BAM" > "$ALIGN_BAM" -t "$THREADS"

# Step 4: Methylation calling using Remora
echo "Running Remora methylation inference..."
remora infer from_pod5_and_bam "$POD5_OUTPUT_FILE" "$ALIGN_BAM" \
    --model "$REMORA_MODEL" \
    --out-bam "$REMORA_OUTPUT_BAM" \
    --device "$REMORA_DEVICE"

# Step 5: Sort and index BAM file
echo "Sorting and indexing the BAM file..."
samtools sort -o "$SORTED_BAM" "$REMORA_OUTPUT_BAM"
samtools index "$SORTED_BAM"
