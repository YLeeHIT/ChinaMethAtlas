#!/bin/bash

# This script extracts insertion-specific BAM files based on an input BED file.

# Input arguments
sample=$1        # Sample ID (e.g., HG002)
inVcf=$2        # BED file containing INS positions
inBam=$3         # Original BAM file
#insBed=$2
# Create INS Bed with chr, start, end, INFO
insBed=${sample}_invcf.bed
#bcftools view -H ${inVcf} |awk '{print $1,$2,$2+1,$10}' OFS="\t" |head -n 10 > ${insBed}
bcftools view -H ${inVcf} |awk '{print $1,$2,$2+1,$10}' OFS="\t" > ${insBed}

# Output BAM file containing all reads overlapping insertion sites
inInsBam=${sample}_INS.bam

# Extract all reads from input BAM that overlap INS regions
samtools view -L ${insBed} ${inBam} -o ${inInsBam} -@ 16
samtools index ${inInsBam}

# Create subdirectories to store intermediate BED and BAM files
test -d bed || mkdir bed
test -d bam || mkdir bam
test -d out || mkdir out

# Iterate through each line of the input INS BED file
while read line
do
    # Construct a unique region name using the first two columns: chr and start
    lineName=$(echo "$line" | cut -f 1,2 --output-delimiter=_)
    echo -e "Processing region: $lineName"

    # Save the current line as a temporary region-specific BED file
    echo -e "${line}" > bed/${lineName}.txt

    # Extract reads overlapping this region from the pre-filtered INS BAM
    samtools view -L bed/${lineName}.txt ${inInsBam} -o bam/${sample}_${lineName}.bam -@ 16
    
    # Build consensus sequence and input ins meth
    bash meth_INS.sh bam/${sample}_${lineName}.bam bed/${lineName}.txt out/${lineName}

done < ${insBed}

