#!/bin/bash

# DMR Calculation and Filtering Script
# This script calculates and filters DMRs and DMCs between two haplotypes (HP1 and HP2) using Metilene.
# - Step 1: Construct input files for Metilene from HP1 and HP2 methylation frequency data.
# - Step 2: Calculate DMRs and DMCs using Metilene.
# - Step 3: Filter DMRs based on length and q-value.

# Default Parameters
threads=16  # Number of threads for parallel processing
maxdist=1000  # Maximum distance between CpGs for DMRs
mincpgs=5  # Minimum number of CpGs in a DMR
minMethDiff=0.3  # Minimum methylation difference for DMRs and DMCs
minlen=100  # Minimum length for filtered DMRs

# Define Input and Output Files
sampleID="sample_name"  # Sample ID for output naming
HP1="haplotype1_methylation.tsv"  # Input methylation frequency file for HP1
HP2="haplotype2_methylation.tsv"  # Input methylation frequency file for HP2
inHP1="${sampleID}_HP1.bed"
inHP2="${sampleID}_HP2.bed"
OneOut="${sampleID}_hap1_hap2.bed"

# Step 1: Construct Metilene input file
echo "Step 1: Preparing Metilene input files..."

# Extract relevant columns and sort data for HP1 and HP2
cut -f 1,2,3,6 ${HP1} | sed '1d' | sort -k1,1 -k2,2n > ${inHP1}
cut -f 1,2,3,6 ${HP2} | sed '1d' | sort -k1,1 -k2,2n > ${inHP2}

# Create Metilene-compatible input file
metilene_input.pl --in1 ${inHP1} --in2 ${inHP2} --h1 hap1 --h2 hap2 --out ${OneOut}
echo "Step 1: Metilene input preparation completed."

# Step 2: Calculate DMRs and DMCs
echo "Step 2: Calculating DMRs and DMCs..."

# Define output files for DMRs and DMCs
DMRs_output="${sampleID}_DMRs.txt"
DMCs_output="${sampleID}_DMCs.txt"

# Calculate DMRs (Differentially Methylated Regions)
metilene -M ${maxdist} -m ${mincpgs} -d ${minMethDiff} -t ${threads} -f 1 \
        -a hap1 -b hap2 ${OneOut} | sort -V -k1,1 -k2,2n > ${DMRs_output}
echo "Step 2.1: DMR calculation completed."

# Calculate DMCs (Differentially Methylated CpGs)
metilene -d ${minMethDiff} -t ${threads} -f 3 -a hap1 -b hap2 ${OneOut} | \
        sort -V -k1,1 -k2,2n > ${DMCs_output}
echo "Step 2.2: DMC calculation completed."

# Step 3: Filter DMRs based on length and q-value
echo "Step 3: Filtering DMRs..."

# Define filtered DMR output file
filtered_DMRs="${sampleID}_filtered_DMRs.txt"

# Filter DMRs by minimum length and q-value
awk -v inlen=${minlen} 'BEGIN{print "chr\tstart\tstop\tq-value\tdelta\tnum\tpMWU\tp2D\tmeang1\tmeang2"} \
        {len=$3-$2+1; if($4<0.05 && len>inlen){print $0}}' ${DMRs_output} > ${filtered_DMRs}

echo "Step 3: DMR filtering completed."
 echo "DMR and DMC analysis, and filtering completed successfully."
