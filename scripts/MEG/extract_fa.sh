#!/bin/bash

# Script Purpose:
# This script is designed to annotate INS (insertion) sequences by aligning them to MEGs (Mobile Element Genes) sequences.
# It constructs the population-level MEGs input file, extracts consensus sequences, performs alignment using minimap2,
# and filters results based on quality scores.

# Input parameters
inSample=$1             # Input sample name
inPop=$2                # Population name
threshold=16            # Number of threads for minimap2
ref="<path_to_super_TE_reference>"  # Path to the MEGs reference sequence (FASTA format)
rootDir="<path_to_INS_root_directory>"  # Root directory containing INS data
sample=$(echo -e "${inSample}" | cut -d"_" -f2)  # Extract sample name from input
outDir="<path_to_output_directory>/seqAnno"  # Output directory for annotation results

# Create output directory if it does not exist
test -d ${outDir} || mkdir -p ${outDir}
cd ${outDir}

# Population-level INS input and MEGs input file
popIns="<path_to_population_INS_file>"  # Population INS result file
popMEGs="${outDir}/${inPop}_MEGs.infile"  # Population MEGs input file

# Step 1: Construct the population MEGs input file
if [ ! -f ${popMEGs} ]; then
    awk 'NR>1{print $1"_"$2}' ${popIns} > ${popMEGs}
fi

# Step 2: Build consensus FASTA file from consensus reads
test -d ${sample} || mkdir -p ${sample}
cd ${sample}
outFile="${outDir}/${sample}/${sample}_fromCons.fa"  # Output FASTA file for consensus sequences

# Loop through each INS and extract consensus sequences
for i in $(cat ${popMEGs}); do
    echo "Processing INS: ${i}"
    inCons="${rootDir}/${inPop}/${inSample}/result/${i}/${sample}_${i}.cons"  # Consensus sequence file
    if [ -f ${inCons} ]; then
        echo "Found consensus sequence for INS: ${i}"
        echo ">${i}" >> ${outFile}
        awk 'NR==2' ${inCons} >> ${outFile}  # Append consensus sequence to FASTA
    else
        echo "WARNING: Consensus sequence for INS: ${i} does not exist"
    fi
done

# Step 3: Align the extracted consensus sequences to MEGs reference using minimap2
outSam="${outDir}/${sample}/${sample}_fromCons.sam"  # Output SAM file for alignment
minimap2 -t ${threshold} -ax map-ont ${ref} ${outFile} > ${outSam}

# Step 4: Filter alignment results based on quality scores
outMEGs="${outDir}/${sample}/${sample}_pos_MEGs.txt"  # Output file for annotated MEGs positions
samtools view ${outSam} | awk '$5 > 10 {print $1"\t"$3}' | sort -u -k1 > ${outMEGs}

# Completion message
echo "Annotation of INS to MEGs completed. Results are saved in ${outMEGs}"