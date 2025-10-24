#!/bin/bash

# Description:
# Annotates INS (insertion) sequences by aligning sample-level consensus sequences to MEGs (Mobile Element Genes).
# Steps:
#   1. Generate a population-wide MEGs ID list.
#   2. Extract sample consensus sequences for all population INS.
#   3. Align these sequences to MEGs reference using minimap2.
#   4. Filter alignment results based on quality threshold.

# ==============================
# Input parameters and config
# ==============================

inSample=$1                     # Sample ID (e.g. SAM_XYZ)
inPop=$2                        # Population name (e.g. North, South, Tibet)
threads=${THREADS:-16}          # Number of threads for minimap2 (default: 16)

# Abstracted paths (can be overridden by environment variables)
REF_MEGS="${REF_MEGS:-../../data/super_TE.fa}"                                 # MEGs reference FASTA
INS_ROOT_DIR="${INS_ROOT_DIR:-../../Demo/for_MEG}"                                      # Root directory for INS results
POP_INS_FILE="${POP_INS_FILE:-${INS_ROOT_DIR}/${inPop}_reAlign.result}"    # Population-wide INS file
OUT_ROOT="${OUT_ROOT:-../../Demo/for_MEG}"                                            # Output directory

# Derived sample-specific variables
#sample=$(echo "${inSample}" | cut -d"_" -f2)
sample="${inSample}"
sampleOutDir="${OUT_ROOT}"
popMEGsFile="${OUT_ROOT}/${inPop}_MEGs.infile"
consensusFasta="${sampleOutDir}/${sample}_fromCons.fa"
alignmentSam="${sampleOutDir}/${sample}_fromCons.sam"
annotatedOutput="${sampleOutDir}/${sample}_pos_MEGs.txt"

# ==============================
# Setup and checks
# ==============================

# Create output directories if needed
mkdir -p "${sampleOutDir}"

# Check required input files
if [ ! -f "${POP_INS_FILE}" ]; then
    echo "Error: Population INS file not found at ${POP_INS_FILE}"
    exit 1
fi

if [ ! -f "${REF_MEGS}" ]; then
    echo "Error: MEGs reference file not found at ${REF_MEGS}"
    exit 1
fi


# ==============================
# Step 1: Build population-level MEGs ID list
# ==============================

if [ ! -f "${popMEGsFile}" ]; then
    echo "Building MEGs ID file for population: ${inPop}"
    awk 'NR > 1 {print $1"_"$2}' "${POP_INS_FILE}" > "${popMEGsFile}"
fi

# ==============================
# Step 2: Extract sample consensus sequences
# ==============================

echo "Extracting consensus sequences for sample: ${inSample}"
> "${consensusFasta}"  # Initialize output FASTA file

while read -r insID; do
    echo "Processing INS: ${insID}"
    consFile="${INS_ROOT_DIR}/${insID}/${sample}_${insID}.cons"

    if [ -f "${consFile}" ]; then
        echo ">${insID}" >> "${consensusFasta}"
        awk 'NR==2' "${consFile}" >> "${consensusFasta}"
    else
        echo "Warning: Consensus file not found for ${insID}"
    fi
done < "${popMEGsFile}"

# ==============================
# Step 3: Align consensus sequences to MEGs reference using minimap2
# ==============================

echo "Aligning consensus sequences using minimap2..."
minimap2 -t "${threads}" -ax map-ont "${REF_MEGS}" "${consensusFasta}" > "${alignmentSam}"

# ==============================
# Step 4: Filter alignments by quality and extract mapped MEGs
# ==============================

echo "Filtering alignments with MAPQ > 10..."
samtools view "${alignmentSam}" | \
        awk '$5 > 10 {print $1"\t"$3}' | sort -u -k1,1 > "${annotatedOutput}"

echo "INS to MEGs annotation completed."
echo "Results saved to: ${annotatedOutput}"


