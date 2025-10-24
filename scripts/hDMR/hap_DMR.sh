#!/bin/bash

# ==========================================
# hDMR Detection Pipeline (Metilene)
# Steps:
#   1. Prepare haplotype input files
#   2. Identify DMRs and DMCs using Metilene
#   3. Filter DMRs based on statistical and length thresholds
# ==========================================

# -------- Input Parameters --------
sampleID=$1
default_indir=$2
threads=${3:-8}
HP1="${sampleID}_HP1_MethylFrequency.tsv"
HP2="${sampleID}_HP2_MethylFrequency.tsv"

# -------- Create Output Directory --------
mkdir -p "${default_indir}/DMR"
cd "${default_indir}/DMR"

# -------- Step 1: Prepare Input Files --------
# Convert haplotype methylation tables to BED format
cut -f1,2,3,6 "${default_indir}/${HP1}" | sed '1d' | sort -k1,1 -k2,2n > "${sampleID}_HP1.bed"
cut -f1,2,3,6 "${default_indir}/${HP2}" | sed '1d' | sort -k1,1 -k2,2n > "${sampleID}_HP2.bed"

# Build Metilene input file
metilene_input.pl \
    --in1 "${sampleID}_HP1.bed" \
    --in2 "${sampleID}_HP2.bed" \
    --h1 hap1 --h2 hap2 \
    --out "${sampleID}_hap1_hap2.bed"

echo "Step 1: Input preparation finished."

# -------- Step 2: Identify DMRs and DMCs --------
maxdist=1000     # Maximum distance between CpGs
mincpgs=5        # Minimum number of CpGs per region
minMethDiff=0.1  # Minimum methylation difference

# DMR detection
metilene -M ${maxdist} -m ${mincpgs} -d ${minMethDiff} -t ${threads} \
    -f 1 -a hap1 -b hap2 "${sampleID}_hap1_hap2.bed" \
    | sort -V -k1,1 -k2,2n > "${sampleID}_DMRs.txt"

echo "Step 2: DMR detection finished."

# -------- Step 3: Filter DMRs --------
# Keep regions with q < 0.05, length ≥ 100 bp
minlen=100

awk -v minlen=${minlen} 'BEGIN{
  print "chr\tstart\tstop\tq-value\tdelta\tnum\tpMWU\tp2D\tmeang1\tmeang2"
  
}{
    len=$3-$2+1;
    if($4<0.05 && len>=minlen){print $0}
        
}' "${sampleID}_DMRs.txt" > "${sampleID}_filter.DMRs.txt"

echo "Step 3: Filtering finished."

# -------- Done --------
echo "All steps completed successfully."

