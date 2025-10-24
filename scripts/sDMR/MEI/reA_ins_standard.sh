#!/bin/bash

# Usage: bash this_script.sh <sample_name> <population> <stepFreq> <stepLen>
# Description: Extracts methylation signals around structural variants (INS),
# computes summary statistics in step-wise windows, and performs normalization.

# Input parameters
sample=$1
pop=$2
stepFreq=$3
stepLen=$4
stepRange=2000

# Extract sample ID (e.g., from "sam_001" get "001")
#ID=$(echo "$sample" | cut -d"_" -f2)
ID="${sample}"

# Configurable paths (abstracted via environment variables or defaults)
GENOME_POS="${GENOME_POS:-../../data/genome.pos}"
BED_DIR="${BED_DIR:-../../Demo/for_MEG}"
PROJECT_ROOT="${PROJECT_ROOT:-../../Demo/for_MEG}"
NORMALIZE_SH="${NORMALIZE_SH:-./re_normalize.sh}"

# Derived paths
NO_IMPRINT_BED="${BED_DIR}/${ID}_noIm.bed"
INDIR="${PROJECT_ROOT}"
INFILE="${INDIR}/${sample}_reAlign.result"
INFILE_AUTO="${INDIR}/${sample}_reAlign_chrauto.result"
POSFILE="${INDIR}/${sample}_source.pos"

# Create output directory if not exists
test -d "$INDIR" || mkdir -p "$INDIR"

# Step 1: Filter autosomal chromosomes (chr1 to chr22)
awk 'NR>1 && $7 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "$INFILE" > "$INFILE_AUTO"

# Step 2: Extract relevant SV position columns
awk -v OFS='\t' '{print $7,$8,$9,$10,$11,$12}' "$INFILE_AUTO" | \
    grep -E '^chr([1-9]|1[0-9]|2[0-2])\b' > "$POSFILE"

# Step 3: Calculate methylation statistics in left and right windows
for i in $(seq 1 "$stepLen"); do
    j=$((i - 1))
    stepS=$(echo "$j * $stepRange + $stepFreq" | bc)
    stepE=$(echo "$i * $stepRange + $stepFreq" | bc)

    outLeft="${INDIR}/left_${i}_${stepS}_${stepE}.bed"
    outRight="${INDIR}/right_${i}_${stepS}_${stepE}.bed"

    # Left window: upstream of insertion site
    awk -v st=$stepS -v ed=$stepE -v OFS='\t' '{
        start = $2 - ed; end = $2 - st;
        if (start <= 0) start = 1;
        print $1, start, end, "left" NR
    }' "$POSFILE" | \
    bedtools intersect -a stdin -b "$NO_IMPRINT_BED" -wa -wb -loj | \
    sort -k4,4V -k1,1V -k2,2n -k3,3n | \
    datamash -g 4 count 8 mean 9 -R 2 | \
    awk '{if($2<5){value=-1}else{value=$3}; print $1"\t"$2"\t"value}' > "$outLeft"

    # Right window: downstream of insertion site
    awk -v st=$stepS -v ed=$stepE -v OFS='\t' '{
        start = $3 + st; end = $3 + ed;
        print $1, start, end, "right" NR
    }' "$POSFILE" | \
    bedtools intersect -a stdin -b "$NO_IMPRINT_BED" -wa -wb -loj | \
    sort -k4,4V -k1,1V -k2,2n -k3,3n | \
    datamash -g 4 count 8 mean 9 -R 2 | \
    awk '{if($2<5){value=-1}else{value=$3}; print $1"\t"$2"\t"value}' > "$outRight"
done

# Step 4: Combine left and right window results into range matrix
paste ${INDIR}/left_* ${INDIR}/right_* | \
    awk -v OFS='\t' '{print $3,$6,$9,$12,$15,$18,$21,$24,$27,$30}' \
    > "${INDIR}/${sample}_range.txt"

# Step 5: Merge original SV info with window stats
paste "$INFILE_AUTO" "${INDIR}/${sample}_range.txt" > "${INDIR}/${sample}_merge.txt"

# Step 6: Normalize the result using external normalization script
norIn="${INDIR}/${sample}_merge.txt"
norDir="$INDIR"
norOut="${INDIR}/${sample}_reAlign_new.result"
bash "$NORMALIZE_SH" "$norIn" "$norDir" "$norOut"

echo "Finished: Final result written to $norOut"


