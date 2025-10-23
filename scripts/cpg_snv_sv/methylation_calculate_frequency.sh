#!/bin/bash
# =========================================================
# CpG processing pipeline (simple version)
# Inputs:
#   1) SORTED_BAM : coordinate-sorted modBAM
#   2) HAP_HDMR   : hDMR regions (BED)
#   3) DEPTH      : average sequencing depth (integer)
#
# Steps:
#   step1) samtools: keep chr1–chr22, primary alignments, MAPQ>=5
#   step2) modkit:   make CpG bed (pileup + extract)
#   step3) filter:   keep CpGs with coverage >= (DEPTH/3)
#   step4) bedtools: subtract hDMR regions
#
# Notes:
#   - Requires: samtools, modkit, bedtools
#   - If your chromosomes are named "1..22" (no 'chr'), update the region list.
# =========================================================

set -euo pipefail

SORTED_BAM="$1"
HAP_HDMR="$2"
DEPTH="$3"
threads="$4"
REFERENCE_GENOME="hg38.fa"
MOTIF_BED="hg38.motif"
PILEUP_LOG="pileup.log"

CUTOFF=$(( (DEPTH + 2) / 3  ))

# Output prefixes
BASE_PREFIX="${SORTED_BAM%.bam}"
FILT_BAM="${BASE_PREFIX}.chr1_22.primary.mapq5.bam"
PILEUP_BED="${BASE_PREFIX}.pileup.bed"                # CpG bed from modkit pileup
EXTRACT_BED="${BASE_PREFIX}.extract.bed"              # CpG bed from modkit extract
PILEUP_BED_KEEP="${BASE_PREFIX}.pileup.depth${CUTOFF}.bed"
EXTRACT_TSV="${BASE_PREFIX}.extract.tsv"
PILEUP_BED_FINAL="${BASE_PREFIX}.pileup.depth${CUTOFF}.no_hDMR.bed"

# ---------------------------
# step1) samtools filtering
# ---------------------------
# -F 0x904 removes: unmapped(0x4) + secondary(0x100) + supplementary(0x800)
samtools view -b -q 5 -F 0x904 "$SORTED_BAM" \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 > "$FILT_BAM"

samtools index "$FILT_BAM"

# ---------------------------
# step2) modkit → CpG bed
# ---------------------------
echo "Generating methylation pileup and saving to BED format..."
modkit pileup "$FILT_BAM" "$PILEUP_BED" \
    --cpg --ref "$REFERENCE_GENOME" \
    --combine-strands \
    --log-filepath "$PILEUP_LOG" \
    -t "$THREADS"

echo "Generating motif BED file for CpG sites..."
modkit motif-bed "$REFERENCE_GENOME" CG 0 > "$MOTIF_BED"

echo "Extracting methylation data into TSV format..."
modkit extract "$FILT_BAM" "$EXTRACT_TSV" 
        --ref "$REFERENCE_GENOME" -t "$THREADS" --include "$MOTIF_BED"

# ---------------------------
# step3) filter by coverage ≥ DEPTH/3
# ---------------------------
awk -v c="$CUTOFF" 'BEGIN{OFS="\t"} $10 >= c{print $1,$2,$3,$10,$11}' "$PILEUP_BED"  > "$PILEUP_BED_KEEP"

# ---------------------------
# step4) remove hDMR regions
# ---------------------------
bedtools subtract -a "$PILEUP_BED_KEEP"  -b "$HAP_HDMR" > "$PILEUP_BED_FINAL"

# Done
echo "Filtered BAM:              $FILT_BAM"
echo "Pileup bed (kept):         $PILEUP_BED_KEEP"
echo "Extract bed (kept):        $EXTRACT_TSV"
echo "Pileup bed (no hDMR):      $PILEUP_BED_FINAL"
