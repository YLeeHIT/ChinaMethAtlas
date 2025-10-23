#!/bin/bash

# ==========================================
# SNV detection (Clair3, fixed ONT mode) + filtering (bcftools)
# Tools:
#   - Clair3 v1.0.4
#   - bcftools v1.17
#
# Usage:
#   snv_ont_clair3.sh <BAM> <REF_FASTA> <CLAIR3_MODEL_DIR> <SAMPLE_ID> <OUTDIR> [THREADS]
#
# Example:
#   snv_ont_clair3.sh sample.sorted.bam hg38.fa /path/to/clair3_models/ont r10_441_HAC HG002 results 16
#
# Notes:
#   - Platform is FIXED to ONT (--platform ont).
#   - Final output: SNPs on autosomes, FILTER=PASS, DP>3, GQ>10.
# ==========================================

BAM="${1}"
REF="${2}"
MODEL_DIR="${3}"
MODEL_NAME="${4}"
SAMPLE="${5}"
OUTDIR="${6}"
THREADS="${7:-8}"

# ---------------------------
# STEP 1. Run Clair3 (ONT platform fixed)
# Output: <OUTDIR>/clair3/merge_output.vcf.gz
# ---------------------------
mkdir -p "${OUTDIR}/clair3"

bash run_clair3.sh \
    --bam_fn "${BAM}" \
    --ref_fn "${REF}" \
    --threads "${THREADS}" \
    --platform "ont" \
    --model_path "${MODEL_DIR}" \
    --output "${OUTDIR}/clair3"

INVCF="${OUTDIR}/clair3/merge_output.vcf.gz"

# ---------------------------
# STEP 2. bcftools filtering
# Keep SNPs only; FILTER=PASS; per-genotype DP>3 and GQ>10; autosomes only.
# ---------------------------
mkdir -p "${OUTDIR}/filter"

bcftools view \
    -r "${REGIONS}" \
    -i 'TYPE="snp" && FILTER="PASS" && (FMT/DP>3) && (FMT/GQ>10)' \
    "${INVCF}" \
    -Oz -o "${OUTDIR}/filter/${SAMPLE}.snv.pass.dp3.gq10.autosomes.vcf.gz"

bcftools index -t "${OUTDIR}/filter/${SAMPLE}.snv.pass.dp3.gq10.autosomes.vcf.gz"

echo "Done."
echo "Final VCF: ${OUTDIR}/filter/${SAMPLE}.snv.pass.dp3.gq10.autosomes.vcf.gz"
