#!/bin/bash

# ===========================
# SV calling & merging script
# Tools (versions used in paper):
#   cuteSV v1.0.13
#   Sniffles v1.0.12
#   SVIM v2.0.0
#   NanoVar v1.3.9
#   SURVIVOR v1.0.6
#   bcftools >=1.10
#
# Usage:
#   sv_call_merge.sh <BAM> <REF_FASTA> <SAMPLE_ID> <OUTDIR> [THREADS]
#
# Example:
#   sv_call_merge.sh sample.sorted.bam hg38.fa HG002 results 16
# ===========================

BAM="${1}"
REF="${2}"
SAMPLE="${3}"
OUTDIR="${4}"
THREADS="${5:-8}"

mkdir -p "${OUTDIR}"

# ---------------------------
# STEP 1. cuteSV
# - Input: sorted BAM + reference
# - Output: VCF to OUTDIR/cutesv/<sample>.cutesv.vcf
# ---------------------------
mkdir -p "${OUTDIR}/cutesv_tmp" "${OUTDIR}/cutesv"
cuteSV "${BAM}" "${REF}" \
        "${OUTDIR}/cutesv/${SAMPLE}.cutesv.vcf" "${OUTDIR}/cutesv_tmp" \
        --threads "${THREADS}" --genotype

# ---------------------------
# STEP 2. Sniffles (v1.x syntax uses -m for input)
# - Input: BAM
# - Output: OUTDIR/sniffles/<sample>.sniffles.vcf
# ---------------------------
mkdir -p "${OUTDIR}/sniffles"
sniffles -m "${BAM}" -v "${OUTDIR}/sniffles/${SAMPLE}.sniffles.vcf" -t "${THREADS}"

# ---------------------------
# STEP 3. SVIM
# - Mode: alignment (BAM + reference)
# - Output VCF appears at OUTDIR/svim/variants.vcf
# ---------------------------
mkdir -p "${OUTDIR}/svim"
svim alignment "${OUTDIR}/svim" "${BAM}" "${REF}" --threads "${THREADS}"
mv "${OUTDIR}/svim/variants.vcf" "${OUTDIR}/svim/${SAMPLE}.svim.vcf"

# ---------------------------
# STEP 4. NanoVar
# - Input can be BAM; also requires reference and a working dir
# - Produces ${sample}.nanovar.pass.vcf in working dir
# ---------------------------
mkdir -p "${OUTDIR}/nanovar"
nanovar -t "${THREADS}" "${BAM}" "${REF}" "${OUTDIR}/nanovar"
# Standardize filename
if [ -f "${OUTDIR}/nanovar/${SAMPLE}.nanovar.pass.vcf"  ]; then
    mv "${OUTDIR}/nanovar/${SAMPLE}.nanovar.pass.vcf" "${OUTDIR}/nanovar/${SAMPLE}.nanovar.vcf"
else
    # fallback: some versions may write NanoVar.pass.vcf
    mv "${OUTDIR}/nanovar/"*.nanovar.pass.vcf "${OUTDIR}/nanovar/${SAMPLE}.nanovar.vcf"
fi

# ---------------------------
# STEP 5. SURVIVOR merge per sample
# - Merge the four caller VCFs with typical settings:
#   max distance = 1000 bp
#   min caller support = 2
#   require same SVTYPE = yes
#   require same strand = yes
#   distance not scaled by SV size = no
#   minimum SV size considered by SURVIVOR = 50 bp
# ---------------------------
VCFLIST="${OUTDIR}/${SAMPLE}.vcf_list.txt"
cat > "${VCFLIST}" <<EOF
${OUTDIR}/cutesv/${SAMPLE}.cutesv.vcf
${OUTDIR}/sniffles/${SAMPLE}.sniffles.vcf
${OUTDIR}/svim/${SAMPLE}.svim.vcf
${OUTDIR}/nanovar/${SAMPLE}.nanovar.vcf
EOF

# Params: <list> <maxdist> <minsupport> <type> <strand> <use_est_dist> <minsize> <outvcf>
SURVIVOR merge "${VCFLIST}" 1000 2 1 1 0 50 "${OUTDIR}/${SAMPLE}.merged.vcf"

# ---------------------------
# STEP 6. bcftools filter to keep:
#   - FILTER == PASS
#   - SVLEN within [50, 100000] (use abs for deletions)
#   - SVTYPE in {DEL, INS, DUP, INV}
# Output: bgzipped + indexed VCF
# ---------------------------
bcftools view \
    -i 'FILTER="PASS" && (INFO/SVTYPE ~ "^(DEL|INS|DUP|INV)$") && (abs(INFO/SVLEN)>=50 && abs(INFO/SVLEN)<=100000)' \
    "${OUTDIR}/${SAMPLE}.merged.vcf" \
    -Oz -o "${OUTDIR}/${SAMPLE}.sv.pass.len50_100k.DEL_INS_DUP_INV.vcf.gz"

bcftools index -t "${OUTDIR}/${SAMPLE}.sv.pass.len50_100k.DEL_INS_DUP_INV.vcf.gz"

echo "Done."
echo "Final VCF: ${OUTDIR}/${SAMPLE}.sv.pass.len50_100k.DEL_INS_DUP_INV.vcf.gz"
