#!/bin/bash

# ==========================================
# SNP Filtering Script (using VCFtools)
# Filters SNPs by MAF, HWE, and variant type
# ==========================================

# -------- Input Parameters --------
input_vcf=$1      # Input VCF file (.vcf.gz)
output_dir=$2     # Output directory

# -------- Create Output Directory --------
mkdir -p "${output_dir}"

# -------- Output Prefix --------
prefix=$(basename "${input_vcf%.vcf.gz}")
out_prefix="${output_dir}/${prefix}.SNP.filtered"

# -------- Step 1: Retain Only SNPs --------
vcftools --gzvcf "${input_vcf}" \
         --remove-indels \
         --maf 0.05 \
         --hwe 0.000001 \
         --recode --recode-INFO-all \
         --out "${out_prefix}"

# -------- Step 2: Compress and Index --------
bgzip -f "${out_prefix}.recode.vcf"
tabix -f -p vcf "${out_prefix}.recode.vcf.gz"

echo "SNP filtering completed for ${input_vcf}"
echo "Output: ${out_prefix}.recode.vcf.gz"
