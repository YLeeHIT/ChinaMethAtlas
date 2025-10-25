#!/bin/bash

# ===================================================
# Compare allele frequencies (AF) between two populations
# Identify SNPs with:
#   - MAF > 0.05 in both populations
#   - AF difference > 0.1
#   - Mark if minor alleles differ between populations
#
# Example:
#   bash compare_af_one_chr.sh north xizang chr1 /path/to/base_dir
# ===================================================

# -------- Input Parameters --------
pop1=$1       # Population 1 (e.g., north)
pop2=$2       # Population 2 (e.g., xizang)
chr=$3        # Chromosome (e.g., chr1)
base_dir=$4   # Base directory containing input/output folders

# -------- Directory Paths --------
vcf1="${base_dir}/filter/${pop1}/${chr}_${pop1}.SNP.filtered.vcf.gz"
vcf2="${base_dir}/filter/${pop2}/${chr}_${pop2}.SNP.filtered.vcf.gz"
out_dir="${base_dir}/diff_SNP/${pop1}_${pop2}"
mkdir -p "$out_dir"

# -------- Calculate Allele Frequencies --------
echo "Calculating allele frequencies..."
vcftools --gzvcf "$vcf1" --freq --out "${out_dir}/${chr}_${pop1}" >/dev/null
vcftools --gzvcf "$vcf2" --freq --out "${out_dir}/${chr}_${pop2}" >/dev/null

# -------- File Definitions --------
frq1="${out_dir}/${chr}_${pop1}.frq"
frq2="${out_dir}/${chr}_${pop2}.frq"
output="${out_dir}/${chr}.AF_diff.txt"

# -------- Threshold Settings --------
af_diff_threshold=0.1
maf_threshold=0.05

# -------- Compare AF Between Populations --------
awk -v af_diff="$af_diff_threshold" -v maf="$maf_threshold" '
FNR == 1 { next  }

# --- Load AF from Population 1 ---
NR == FNR {
    key = $1":"$2
    split($5, a1, ":"); split($6, a2, ":")
    allele1 = a1[1]; freq1 = a1[2] + 0
    allele2 = a2[1]; freq2 = a2[2] + 0

    if (freq1 <= freq2) {
        maf1 = freq1; minor1 = allele1; major1 = allele2                            
    } else {
        maf1 = freq2; minor1 = allele2; major1 = allele1                            
    }

    if (maf1 > maf) {
        maf_map[key] = maf1
        minor_map[key] = minor1
        major_map[key] = major1                                            
    }
    next
                        
}

# --- Compare with Population 2 ---
{
    key = $1":"$2
    if (!(key in maf_map)) next

    split($5, a1, ":"); split($6, a2, ":")
    allele1 = a1[1]; freq1 = a1[2] + 0
    allele2 = a2[1]; freq2 = a2[2] + 0

    if (freq1 <= freq2) {
        maf2 = freq1; minor2 = allele1; major2 = allele2                                    
    } else {
        maf2 = freq2; minor2 = allele2; major2 = allele1                                        
    }

    if (maf2 <= maf) next

    maf1 = maf_map[key]
    diff = (maf1 > maf2 ? maf1 - maf2 : maf2 - maf1)

    if (diff > af_diff) {
        status = (minor2 == minor_map[key]) ? "SAME" : "DIFF"
        geno1 = minor_map[key] "/" major_map[key]
        geno2 = minor2 "/" major2
        printf "%s\t%s\t%.4f\t%.4f\t%.4f\t%s\t%s\t%s\t%s\n", \
            $1, $2, maf1, maf2, diff, minor_map[key], minor2, status, \
            (status == "DIFF" ? geno1 " vs " geno2 : "-")                                                                            
    }
                                                
}' "$frq1" "$frq2" > "$output"

# -------- Completion Message --------
echo "AF comparison completed for ${chr} (${pop1} vs ${pop2})"
echo "Output saved to: ${output}"

