#!/bin/bash

# Population VCF Merging, Standardization, and Filtering Script
# This script performs the following steps:
# 1. Merge VCF files from all samples into a single file.
# 2. Standardize the merged VCF file by retaining only necessary information fields.
# 3. Split the merged VCF into different population files.
# 4. Apply specific filters, such as Hardy-Weinberg Equilibrium (HWE) and genotype filtering, to each population.
# 5. Further filter by heterozygosity and structural variant (SV) types.

# Requirements: jasmine, bcftools, vcftools

# Constants
input_list=$1            # Input list of VCF files
min_support=$2
pop=$3
merged_vcf="merged_population.vcf"      # Output merged VCF file
sorted_vcf="merged_population_sorted.vcf" # Sorted merged VCF
threads=8                               # Number of threads for processing
#input_list="inlist.txt"            # Input list of VCF files
#min_support=5                           # Minimum support for jasmine merging

# Step 1: Merge VCF files
# - The `inlist.txt` file should contain the paths to all VCF files to be merged.
echo "Step 1: Merging VCF files..."
jasmine file_list=${input_list} out_file=${merged_vcf} min_support=${min_support} threads=${threads} --output_genotypes > merge_out.log 2> merge_error.log

# Step 2: Annotate and standardize the merged VCF
# - Only retain necessary fields (SVTYPE, SVLEN, END, GT, DR, DV) to create a standardized VCF.
echo "Step 2: Annotating and standardizing the merged VCF..."
standardized_vcf="population_stand.vcf"
bcftools annotate -x ^INFO/SVTYPE,^INFO/SVLEN,^INFO/END,^FORMAT/GT,^FORMAT/DR,^FORMAT/DV -o ${standardized_vcf} --threads ${threads} ${merged_vcf}

# Step 3: Sort the standardized VCF
echo "Step 3: Sorting the standardized VCF..."
sorted_standardized_vcf="population_stand_sorted.vcf"
bcftools sort ${standardized_vcf} -o ${sorted_standardized_vcf}

# Step 4: Split by population
# - Example: Creating a sample list for the each population and extracting their VCF data.
echo "Step 4: Splitting VCF by population..."
pop_list="${pop}_samples.txt"
pop_vcf="${pop}_samples.vcf"
bcftools view -S ${pop_list} -o ${pop_vcf} ${sorted_standardized_vcf}

# Step 5: Filter population VCF by HWE and genotype
# - Filter based on HWE threshold and genotype quality.
pop="${pop}"  # Replace with actual population name (e.g., "north", "south", and "xizang")
AC=$(( $(wc -l < "${pop_list}") / 2 )) # Ensure high frequency of heterozygous DEL in the population
hwe_threshold=1e-6
geno_threshold=0.1
echo "Step 5: Filtering VCF by HWE and genotype for ${pop}..."

vcftools --vcf ${pop}_samples.vcf --max-missing 0.5 --hwe ${hwe_threshold} --recode --recode-INFO-all --out ${pop}_filtered_hwe

# Step 6: Split by specific SV types (INS, DEL, DUP, INV)
echo "Step 6: Splitting VCF by specific SV types for ${pop}..."
invcf="${pop}_filtered_hwe.recode.vcf"
bcftools view -i 'INFO/SVTYPE="INS"' ${invcf} -o ${pop}_filtered_INS.vcf
bcftools view -i 'INFO/SVTYPE="DEL"' ${invcf} -o ${pop}_filtered_DEL.vcf
bcftools view -i "INFO/AN-INFO/AC>=${AC}" ${pop}_filtered_DEL.vcf -o ${pop}_filtered_DEL_AC${AC}.vcf
bcftools view -i 'INFO/SVTYPE="DUP"' ${invcf} -o ${pop}_filtered_DUP.vcf
bcftools view -i 'INFO/SVTYPE="INV"' ${invcf} -o ${pop}_filtered_INV.vcf

echo "VCF merging, filtering, and SV type splitting completed successfully."
