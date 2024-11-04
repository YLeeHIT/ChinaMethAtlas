#!/bin/bash
# Structural Variant (SV) Filtering and Standardization for each sample
# This script filters and standardizes VCF files for specific structural variants (INS, DEL, DUP, INV).
# - Step 1: Filter VCF entries based on quality, depth, length, and type.
# - Step 2: Standardize the filtered VCF format by retaining essential fields and reheadering.

# Constants
inID=${sample}                                                       # Sample ID as input parameter
inSample="${inID}.cutesv.force_calling.genotype.vcf"                 # Input sample VCF file
outFilter="${inID}.Filter.vcf"                                       # Output file after filtering
outStandard="${inID}.Filter.Stand.vcf"                               # Final standardized VCF output file
threads=16

# Step 1: Filter VCF File
# Filters the VCF based on chromosome, variant type (INS, DEL, DUP, INV), quality, length, and depth.
echo "Processing ${inID} - Filtering Step"
cat <(bcftools view -h "${inSample}") <(bcftools view -H "${inSample}" |
    awk '$1~/^chr[1-9]$|^chr1[0-9]$|chr2[0-2]$/ {
                split($8, y, ";"); 
                type = substr(y[2], length(y[2])-2);
                len = int(substr(y[3], 7));
                end = int(substr(y[4], 5));
                start = int($2);
                quality = $7;
                if (quality == "PASS") {
                    if (type == "INS" || type == "DEL" || type == "DUP" || type == "INV") {
                        split($10, z, ":")
                        depth = z[2] + z[3];
                        if (len < 0) {len = -len}
                        if (depth >= 3) {
                            if (len >= 50 && len < 100000) {
                                print $0
                            }
                        }
                    }
                }
            }') > "${outFilter}"

# Step 2: Standardize VCF File
# Standardizes the VCF format by reheadering and retaining essential INFO and FORMAT fields.
echo "Standardizing ${inID}"
echo -e "NULL\t${inID}" > "${inID}.name"
bcftools reheader -s "${inID}.name" -o "tmp.${outStandard}" "${outFilter}"
bcftools annotate -x ^INFO/SVTYPE,^INFO/SVLEN,^INFO/END,^FORMAT/GT,^FORMAT/DR,^FORMAT/DV "tmp.${outStandard}" -o "${outStandard}" --threads ${threads}

# Cleanup temporary files
rm "tmp.${outStandard}" "${inID}.name"

echo "Filtering and standardization completed for ${inID}."
