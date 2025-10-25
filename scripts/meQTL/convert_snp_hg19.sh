#!/bin/bash

# ==========================================
# Convert genome coordinates using CrossMap
# Input file format: chr, pos, <annotations...>
# Output: converted coordinates with original annotations restored
# Example:
#   bash convert_coordinates.sh input.txt hg38ToHg19.over.chain.gz pop_name
# ==========================================

# -------- Input Parameters --------
input_file=$1     # Input annotation file
chain_file=$2     # Chain file for coordinate conversion
pop=$3            # Population or sample tag

# -------- Temporary and Output Files --------
tmp_bed="tmp_${pop}.bed"
mapped_bed="tmp_${pop}.mapped.bed"
output_file="${pop}.converted.hg19.txt"

# -------- Step 2: Build BED File --------
# Combine all annotation columns into a single string separated by "_"
awk -F'\t' -v OFS='\t' '{
    ann = $3;
    for (i = 4; i <= NF; i++) ann = ann "_" $i;
    print $1, $2 - 1, $2, ann            
}' "$input_file" > "$tmp_bed"

# -------- Step 3: Run CrossMap --------
CrossMap bed "$chain_file" "$tmp_bed" "$mapped_bed"

# -------- Step 4: Restore Annotations --------
# Output format: chr  pos  <original annotation fields>
awk -F'\t' -v OFS='\t' '{
    split($4, a, "_");
    if (a[6] == "DIFF") {
        print $1, $3, a[1], a[2], a[3], a[4], a[5], a[6], a[7], $6          
    } else {
        print $1, $3, a[1], a[2], a[3], a[4], a[5], a[6], a[7], "-"                
    }    
}' "$mapped_bed" > "$output_file"

# -------- Step 5: Cleanup --------
rm -f "$tmp_bed" "$mapped_bed"

# -------- Done --------
echo "Conversion completed. Output written to: $output_file"
