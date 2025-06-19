#!/bin/bash

#!/bin/bash
# Usage: bash calc_region_methylation.sh <chr> <start> <end> <noImBed>
# Description: Calculate average methylation level for a given genomic region.

chr=$1         # Chromosome name (e.g. chr1)
start=$2       # Region start (0-based)
end=$3         # Region end (exclusive)
noImBed=$4     # BED file with methylation values (5th column)

# Optional genome annotation path (not used here, but reserved)
GENOME_POS="${GENOME_POS:-../../genome.pos}"

# Input validation
if [[ -z "$chr" || -z "$start" || -z "$end" || -z "$noImBed" ]]; then
    echo "Error: Missing arguments."
    echo "Usage: bash $0 <chr> <start> <end> <noImBed>"
    exit 1
fi

# Compute average methylation value in the region
echo -e "${chr}\t${start}\t${end}" | \
    bedtools intersect -a "$noImBed" -b - | \
    datamash mean 5 -R 2

