#!/bin/bash

# Usage: bash normalize.sh <input_filename> <input_directory> <output_filename>
# Description: This script normalizes methylation values (body, up, down) to [0,1] scale using min-max scaling.

infile=$1
indir=$2
outfile=$3

awk -v OFS='\t' '
{
    n = 0

    # Collect positive values from body/up/down columns
    for (i = 10; i <= 12; i++) {
        if ($i > 0) values[++n] = $i
    }

    # Collect positive values from extended range (percentage-type values), scaled to [0,1]
    for (i = 16; i <= 25; i++) {
        if ($i > 0) values[++n] = $i / 100
    }

    # Handle case where no valid values are found
    if (n == 0) {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"-1","-1","-1"
        next
    }

    # Find min and max among collected values
    min = max = values[1]
    for (i = 2; i <= n; i++) {
        if (values[i] < min) min = values[i]
        if (values[i] > max) max = values[i]
    }

    # If all values are the same, assign "-1" to indicate invalid normalization
    if (max == min) {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"-1","-1","-1"
    } else {
        body = ($10 > 0) ? ($10 - min) / (max - min) : -1
        up   = ($11 > 0) ? ($11 - min) / (max - min) : -1
        down = ($12 > 0) ? ($12 - min) / (max - min) : -1
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,body,up,down
    }
}
' "${indir}/${infile}" > "${outfile}"


