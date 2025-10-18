#!/bin/bash

# Description:
# This script merges individual methylation data files into a population-level file.
# Dependencies: datamash, awk

# Input parameters
pop=$1           # Population group
indir=<input_directory>         # Directory containing individual sample data
inlab=lab.txt                   # List of individual samples
outdir=<output_directory>       # Output directory for the merged result

# Merge individual files into a population-level file
cd ${outdir}

# Create symbolic links for each individual sample's data
for sample in $(cat ${inlab}); do 
    ln -s ${indir}/${sample}_cpg5.cpg ${sample}.cpg5;  # Link cpg5 file for each sample
done

popfile=${pop}.result
sampleNum=$(wc -l < ${inlab})  # Count the number of samples in the list

# Combine individual files and calculate population-level methylation data
cat ./*cpg5 | sort -k1,1V -k2,2n -k3,3n | datamash -g1,2,3 mean 4,5,6,8,9 count 7 -R 4 |
    awk -v num=${sampleNum} -v OFS="\t" '{$9=$9/num; print $0}' > ${popfile}
    echo "Merging completed. Results saved to: ${popfile}"
done
