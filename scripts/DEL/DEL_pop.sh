#!/bin/bash

# DEL Methylation Level Calculation and Normalization Script
# This script performs a multi-step process to calculate methylation levels in DEL regions:
# - Step 0: Extract 0/1 DELs from VCF.
# - Step 1: Construct DEL body methylation data.
# - Step 2: Generate upstream and downstream regions for DEL.
# - Step 3: Merge body, upstream, and downstream methylation data.
# - Step 4: Normalize methylation levels over a specified range.
# - Step 5: merge all cpg files

# Define input parameters and file paths
sample="sample_name"           # Sample identifier
population="population_name"   # Population name
indir="input_directory"        # Input directory for population-specific DEL data
infile="input_pop.vcf"         # Input VCF file for DEL information
noImBed="input_meth.bed"       # BED file of CpGs
genome="genome_position.txt"   # Genome file for bedtools flank operations
threads=16                     # CPU threads for parallel processing

# Step 0: Extract 0/1 and 1/1 DELs
echo "Step 0: Extracting high-frequency 0/1 DELs"
outfile1="${population}_del_01.vcf"
outfile2="${population}_del_11.vcf"

bcftools view -h ${infile} > ${outfile1}
bcftools view -H ${infile} | grep "0\/1" >> ${outfile1}
bcftools view -h ${infile} > ${outfile2}
bcftools view -H ${infile} | grep -v "0\/1" >> ${outfile2}

# Step 1: Construct DEL body methylation data
echo "Step 1: Constructing DEL body methylation data"
invcf="${outfile1}"
bcftools view -s ${sample} ${invcf} -o "${sample}_01.vcf"
bcftools view -H "${sample}_raw_01.vcf" | grep "0\/1" | \
        awk 'BEGIN{OFS="\t"}{split($8,x,";"); split(x[2],y,"="); len=(y[2]<0) ? -y[2] : y[2]; split(x[3],z,"="); end=z[2];
        print $1,$2,end,len}' > "${sample}_01.bed"

# Intersect with non-imprint regions and calculate CpG count and mean
inDel="${sample}_raw_01.bed"
bedtools intersect -a ${noImBed} -b ${inDel} -wa -wb | sort -k6,6V -k7,7n -k9,9V | \
        datamash -g 6,7,8,9 count 4 mean 5 | \
        datamash -g 1,2,3 absmax 5,6 | awk '$4>=5' > "${sample}_01_CpG_5.seg"
cut -f 1,2,3 "${sample}_01_CpG_5.seg" > "${sample}_01_CpG_5.pos"

# Step 2: Construct upstream and downstream methylation data
echo "Step 2: Generating upstream and downstream DEL regions"
bodyPos="${sample}_01_CpG_5.pos"

# Upstream (2kb)
bedtools flank -i ${bodyPos} -g ${genome} -l 2000 -r 0 | \
        awk '{print $0"\tUp"NR}' | \
        bedtools intersect -a stdin -b ${noImBed} -wa -wb -loj | \
        sort -k4,4V -k1,1V -k2,2n -k3,3n | \
        datamash -g 4 count 8 mean 9 > "${sample}_up.cpg"

# Downstream (2kb)
bedtools flank -i ${bodyPos} -g ${genome} -l 0 -r 2000 | \
        awk '{print $0"\tDown"NR}' | \
        bedtools intersect -a stdin -b ${noImBed} -wa -wb -loj | \
        sort -k4,4V -k1,1V -k2,2n -k3,3n | \
        datamash -g 4 count 8 mean 9 > "${sample}_down.cpg"

# Step 3: Merge body, upstream, and downstream methylation data
echo "Step 3: Merging body, upstream, and downstream methylation data"
paste "${sample}_01_CpG_5.seg" "${sample}_up.cpg" "${sample}_down.cpg" | \
        awk -v name=${sample} 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,$8,$11,$4,$7,$10,name}' > "${sample}_all.cpg"

# Step 4: Normalize methylation levels over range
echo "Step 4: Normalizing methylation levels"
stepFeq=1000
stepLen=5
stepRange=2000
infile="${sample}_01_CpG_5.pos"

for i in $(seq 1 ${stepLen}); do
        echo -e "Step frequency is ${i}"
        j=$((${i} - 1))
        stepS=$(echo "$j * ${stepRange} + ${stepFeq}" | bc)
        stepE=$(echo "$i * ${stepRange} + ${stepFeq}" | bc)
        echo -e "Step start is ${stepS} and Step end is ${stepE}"
        outLeft="left_${i}_${stepS}_${stepE}.bed"
        outRight="right_${i}_${stepS}_${stepE}.bed"
        awk -v st=${stepS} -v ed=${stepE} -v OFS='\t' '{start=$2-ed; end=$2-st; if(start<=0){start=1}; print $1, start, end, "left" NR}' ${infile} | \
                bedtools intersect -a stdin -b ${noImBed} -wa -wb -loj | sort -k4,4V -k1,1V -k2,2n -k3,3n | \
                datamash -g 4 count 8 mean 9 | awk '{if($2<5){value=-1}else{value=$3}; print $1 "\t" $2 "\t" value}' > ${outLeft}
        awk -v st=${stepS} -v ed=${stepE} -v OFS='\t' '{start=$3+st; end=$3+ed; print $1, start, end, "right" NR}' ${infile} | \
                bedtools intersect -a stdin -b ${noImBed} -wa -wb -loj | sort -k4,4V -k1,1V -k2,2n -k3,3n | \
                datamash -g 4 count 8 mean 9 | awk '{if($2<5){value=-1}else{value=$3}; print $1 "\t" $2 "\t" value}' > ${outRight}
done

# Combine results
rangeFile="${sample}_range.txt"
paste left* right* | awk -v OFS='\t' '{print $3,$6,$9,$12,$15,$18,$21,$24,$27,$30}' > ${rangeFile}
allFile="${sample}_all.cpg"
outFile="${sample}_merge.cpg"
paste ${allFile} ${rangeFile} > ${outFile}

# Step 5: Normalize relative methylation levels
echo "Step 5: Normalizing methylation levels"
finalOutput="${sample}_normalize.cpg"

awk -v OFS='\t' '{
        n=0
        for (i = 4; i <= 6; i++) { if ($i > 0) values[++n] = $i }
        for (i = 11; i <= 20; i++) { if ($i > 0) values[++n] = $i }
        min = values[1]; max = values[1]
        for (i = 1; i <= n; i++) {
                if (values[i] < min) min = values[i]
                if (values[i] > max) max = values[i]
        }
        if (max == min) {
                print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, "-1", "-1", "-1"
        } else {
                body = (values[1] - min) / (max - min)
                up = (values[2] - min) / (max - min)
                down = (values[3] - min) / (max - min)
                print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, body, up, down
        }
}' ${outFile} > ${finalOutput}

# Step 6: Merge all sample CpG files
cat ./*cpg |sort -k1,1V -k2,2n -k3,3n |datamash -g 1,2,3 -R 2 mean 4,5,6 count 10 |awk '{print $0"\t"$3-$2}' > ../${population}.cpg

echo "DEL methylation level calculation, normalization, and integration all cpgs completed successfully."
