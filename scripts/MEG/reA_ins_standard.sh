#!/bin/bash
genome=~/new/data/gtf/comprehensive_gene_annotation/genome.pos
sample=$1
pop=$2
stepFeq=$3
stepLen=$4
stepRange=2000
ID=$(echo $sample |cut -d"_" -f2)
noImBed=~/fast5/Sample/${ID}/imprint/${ID}_noIm.bed
rootDir=/home/user/liyang/shell/DNAmeth/python/newSV/INS/data/
indir=${rootDir}/${pop}/${sample}/reAlign
infile=${indir}/${sample}_reAlign.result
infileChrauto=${indir}/${sample}_reAlign_chrauto.result
posfile=${indir}/${sample}_source.pos
awk 'NR>1 && $7 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print }' ${infile} > ${infileChrauto}
awk -v OFS='\t' '{print $7,$8,$9,$10,$11,$12}' ${infileChrauto}| grep -E '^chr([1-9]|1[0-9]|2[0-2])\b' > ${posfile}
for i in $(seq 1 ${stepLen}); do
    echo -e "Step frequency is ${i}"
    j=$((${i} - 1))
    stepS=$(echo "$j * ${stepRange} + ${stepFeq}" |bc)
    stepE=$(echo "$i * ${stepRange} + ${stepFeq}" |bc)
    echo -e "Step start is ${stepS} and Step end is ${stepE}"
    outLeft="left_${i}_${stepS}_${stepE}.bed"
    outRight="right_${i}_${stepS}_${stepE}.bed"
    awk -v st=${stepS} -v ed=${stepE} -v OFS='\t' '{start=$2-ed;end=$2-st;if(start<=0){start=1};print $1,start,end,"left"NR}' ${posfile}| 
        bedtools intersect -a stdin -b ${noImBed} -wa -wb -loj |sort -k4,4V -k1,1V -k2,2n -k3,3n |
        datamash -g 4 count 8 mean 9 -R 2 |awk '{if($2<5){value=-1}else{value=$3};print $1"\t"$2"\t"value}' > ${indir}/${outLeft}
    awk -v st=${stepS} -v ed=${stepE} -v OFS='\t' '{start=$3+st;end=$3+ed;print $1,start,end,"right"NR}' ${posfile}| 
        bedtools intersect -a stdin -b ${noImBed} -wa -wb -loj |sort -k4,4V -k1,1V -k2,2n -k3,3n |
        datamash -g 4 count 8 mean 9 -R 2 |awk '{if($2<5){value=-1}else{value=$3};print $1"\t"$2"\t"value}' > ${indir}/${outRight}
done

paste ${indir}/left* ${indir}/right* |awk -v OFS='\t' '{print $3,$6,$9,$12,$15,$18,$21,$24,$27,$30}' > ${indir}/${sample}_range.txt
paste ${infileChrauto} ${indir}/${sample}_range.txt > ${indir}/${sample}_merge.txt

normalizeSh=/home/user/liyang/shell/DNAmeth/python/newSV/INS/re_normalize.sh
norIn=${sample}_merge.txt
norDir=${indir}
norOut=${indir}/${sample}_reAlign_new.result
bash ${normalizeSh} ${norIn} ${norDir} ${norOut}







