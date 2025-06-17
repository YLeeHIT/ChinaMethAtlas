#!/bin/bash
# Function: annotate mobile INS
# Input:
#   1) pop: [xizang, north, south] population
#   2）outDir: [best named anno, absolute path] output directory
# Output
#   1) all result located in anno directory
# P.S.
#   1) If you want change the output directory, you should modify the rootDir 
#   2) You must run the population integration shell before doing this process
# Example:
#   bash annoINSSource.sh xizang xizang_reAlign.result anno
#

pop=$1
#inIns=$2
outPath=$2
rootDir=/home/user/liyang/shell/DNAmeth/python/newSV/INS/data/repop/${pop}/out
inIns=${rootDir}/${pop}_reAlign.result
outDir=${rootDir}/${outPath}

# only for find pop_num.cpg files and do not perform frequency processing
if [ "${pop}" == "xizang" ];then
    num=10
else
    num=48
fi

# basic annotation files
gene="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/all-gene/gene_2kb"
exon="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/genic/exon"
promoter="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/genic/promoter"
CGI="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/CGI/CGI"
DEL="/home/user/liyang/shell/DNAmeth/python/newSV/DEL/Pop/${pop}/${pop}_${num}.cpg"
DEL_all="/home/user/liyang/shell/DNAmeth/python/newSV/data/pop_merge/${pop}/fromDEL2/${pop}_del2.bed"
SINE="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/repeat/SINE"
LINE="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/repeat/LINE"
LTR="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/repeat/LTR"
Retroposon="/home/user/liyang/new/data/gtf/comprehensive_gene_annotation/repeat/Retroposon"

test -d ${outDir} || mkdir ${outDir}
cd ${outDir}

inInsBed="${pop}_ins.bed"
inInsSourceBed="${pop}_ins_source.bed"
annoGene="${pop}_ins.gene"
annoCGI="${pop}_ins.CGI"
annoExon="${pop}_ins.exon"
annoPromoter="${pop}_ins.promoter"
annoDEL="${pop}_ins.Del"
annoDELAll="${pop}_ins.DelAll"
awk 'NR>1{print $1"\t"$2-100"\t"$2+100"\t"$14}' ${inIns} > ${inInsBed}
awk 'NR>1{print $7"\t"$8"\t"$8+$9"\t"$14}' ${inIns} > ${inInsSourceBed}

bedtools intersect -a ${inInsBed} -b ${gene} -wa -wb > ${annoGene}
bedtools intersect -a ${inInsBed} -b ${exon} -wa -wb > ${annoExon}
bedtools intersect -a ${inInsBed} -b ${promoter} -wa -wb > ${annoPromoter}
bedtools intersect -a ${inInsBed} -b ${CGI} -wa -wb > ${annoCGI}
bedtools intersect -a ${inInsSourceBed} -b ${DEL} -wa -wb -f 0.1 > ${annoDEL}
bedtools intersect -a ${inInsSourceBed} -b ${DEL_all} -wa -wb -f 0.1 > ${annoDELAll}
geneNum=$(wc -l ${annoGene} |cut -f1 -d" ")
exonNum=$(wc -l ${annoExon} |cut -f1 -d" ")
promNum=$(wc -l ${annoPromoter} |cut -f1 -d" ")
CGINum=$(wc -l ${annoCGI} |cut -f1 -d" ")
DELNum=$(wc -l ${annoDEL} |cut -f1 -d" ")
DELAllNum=$(wc -l ${annoDELAll} |cut -f1 -d" ")
#bcftools view -H xizang_hwe-6_heter-2.del.vcf |awk 'BEGIN{OFS="\t"}{split($8,x,";");split(x[3],y,"=");split(x[4],z,"=");print $1,$2,y[2],z[2]}' > fromDEL2/xizang_del2.bed

# MGEs
annoLINE="${pop}_ins.LINE"
annoSINE="${pop}_ins.SINE"
annoLTR="${pop}_ins.LTR"
annoRetroposon="${pop}_ins.Retroposon"
bedtools intersect -a ${inInsBed} -b ${LINE} -wa -wb > ${annoLINE}
bedtools intersect -a ${inInsBed} -b ${SINE} -wa -wb > ${annoSINE}
bedtools intersect -a ${inInsBed} -b ${LTR} -wa -wb > ${annoLTR}
bedtools intersect -a ${inInsBed} -b ${Retroposon} -wa -wb > ${annoRetroposon}
LINENum=$(wc -l ${annoLINE} |cut -f1 -d" ")
SINENum=$(wc -l ${annoSINE} |cut -f1 -d" ")
LTRNum=$(wc -l ${annoLTR} |cut -f1 -d" ")
RetroposonNum=$(wc -l ${annoRetroposon} |cut -f1 -d" ")

# statistic
outSta="${pop}_anno.sta"
echo -e "Type\tNum" > ${outSta}
echo -e "Gene\t${geneNum}" >> ${outSta}
echo -e "Exon\t${exonNum}" >> ${outSta}
echo -e "Promoter\t${promNum}" >> ${outSta}
echo -e "CGI\t${CGINum}" >> ${outSta}
echo -e "DEL\t${DELNum}" >> ${outSta}
echo -e "DEL_all\t${DELAllNum}" >> ${outSta}
echo -e "LINE\t${LINENum}" >> ${outSta}
echo -e "SINE\t${SINENum}" >> ${outSta}
echo -e "LTR\t${LTRNum}" >> ${outSta}
echo -e "Retroposon\t${RetroposonNum}" >> ${outSta}
echo -e "Annotation of ${pop} has finished"


