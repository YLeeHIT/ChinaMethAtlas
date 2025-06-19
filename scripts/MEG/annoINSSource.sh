#!/bin/bash

# Description:
# Annotates INS regions with various genomic features including gene, exon, promoter,
# CpG islands (CGI), deletion overlap, and repeat elements (LINE, SINE, LTR, Retroposon).
# Outputs per-feature annotation and summary statistics.

# ----------------------
# Input parameters
# ----------------------
pop=$1               # Population name
outPath=$2           # Subdirectory name for output
num=$3               # Identifier for DEL group (e.g., 1, 2, 3)

# ----------------------
# Configurable paths
# ----------------------
ROOT_DIR="${ROOT_DIR:-../../Demo/for_MEG}"
IN_INS="${ROOT_DIR}/${pop}_reAlign.result"
OUT_DIR="${ROOT_DIR}/${outPath}"

# Genomic annotation references
GENE="${GENE:-../../data/gene_2kb}"
EXON="${EXON:-../../data/exon}"
PROMOTER="${PROMOTER:-../../data/promoter}"
CGI="${CGI:-../../data/CGI}"
DEL="${DEL:-../../data/${pop}_${num}.cpg}"
DEL_ALL="${DEL_ALL:-../../data/${pop}_del2.bed}"
SINE="${SINE:-../../data/SINE}"
LINE="${LINE:-../../data/LINE}"
LTR="${LTR:-../../data/LTR}"
RETROPOSON="${RETROPOSON:-../../data/Retroposon}"

# ----------------------
# Create output directory
# ----------------------
test -d "${OUT_DIR}" || mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# ----------------------
# Prepare INS BED files
# ----------------------
INS_BED="${pop}_ins.bed"
INS_SOURCE_BED="${pop}_ins_source.bed"

awk 'NR > 1 {print $1, $2 - 100, $2 + 100, $14}' OFS='\t' "${IN_INS}" > "${INS_BED}"
awk 'NR > 1 {print $7, $8, $8 + $9, $14}' OFS='\t' "${IN_INS}" > "${INS_SOURCE_BED}"

# ----------------------
# Perform annotation using bedtools intersect
# ----------------------
bedtools intersect -a "${INS_BED}" -b "${GENE}"      -wa -wb > "${pop}_ins.gene"
bedtools intersect -a "${INS_BED}" -b "${EXON}"      -wa -wb > "${pop}_ins.exon"
bedtools intersect -a "${INS_BED}" -b "${PROMOTER}"  -wa -wb > "${pop}_ins.promoter"
bedtools intersect -a "${INS_BED}" -b "${CGI}"       -wa -wb > "${pop}_ins.CGI"
bedtools intersect -a "${INS_SOURCE_BED}" -b "${DEL}"     -wa -wb -f 0.1 > "${pop}_ins.Del"
bedtools intersect -a "${INS_SOURCE_BED}" -b "${DEL_ALL}" -wa -wb -f 0.1 > "${pop}_ins.DelAll"

bedtools intersect -a "${INS_BED}" -b "${LINE}"       -wa -wb > "${pop}_ins.LINE"
bedtools intersect -a "${INS_BED}" -b "${SINE}"       -wa -wb > "${pop}_ins.SINE"
bedtools intersect -a "${INS_BED}" -b "${LTR}"        -wa -wb > "${pop}_ins.LTR"
bedtools intersect -a "${INS_BED}" -b "${RETROPOSON}" -wa -wb > "${pop}_ins.Retroposon"

# ----------------------
# Count annotated entries
# ----------------------
count_lines() {
    wc -l "$1" | cut -f1 -d" "
}

geneNum=$(count_lines "${pop}_ins.gene")
exonNum=$(count_lines "${pop}_ins.exon")
promoterNum=$(count_lines "${pop}_ins.promoter")
cgiNum=$(count_lines "${pop}_ins.CGI")
delNum=$(count_lines "${pop}_ins.Del")
delAllNum=$(count_lines "${pop}_ins.DelAll")
lineNum=$(count_lines "${pop}_ins.LINE")
sineNum=$(count_lines "${pop}_ins.SINE")
ltrNum=$(count_lines "${pop}_ins.LTR")
retroNum=$(count_lines "${pop}_ins.Retroposon")

# ----------------------
# Output summary statistics
# ----------------------
STAT_OUT="${pop}_anno.sta"
{
    echo -e "Type\tNum"
    echo -e "Gene\t${geneNum}"
    echo -e "Exon\t${exonNum}"
    echo -e "Promoter\t${promoterNum}"
    echo -e "CGI\t${cgiNum}"
    echo -e "DEL\t${delNum}"
    echo -e "DEL_all\t${delAllNum}"
    echo -e "LINE\t${lineNum}"
    echo -e "SINE\t${sineNum}"
    echo -e "LTR\t${ltrNum}"
    echo -e "Retroposon\t${retroNum}"
} > "${STAT_OUT}"

echo "INS annotation for population ${pop} completed. Summary saved to ${STAT_OUT}"



