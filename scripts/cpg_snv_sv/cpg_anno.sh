#!/bin/bash

# ==========================================
# CpG region annotation counter (single base dir for annotations)
# Usage:
#   annotate_cpg.sh <SAMPLE_ID> <IN_DIR> <OUT_DIR> <ANN_BASE_DIR>
# ==========================================

SAMPLE_ID="$1"
IN_DIR="$2"
OUT_DIR="$3"
ANN_BASE_DIR="$4"

# Derive sub-directories from the base dir
GENIC_DIR="${ANN_BASE_DIR}/genic"
REPEAT_DIR="${ANN_BASE_DIR}/repeat"
CGI_DIR="${ANN_BASE_DIR}/CGI"
ENHANCER_DIR="${ANN_BASE_DIR}/enhancer"
TELOMERE_DIR="${ANN_BASE_DIR}/Telomere"

sampleName="${IN_DIR}/${SAMPLE_ID}_chrauto.bed"
outSample="${OUT_DIR}/${SAMPLE_ID}.anno"
outFile="${OUT_DIR}/${SAMPLE_ID}_transpose.anno"

mkdir -p "${OUT_DIR}"
: > "${outSample}"

gene_list=("gene" "Intergenic" "intron" "exon" "CDS" "UTR" "TSS200" "TSS1500" "promoter")
repeat_list=("SINE" "LINE" "Simple_repeat" "LTR" "Low_complexity" "Satellite" "Retroposon" "RC")
CGI_list=("CGI" "shore" "shelf" "interCGI")
enhancer_list=("DNase" "H3k27ac" "H3k27me3" "H3k36me3" "H3k4me1" "H3k4me2" "H3k4me3" "H3k9ac" "H3k9me3" "H4k20me1")
Telomere_list=("telomere_1wb" "telomere_10wb" "telomere_100wb")

meth_type=("gene" "repeat" "CGI" "enhancer" "Telomere")

for i in "${meth_type[@]}"; do
    case "${i}" in
        gene)
            echo "gene"
            for gene_type in "${gene_list[@]}"; do
                echo "Doing ${gene_type}"
                num=$(bedtools intersect -a "${sampleName}" -b "${GENIC_DIR}/${gene_type}" -wa | cut -f1-3 | sort -u | wc -l)
                echo -e "${gene_type}\t${num}" >> "${outSample}"
            done
            ;;

        repeat)
            echo "repeat"
            for repeat_type in "${repeat_list[@]}"; do
                echo "Doing ${repeat_type}"
                num=$(bedtools intersect -a "${sampleName}" -b "${REPEAT_DIR}/${repeat_type}" -wa | cut -f1-3 | sort -u | wc -l)
                echo -e "${repeat_type}\t${num}" >> "${outSample}"
            done
            ;;

        CGI)
            echo "CGI"
            for CGI_type in "${CGI_list[@]}"; do
                echo "Doing ${CGI_type}"
                if [ "${CGI_type}" = "interCGI"  ]; then
                    # interCGI: sites NOT intersecting with the provided interCGI set
                    num=$(bedtools intersect -a "${sampleName}" -b "${CGI_DIR}/${CGI_type}" -wa -v | cut -f1-3 | sort -u | wc -l)
                else
                    num=$(bedtools intersect -a "${sampleName}" -b "${CGI_DIR}/${CGI_type}" -wa | cut -f1-3 | sort -u | wc -l)
                fi
                echo -e "${CGI_type}\t${num}" >> "${outSample}"
            done
            ;;

        enhancer)
            echo "enhancer"
            for enhancer_type in "${enhancer_list[@]}"; do
                echo "Doing ${enhancer_type}"
                num=$(bedtools intersect -a "${sampleName}" -b "${ENHANCER_DIR}/${enhancer_type}" -wa | cut -f1-3 | sort -u | wc -l)
                echo -e "${enhancer_type}\t${num}" >> "${outSample}"
            done
            ;;

        Telomere)
            echo "Telomere"
            for Telomere_type in "${Telomere_list[@]}"; do
                echo "Doing ${Telomere_type}"
                num=$(bedtools intersect -a "${sampleName}" -b "${TELOMERE_DIR}/${Telomere_type}" -wa | cut -f1-3 | sort -u | wc -l)
                echo -e "${Telomere_type}\t${num}" >> "${outSample}"
            done
            ;;

        *) echo "Unknown type: ${i}" ;;
    esac
done


datamash transpose < "${outSample}" > "${outFile}"

echo "Done."
echo "Annotation table  : ${outSample}"
echo "Transposed output : ${outFile}"
