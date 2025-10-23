#!/bin/bash

# BAM -> FASTQ + NanoStat N50 & estimated error rate from mean Q
# Usage: fastq_qc_nanostat.sh <BAM> <SAMPLE_ID> <OUTDIR> [THREADS]
# Requires: samtools, NanoStat, gzip

BAM="$1"
SAMPLE="$2"
OUTDIR="$3"
THREADS="${4:-8}"

mkdir -p "${OUTDIR}/nanostat"

# 1) BAM -> FASTQ (gz)
samtools fastq -@ "${THREADS}" "${BAM}" | gzip -c > "${OUTDIR}/${SAMPLE}.fastq.gz"

# 2) NanoStat
NanoStat --fastq "${OUTDIR}/${SAMPLE}.fastq.gz" --outdir "${OUTDIR}/nanostat" --name "NanoStats"

REPORT="${OUTDIR}/nanostat/NanoStats.txt"

# 3) Parse N50 & Mean Q, then estimate error rate = 10^(-Q/10)
N50=$(grep -i -m1 "^N50 read length" "${REPORT}" | awk '{print $NF}')
MEAN_Q=$(grep -i -m1 "^Mean read quality" "${REPORT}" | awk '{print $NF}')
ERR_RATE=$(awk -v q="${MEAN_Q}" 'BEGIN{printf "%.6f", 10^(-q/10)}')

# 4) Summary
{
    echo -e "Sample\tN50\tMeanQ\tErrorRate"
    echo -e "${SAMPLE}\t${N50}\t${MEAN_Q}\t${ERR_RATE}"
        
} > "${OUTDIR}/summary_nanostat.txt"

echo "Done."
echo "FASTQ: ${OUTDIR}/${SAMPLE}.fastq.gz"
echo "NanoStat report: ${REPORT}"
echo "Summary: ${OUTDIR}/summary_nanostat.txt"
