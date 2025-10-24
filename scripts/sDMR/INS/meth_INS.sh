#!/bin/bash
inbam=$1
inbed=$2
outdir=$3

python extractReadFromINS.py --bam ${inbam} --vcf ${inbed} --outdir ${outdir} 
