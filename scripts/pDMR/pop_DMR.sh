#!/bin/bash

# ==========================================
# Group-wise DMR and DMC Detection (Metilene)
# Steps:
#   1. Prepare group input files
#   2. Identify DMRs between two groups
#   3. Filter and summarize results
# ==========================================

# -------- Input Parameters --------
g1=$1             # Group 1 name
g2=$2             # Group 2 name
rootdir=$3        # Project root directory
threads=16

indir="${rootdir}/input"
outdir="${rootdir}/output"
labdir="${rootdir}/group"

# -------- Step 1: Construct Input --------
echo "### Step 1: Construct input file ###"

g1_ID=$(cat "${labdir}/${g1}")
g2_ID=$(cat "${labdir}/${g2}")

mkdir -p "${indir}"
cd "${rootdir}/raw"

metilene_input.pl \
    --in1 ${g1_ID} \
    --in2 ${g2_ID} \
    --h1 ${g1} \
    --h2 ${g2} \
    --out "${indir}/${g1}_vs_${g2}.bed"

echo "Input construction finished."


# -------- Step 2: Identify DMRs and DMCs --------
echo "### Step 2: Identify DMRs and DMCs ###"

mkdir -p "${outdir}"

infile="${indir}/${g1}_vs_${g2}.bed"
outfile="${outdir}/${g1}_vs_${g2}_DMRs.txt"
ffile="${outdir}/${g1}_vs_${g2}_filter_DMRs.bed"

maxdist=1000
mincpgs=5
minMethDiff=0.1

# Detect DMRs
metilene -M ${maxdist} -m ${mincpgs} -d ${minMethDiff} -t ${threads} \
    -f 1 -a ${g1} -b ${g2} ${infile} | sort -V -k1,1 -k2,2n > ${outfile}

echo "DMR and DMC identification finished."

# -------- Step 3: Filter --------
echo "### Step 3: Filter ###"

# Filter DMRs: q < 0.05 and |delta| > 0.1
awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tq-value\tdelta\tnum\tpMWU\tp2D\tmeang1\tmeang2"}
     {if($4<0.05 && ($5>=0.1 || $5<=-0.1)) print $1,$2,$3,$5}' ${outfile} > ${ffile}

