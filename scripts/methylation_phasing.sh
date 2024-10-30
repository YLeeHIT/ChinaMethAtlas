#!/bin/n/bash

# Methylation Phasing Pipeline Script
# This script performs phasing on methylation calling data through multiple steps:
# - Step 1: Phasing BAM file and generating phased VCF using Clair3.
# - Step 2: Processing methylation calling results and phasing them using NanoMethPhase.
# - Step 3: Generating haplotype-specific methylation frequencies.
# 
# Default Parameters:
# clair3: Path to Clair3 executable.
# nanomethphase: Path to NanoMethPhase executable.
# reference: GRch38 Human Reference genome file.
# ont_model: Model file for Oxford Nanopore sequencing.
# threads: Number of CPU threads to use.
# deep_threshold: Threshold for deepsignal methylation calls.

# Define file paths and parameters
meth_file="methylation_calls.tsv"  # Methylation calling file in TSV format
inbam="input.bam"  # Input BAM file
vcf_file="phased_variants.vcf.gz"  # Phased VCF file
sampleID="sample_name"  # Sample ID for output naming
clair3="run_clair3.sh"  # Specify Clair3 executable
nanomethphase="nanomethphase.py"  # Specify NanoMethPhase script
reference="hg38.fa"  # Reference genome file
ont_model="ont_model"  # ONT model file
platform="ont"
threads=16
deep_threshold=0.6  # Threshold for deepsignal methylation calls

# Step 1: Phasing BAM and generating phased VCF using Clair3
echo "Starting Step 1: Phasing BAM and generating phased VCF..."
source activate clair
bash ${clair3} -b ${inbam} -f ${reference} -m ${ont_model} -p ${platform} -t ${threads} \
        -o ${sampleID}_vcf --sample_name=${sampleID} --enable_phasing
vcf="${sampleID}_vcf/phased.vcf.gz"
echo "Step 1: Clair3 phasing completed."

# Step 2: Phasing methylation calling results using NanoMethPhase
echo "Starting Step 2: Phasing methylation calls..."
source activate nanomethphase

# Step 2.1 Construct DeepSignal-compatible file format
echo "Starting Step 2.1: Constructing DeepSignal-compatible TSV..."
deep_format="${sampleID}_deepSignal.tsv"
grep -v "^read_id" ${meth_file} | awk -v OFS='\t' 'BEGIN{print "chrom","pos","strand","pos_in_strand","readname","prob_0","prob_1","called_label","k_mer"} \
        {if($11>0.5){zo=1}else{zo=0}; print $4,$3,$6,$2,$1,"t",1-$11,$11,zo,$15}' > ${deep_format}
echo "Step 2.1: DeepSignal-compatible file created."

# Step 2.2 Pretreatment of methylation calling for phasing
echo "Starting Step 2.2: Preprocessing methylation calls for phasing..."
bed_file="${sampleID}.bed.gz"
python ${nanomethphase} methyl_call_processor -mc ${deep_format} -tc deepsignal:${deep_threshold} -t ${threads} \
        | sort -k1,1 -k2,2n -k3,3n \
        | bgzip > ${bed_file} && tabix -p bed ${bed_file}
echo "Step 2.2: Methylation call preprocessing completed."

# Step 2.3 Phasing methylome data
echo "Starting Step 2.3: Phasing methylome..."
python ${nanomethphase} phase -b ${inbam} -v ${vcf} -mc ${bed_file} -o ${sampleID} \
        --reference ${reference} -of bam,methylcall -t ${threads} -mt cpg -ind -ow
echo "Step 2.3: Methylome phasing completed."

# Step 3: Generating haplotype-specific methylation frequencies
echo "Starting Step 3: Calculating haplotype-specific methylation frequencies..."

# Hap1 methylation frequency
sed '1d' ${sampleID}_NanoMethPhase_HP1_MethylFrequency.tsv | \
    awk -F'\t' '{if ($4=="-") {$2=$2-1;$3=$3-1}; print $1,$2,$3,$5,$6}' OFS='\t' | \
    sort -k1,1 -k2,2n | \
    datamash -g1,2,3 sum 4,5 | \
    awk -F'\t' '{print $0,$5/$4}' OFS='\t' | \
    sed '1i chromosome\tstart\tend\tNumOfAllCalls\tNumOfModCalls\tMethylFreq' > ${sampleID}_HP1_MethylFrequency.tsv

# Hap2 methylation frequency
sed '1d' ${sampleID}_NanoMethPhase_HP2_MethylFrequency.tsv | \
    awk -F'\t' '{if ($4=="-") {$2=$2-1;$3=$3-1}; print $1,$2,$3,$5,$6}' OFS='\t' | \
    sort -k1,1 -k2,2n | \
    datamash -g1,2,3 sum 4,5 | \
    awk -F'\t' '{print $0,$5/$4}' OFS='\t' | \
    sed '1i chromosome\tstart\tend\tNumOfAllCalls\tNumOfModCalls\tMethylFreq' > ${sampleID}_HP2_MethylFrequency.tsv

echo "Step 3: Haplotype-specific methylation frequencies calculated."
echo "Methylation phasing pipeline completed successfully."
