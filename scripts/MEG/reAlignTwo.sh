#!/bin/bash

Chr=$1
Start=$2
End=$3
noImBed=$4
genome="../../genome.pos"

# calculate the methlyation of (chr, start, end)
echo -e "${Chr}\t${Start}\t${End}" |bedtools intersect -a ${noImBed} -b - |datamash mean 5 -R 2
