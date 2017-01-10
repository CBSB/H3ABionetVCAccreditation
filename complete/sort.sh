#!/bin/bash

#This should contains options for sorting. 
# 1. samtools
# 2. picard
# 3. sambamba

module load samtools
picarddir="/usr/src/picard-tools/picard-tools-2.6.0"
samtools sort -o $align_res/$samplename.sorted.bam -@ 4 ${align_res}/$samplename.bam

