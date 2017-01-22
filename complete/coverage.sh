#!/bin/bash

module load bedtools/2.26.0
set -x

inputbam=$1
targetedbed=$2
reports=$3
samplename=$4

sed -n '/chr1:/p' $targetedbed  > $reports/trueseq.targeted.regions.chr1.hg19.bed

bedtools coverage -hist -b $inputbam -a $reports/targeted.regions.chr1.hg19.bed > $reports/$samplename.coverage.hist.txt

# -dz  with zero-based coordinates
# -d with one-based coordinates
