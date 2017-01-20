#!/bin/bash

module load bedtools/2.26.0

inputbam=$1
targetedbed=$2
reports=$3
samplename=$4

bedtools coverage -hist -b $inputbam -a $targetedbed > $reports/$samplename.coverage.hist.txt

# -dz  with zero-based coordinates
# -d with one-based coordinates
