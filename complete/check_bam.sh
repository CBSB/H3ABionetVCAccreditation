#!/bin/bash

## this script checks if a bam file is valid! 

set -x
file=$1
stage=$2
reports=$3
samplename=$4

if [ ! -s $file ]; then
	echo $stage file: $file is empty!
	exit
fi
numalign=$(samtools view -c $file) 
if [  $numalign -eq 0 ]; then
	echo $stage file: $file is empty!
	exit
fi

if [[ $stage == Alignment* ]] ; then
	samtools flagstat $file >> $reports/${samplename}.summary.txt
fi
