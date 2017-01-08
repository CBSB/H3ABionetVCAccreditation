#!/bin/bash

## this script should contain options for doing alignment. Either bwa or novoalign

# the output should be saved in align_res, and specify the tool used for alginment and resources!


read1=$1
read2=$2
bwa_index=$3
novoalign_index=$4
rgheader=$5
align_res=$6

tool=$7

module load bwa/0.7.15
module load novocraft/3.06.04
module load samtools/1.3.1

if [ `expr ${#tool}` -lt 1  ]; then
	## do both bwa and novoalign
	bwa mem -M -t 4 -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ 4 -bS > ${align_res}/$samplename.bam
	novoalign -d $novoalign_index -f $read1 $read2 -o SAM $rgheader -c 4
else
# do novoalign or bwa
fi

exit_code=$?
if [ ! $exit_code -eq 0 ]; then
       	echo 'Alignment did NOT work'
        exit
fi
if [ ! -s ${align_res}/$samplename.bam ]; then
	echo 'Alignment file empty!'
	exit
fi
numalign=$(samtools view -c ${align_res}/$samplename.bam)
if [  $numalign -eq 0 ]; then
	echo 'Empty bam file'
	exit
fi


