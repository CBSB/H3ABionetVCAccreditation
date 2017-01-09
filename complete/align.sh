#!/bin/bash

## this script should contain options for doing alignment. Either bwa or novoalign

# the output should be saved in align_res, and specify the tool used for alginment and resources!
set -x
read1=$1
read2=$2
bwa_index=$3
novoalign_index=$4
rgheader=$5
align_res=$6
samplename=$7
reports=$8
tool=$9

module load bwa/0.7.15
module load novocraft/3.06.04
module load samtools/1.3.1

if [ `expr ${#tool}` -lt 1  ]; then
	## do both bwa and novoalign
	for i in $(seq 1 1 8); do
		bwa mem -M -t $i -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ $i -bS > ${align_res}/$samplename.aligned.bwa.bam
		./check_bam.sh ${align_res}/$samplename.aligned.bwa.bam Alignment_bwa $reports $samplename
		novoalign -c $i -d $novoalign_index -f $read1 $read2 -o SAM $rgheader | samtools view -@ $i -bS  > ${align_res}/$samplename.aligned.novoalign.bam
		./check_bam.sh ${align_res}/$samplename.aligned.novoalign.bam Alignment_novoalign $reports $samplename
	done
else
	# do novoalign or bwa
	if [ $tool == 'bwa' ]; then
		bwa mem -M -t 4 -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ 4 -bS > ${align_res}/$samplename.aligned.bwa.bam
		./check_bam.sh ${align_res}/$samplename.aligned.bwa.bam Alignment_bwa $reports $samplename
	else
		novoalign -c 4 -d $novoalign_index -f $read1 $read2 -o SAM $rgheader | samtools view -@ 4 -bS  > ${align_res}/$samplename.aligned.novoalign.bam
		./check_bam.sh ${align_res}/$samplename.aligned.novoalign.bam Alignment_novoalign $reports $samplename
	fi
fi

