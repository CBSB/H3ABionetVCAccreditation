#!/bin/bash

#Checking deduplication options:
# 1. picard
# 2. sambamba
# 3. samblaster

module load samtools
module load picard-tools/2.6.0

picarddir="/usr/src/picard-tools/picard-tools-2.6.0"

 java -Xmx10g -XX:-UseGCOverheadLimit -jar $picarddir/picard.jar MarkDuplicates \
	 I=$align_res/$samplename.sorted.bam\
	 O=$align_res/$samplename.dedup.bam \
	 M=$reports/$samplename.dedup.txt

exit_code=$?
if [ ! $exit_code -eq 0 ]; then
	echo 'Marking duplicates did NOT work'
	exit
fi

if [ ! -s ${align_res}/$samplename.dedup.bam ]; then
	echo 'Marking duplicates file empty!'
	exit
fi

numalign=$(samtools view -c ${align_res}/$samplename.dedup.bam)

if [  $numalign -eq 0 ]; then
	echo 'Empty bam file'
	exit
fi

