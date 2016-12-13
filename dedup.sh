#!/bin/bash

module load picard-tools/2.6.0

picarddir="/usr/src/picard-tools/picard-tools-2.6.0" #similar in all
align="/home/classrooms/Workshop_low_pass/align"     #depending on your data

module load samtools

samtools sort -o $align/sorted.bam -@ 2 $align/HG00120.lowcoverage.chr20.smallregion_2i_addingRG.bam


java -jar $picarddir/picard.jar MarkDuplicates \
      I=$align/sorted.bam\
      O=$align/marked_duplicates.bam \
      M=$align/marked_dup_metrics.txt

 
