#!/bin/bash

module load bwa/0.7.15

bwa mem -M -R "@RG\tID:Set1\tSM:CBSB_Khartoum\tPL:illumina\tPU:synthetic\tLB:synthetic\tDT:2016-12-12" ~/gatkbundle/hg19/ucsc.hg19.fasta ~/exercise/reads/H3A_NextGen_assessment.Chr1_50X.set1_read1.fq  ~/exercise/reads/H3A_NextGen_assessment.Chr1_50X.set1_read2.fq  > ~/exercise/set1_addingRG.sam



module load samtools
samtools view -bS ~/exercise/alignment/set1_addingRG.sam
 -o ~/exercise/alignment/set1_addingRG.bam

samtools flagstat ~/exercise/alignment/set1_addingRG.bam>> alignment.summary.txt
