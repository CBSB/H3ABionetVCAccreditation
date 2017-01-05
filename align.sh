#!/bin/bash

module load bwa/0.7.15

bwa mem -M -R "@RG\tID:synthetic\tPL:illumina\tPU:synthetic\tLB:synthetic\tDT:2016-7-1\tSM:NA12878-Garvan"  /home/classrooms/Workshop_low_pass/ref/human_g1k_v37_chr20.fa  /home/classrooms/Workshop_low_pass/fastq/HG00120.lowcoverage.chr20.smallregion_1.fastq /home/classrooms/Workshop_low_pass/fastq/HG00120.lowcoverage.chr20.smallregion_2.fastq > /home/classrooms/Workshop_low_pass/align/HG00120.lowcoverage.chr20.smallregion_2i_addingRG.sam


module load samtools
samtools view -bS /home/classrooms/Workshop_low_pass/align/HG00120.lowcoverage.chr20.smallregion_2i_addingRG.sam -o /home/classrooms/Workshop_low_pass/align/HG00120.lowcoverage.chr20.smallregion_2i_addingRG.bam

samtools flagstat /home/classrooms/Workshop_low_pass/align/HG00120.lowcoverage.chr20.smallregion_2i_addingRG.bam >> summary.txt
