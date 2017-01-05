#!/bin/bash

module load bedtools/2.26.0

result=/home/aeahmed/assessment/results/depth
input=~//assessment/results/CBSB_Khartoum/align/CBSB_Khartoum.dedup.bam
targetedbed=/home/aeahmed/assessment/genome//ZachUses_Illumina_truseq_exome_targeted_regions.hg19.chr1.bed
targetedbed2=/home/aeahmed/assessment/genome/targets.numeric.chroms.bed2
#bedtools genomecov -ibam $input > $result/sample.coverage.hist.txt 

#bedtools genomecov -ibam $input -bga > $result/sample.coverage.bedg

bedtools coverage -hist -b $input  -a $targetedbed > $result/sample.exome.coverage.hist.txt

echo "Pipeline (vanialla gatk) run finished! yeeeey!"| mail -s "accreditation pipeline" "azzaea@gmail.com"

# -dz  with zero-based coordinates
# -d with one-based coordinates
