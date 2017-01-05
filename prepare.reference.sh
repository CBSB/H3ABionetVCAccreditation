#!/bin/bash

set -x

projectdir=/home/aeahmed/assessment
reference=$projectdir/genome/hg19/
cd $reference

:<<'comment_done_already'
sudo gunzip $reference/ucsc*
cd $reference

set +x
module load bwa

set -x
bwa index -a bwtsw $reference/ucsc.hg19.fasta
 
comment_done_already

set +x 
module load novocraft

set -x
novoindex ucsc.hg19.nix ucsc.hg19.fasta

cd ../ref/

module load bwa

bwa index -a bwtsw human_g1k_v37_chr20.fa
ll -htr

samtools
module load samtools/1.3.1
samtools
samtools faidx human_g1k_v37_chr20.fa

picard
module avail
module load picard-tools/2.6.0
module list

java -jar /usr/src/picard-tools/picard-tools-2.6.0/picard.jar CreateSequenceDictionary     OUTPUT=human_g1k_v37_chr20.dict     REFERENCE=human_g1k_v37_chr20.fa

