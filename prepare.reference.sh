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


