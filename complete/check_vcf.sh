#!/bin/bash

## You need to add gatk VariantEval with known sites: $dbsnp129 

set -x
file=$1
stage=$2
reports=$3
gatkdir=$4
reference=$5
dbsnp129=$6
email=$7

module load gatk/3.6

java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T VariantEval\
        -R $reference\
        -o $reports/variant.eval.grp\
        --eval:$stage  $file\
        -D $dbsnp129\
        -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary


