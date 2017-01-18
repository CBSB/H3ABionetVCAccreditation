#!/bin/bash

module load bcftools/1.3.1 
gatk_file=$1
freebayes_file=$2

bcftools isec   $gatkvcf_file $freebayesvcf_file  -p -c both

bcftools stats $gatk_file  $freebayes_file -c both > targeted.vcfs.stats
