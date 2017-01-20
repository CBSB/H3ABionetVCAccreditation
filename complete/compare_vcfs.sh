#!/bin/bash

module load bcftools/1.3.1 
gatk_file=$1
freebayes_file=$2
var_res=$3
dbsnp138=$4

bcftools stats $gatk_file  $freebayes_file -c both > $var_res/targeted.vcfs.stats
plot-vcfstats -p $var_res/stats.plots/ $var_res/targeted.vcfs.stats

bcftools stats $gatk_file  $dbsnp138 -c both > $var_res/hc.dbsnp.vcfs.stats
plot-vcfstats -p $var_res/stats.plots.hc.dbsnp/ $var_res/hc.dbsnp.vcfs.stats

bcftools stats $freebayes_file  $dbsnp138 -c both > $var_res/fb.dbsnp.vcfs.stats
plot-vcfstats -p $var_res/stats.plots.fb.dbsnp/ $var_res/fb.dbsnp.vcfs.stats


