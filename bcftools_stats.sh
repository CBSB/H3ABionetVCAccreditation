#!/bin/bash

mdule load bcftools/1.3.1 

bcftools stats CBSB_Khartoum.raw.calls.freebayes.targeted.vcf.gz CBSB_Khartoum.raw.calls.haplotypecaller.targeted.vcf.gz -c both > targeted.vcfs.stats 

plot-vcfstats -p Intersection_of_targettedVCF/ targeted.vcfs.stats 
