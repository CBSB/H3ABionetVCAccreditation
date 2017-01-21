#!bin/bash

Module load freebayes/1.1.0 

freebayes -= -f /home/mirrors/gatk_bundle/2.8/hg19/hg19/ucsc.hg19.fasta CBSB_Khartoum.dedup.bam > var.vcf
