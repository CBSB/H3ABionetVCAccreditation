#!/bin/bash

gatkdir="/usr/src/gatk/gatk-3.6/"
referencedir="/home/classrooms/Workshop_low_pass/ref/"
bams="/home/classrooms/Workshop_low_pass/align/"

module load gatk/3.6

java -jar $gatkdir/GenomeAnalysisTK.jar \
                -T HaplotypeCaller\
                -R $referencedir/human_g1k_v37_chr20.fa\
		-I $bams/recal.default.bam\
		--emitRefConfidence GVCF \
		--dbsnp $referencedir/dbsnp_135.b37.chr20.smallregion.vcf \
		-o /home/classrooms/Workshop_low_pass/variant/output.raw.snps.indels.g.vcf


