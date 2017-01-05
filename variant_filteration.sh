#!/bin/bash

projectdir=/home/aeahmed/assessment
referencedir=${projectdir}/genome/hg19/
reference=${referencedir}/ucsc.hg19.fasta

result=${projectdir}/results

module load gatk/3.6
module load picard-tools/2.6.0

set -x
############################################ Specifying path for java tools:i
picarddir="/usr/src/picard-tools/picard-tools-2.6.0" #similar in all
gatkdir="/usr/src/gatk/gatk-3.6/"

samplename=CBSB_Khartoum
align_res=${result}/${samplename}/align
inputbam=$align_res/$samplename.dedup.bam
rawvcf=${result}/joint_called.vcf 
######################### The measurements

java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T SelectVariants \
        -R $reference \
	-V $rawvcf \
	-L chr1 \
	-selectType SNP \
	-o $result/raw_snps.vcf

java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R $reference \
        -V $rawvcf \
        -L chr1 \
        -selectType INDEL \
        -o $result/raw_indels.vcf


java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $reference \
    	-V $result/raw_snps.vcf \
	--filterExpression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	 --filterName "my_snp_filter" \
    	-o $result/filtered_snps.vcf 

java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R $reference \
    	-V $result/raw_indels.vcf \
    	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
    	--filterName "my_indels_filter" \
    	-o $result/filtered_indels.vcf


echo "splitting vcf done!"| mail -s "vcf split"  "azzaea@gmail.com"
