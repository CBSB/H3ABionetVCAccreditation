#!/bin/bash

set -x
rawvcf=$1
gatkdir=
reference=$2
var_res=

stage=$2
reports=$3
gatkdir=$4
result=
email=$7
samplename=
############ Loading needed modules
set +x
module load gatk/3.6
module load picard-tools/2.6.0
module load bcftools/1.3.1
set -x

mkdir $var_res/filteration

########### The measurements
 
java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T SelectVariants \
        -R $reference \
	-V $rawvcf \
	-selectType SNP \
	-o $var_res/filteration/$stage.snps.vcf

bcftools query -f -H "%CHROM\t%ID\t%QUAL\t[%DP]\t%INFO/DP\t[%GQ]\t%INFO/AB\t[%QD]\t%INFO/MQ\t%INFO/FS\t%INFO/SOR\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\n" $input -o $var_res/filteration/tabulted_annotations.snps.$stage.txt

## also try this:
 bcftools query -H -f "%CHROM\t%ID\t%QUAL\t[%DP]\t%INFO/DP\t[%GQ]\t%INFO/AB\t%QD\t[%MQ]\t[%FS]\t[%SOR]\t[%MQRankSu]m\t[%ReadPosRankSum]\n" /home/assessment/results/CBSB_Khartoum/variants/CBSB_Khartoum.raw.calls.freebayes.vcf -o tabulated.tmp


java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R $reference \
        -V $rawvcf \
        -selectType INDEL \
        -o $var_res/filteration/$stage.indels.vcf



