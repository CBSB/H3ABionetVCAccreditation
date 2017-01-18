#!/bin/bash

set -x
rawvcf=$1
gatkdir=$2
reference=$3
var_res=$4
vcall_tool=$5
samplename=$6

############ Loading needed modules
set +x
module load gatk/3.6
module load picard-tools/2.6.0
module load bcftools/1.3.1
set -x

mkdir $var_res/filteration

############# Extracting snps and tabulating results
java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T SelectVariants \
        -R $reference \
	-V $rawvcf \
	-selectType SNP \
	-o $var_res/filteration/$samplename.$vcall_tool.snps.vcf

if [ $vcall_tool == gatk* ]; then
	bcftools query -H -f "%CHROM\t%ID\t%QUAL\t[%DP]\t%INFO/DP\t[%GQ]\t%QD\t%MQ\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\n" $var_res/filteration/$vcall_tool.snps.vcf -o $var_res/filteration/$samplename.$stage.tabulted_annotations.snps.txt
elif [ $vcall_tool == freebayes* ]; then
	bcftools query -H -f "%CHROM\t%ID\t%QUAL\t[%DP]\t%INFO/DP\t[%GQ]\tAB\t%QD\t%MQ\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\n" $var_res/filteration/$vcall_tool.snps.vcf -o $var_res/filteration/tabulted_annotations.snps.$stage.txt
fi

############ Extracting indels and tabulating results
java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T SelectVariants \
        -R $reference \
        -V $rawvcf \
        -selectType INDEL \
        -o $var_res/filteration/$samplename.$vcall_tool.indels.vcf

if [ $vcall_tool == gatk* ]; then
	bcftools query -H -f "%CHROM\t%ID\t%QUAL\t[%DP]\t%INFO/DP\t[%GQ]\t%QD\t%MQ\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\n" $var_res/filteration/$vcall_tool.snps.vcf -o $var_res/filteration/$samplename.$stage.tabulted_annotations.indels.txt
elif [ $vcall_tool == freebayes* ]; then

fi	
