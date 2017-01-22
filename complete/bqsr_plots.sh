#!/bin/bash

projectdir=/home/CBSB_UofK.Chr1_50X.set1
referencedir=${projectdir}/data/genome/hg19/
reference=${referencedir}/ucsc.hg19.fasta
align_res=/home/CBSB_UofK.Chr1_50X.set1/results_2/CBSB_Khartoum/align
targeted=/home/CBSB_UofK.Chr1_50X.set1/data/TruSeq_exome_targeted_regions.hg19.bed
dbsnp138=/home/CBSB_UofK.Chr1_50X.set1/data/genome/hg19//dbsnp_138.hg19.vcf
Mills=/home/CBSB_UofK.Chr1_50X.set1/data/genome/hg19//Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
TG_1000=/home/CBSB_UofK.Chr1_50X.set1/data/genome/hg19//1000G_phase1.indels.hg19.sites.vcf
reports=/home/CBSB_UofK.Chr1_50X.set1/results_2/tools_reports

module load R/3.2.3
module load gatk/3.6
gatkdir="/usr/src/gatk/gatk-3.6/"
:<<'comment'
java -jar $gatkdir/GenomeAnalysisTK.jar \
		-T BaseRecalibrator\
		-R $reference\
		-I $align_res/CBSB_Khartoum.dedup.picard.bam\
		-L $targeted\
		-knownSites $dbsnp138\
	 	-knownSites $Mills\
		-knownSites $TG_1000\
		-o $reports/CBSB_Khartoum.targeted.recal.table\
		-nct 4


java -jar $gatkdir/GenomeAnalysisTK.jar \
                -T BaseRecalibrator\
                -R $reference\
               -I $align_res/CBSB_Khartoum.dedup.picard.bam\
		 -L $targeted\
                -knownSites $dbsnp138\
                -knownSites $Mills\
                -knownSites $TG_1000\
		-BQSR $reports/CBSB_Khartoum.targeted.recal.table\
		-o $reports/CBSB_Khartoum.targeted.recal.table.after -nct 4
comment

java -jar $gatkdir/GenomeAnalysisTK.jar \
		-T AnalyzeCovariates\
		-R  $reference\
		-L $targeted\
		-before $reports/CBSB_Khartoum.targeted.recal.table\
		-after $reports/CBSB_Khartoum.targeted.recal.table.after\
		-plots $reports/CBSB_Khartoum.recalibration_plots.pdf 
exitcode=$?

echo "####################################################################"
echo "Program exit status $exitcode"



echo "####################################################################"
echo "Program exit status $exitcode"

