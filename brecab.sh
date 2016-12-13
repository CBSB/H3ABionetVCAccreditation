#!/bin/bash

gatkdir="/usr/src/gatk/gatk-3.6/"
referencedir="/home/classrooms/Workshop_low_pass/ref/"
bams="/home/classrooms/Workshop_low_pass/align/"

module load samtools/1.3.1
samtools index $bams/marked_duplicates.bam
module purge

module load gatk/3.6
if [ $done ]; then
java -jar $gatkdir/GenomeAnalysisTK.jar \
		-T BaseRecalibrator\
		-R $referencedir/human_g1k_v37_chr20.fa\
		-I $bams/marked_duplicates.bam\
		-L 20\
		-knownSites $referencedir/dbsnp_135.b37.chr20.smallregion.vcf\
	 	-knownSites $referencedir/1kg.pilot_release.merged.indels.sites.hg19.chr20.vcf \
		-o $bams/recal.table.default

java -jar $gatkdir/GenomeAnalysisTK.jar\
		-T PrintReads\
		-R $referencedir/human_g1k_v37_chr20.fa\
		-I  $bams/marked_duplicates.bam\
		-L 20 \
		-BQSR $bams/recal.table.default\
		-o $bams/recal.default.bam

#: <<'comment'
java -jar $gatkdir/GenomeAnalysisTK.jar \
                -T BaseRecalibrator\
                -R $referencedir/human_g1k_v37_chr20.fa\
                -I $bams/marked_duplicates.bam\
                -L 20\
                -knownSites $referencedir/dbsnp_135.b37.chr20.smallregion.vcf\
                -knownSites $referencedir/1kg.pilot_release.merged.indels.sites.hg19.chr20.vcf \
                -o $bams/recal.table.default.after\
		-BQSR $bams/recal.table.default
exitcode=$?

echo "####################################################################"
echo "Program exit status $exitcode"

fi

module load R/3.2.3
java -jar $gatkdir/GenomeAnalysisTK.jar \
		-T AnalyzeCovariates\
		-R $referencedir/human_g1k_v37_chr20.fa\
		-before $bams/recal.table.default\
		-after $bams/recal.table.default.after\
		-plots $bams/recalibration_plots.default.pdf
#comment

exitcode=$?

echo "####################################################################"
echo "Program exit status $exitcode"

