#!/bin/bash

projectdir=/home/aeahmed/assessment
referencedir=${projectdir}/genome/hg19/
bwa_index=${referencedir}/ucsc.hg19.fasta
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

######################### The measurements
:<<'comment'
java -jar $picarddir/picard.jar  CollectAlignmentSummaryMetrics R=$reference I=$inputbam O=$result/Picard_CollectAlignmentSummaryMetrics.txt
# produces metrics relating to the overall alignment of reads withing a SAM/BAM. ie it works per bam file
	# The reference file needs a companion dictionary here! it means, use the GATK reference, not the bwa indexed file

java -jar $picarddir/picard.jar CollectWgsMetrics I=$inputbam O=$result/Picard_CollectWGsMetrics.txt R=$reference
# produces metrics relating to reads that pass base and mapping quality filters and coverage (read depth) levels (user defined)
# to use this tool with WES, you need to give the coordinates of your genomic region!!

java -jar $picarddir/picard.jar CollectInsertSizeMetrics I=$inputbam O=$result/Picard_Insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5
# produces metrics for validating library construction including the insert size distribution and read 	orientation of paired end libraries
# setting the minimum percentage (M=0.5) is convineint for processing a small file
comment

java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T DiagnoseTargets \
        -R $reference \
	-I $inputbam \
	-L chr1 \
	-o $result/GATK-DiagnoseTargets.vcf
:<<'comment2'
java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T DepthOfCoverage \
        -R $reference \
	-I $inputbam \
	-L chr1 \
	-o $result/GATK-DOC-output
comment2
echo "collecting stats done!"| mail -s "bam checks"  "azzaea@gmail.com"
