#!/bin/bash

########################################### Input files for my analysis
#samplename=HG00120.lowcoverage.chr20.smallregion
#read1=/home/classrooms/Workshop_low_pass/fastq/${samplename}_1.fastq 
#read2=/home/classrooms/Workshop_low_pass/fastq/${samplename}_2.fastq

############################################# Load necessary modules
module load bwa/0.7.15
module load samtools
module load gatk/3.6
module load picard-tools/2.6.0
picarddir="/usr/src/picard-tools/picard-tools-2.6.0" #similar in all
gatkdir="/usr/src/gatk/gatk-3.6/"

set -x
######################################## Preparing for alignment
targeted=20
referencedir=/home/classrooms/Workshop_low_pass/ref
bwa_index=${referencedir}/human_g1k_v37_chr20.fa
reference=${referencedir}/human_g1k_v37_chr20.fa
dbsnp=${referencedir}/dbsnp_135.b37.chr20.smallregion.vcf
TG_1000G=${referencedir}/1kg.pilot_release.merged.indels.sites.hg19.chr20.vcf

############################################## These are output directories. Let's create them!
result=/home/classrooms/Workshop_low_pass/results
rm $result/raw_variants.txt

while read line ; do
	samplename=$(echo "$line" | cut -d ' ' -f1 )
	read1=$(echo "$line" | cut -d ' ' -f2)
	read2=$(echo "$line" | cut -d ' ' -f3)

	align_res=${result}/${samplename}/align
	vars=${result}/${samplename}/variants
	reports=$result/reports

	mkdir -p  ${align_res}
	mkdir -p $reports
	mkdir -p $vars

	###################################################################### Actual start of the pipeline
	rgheader="@RG\tID:${samplename}\tPL:illumina\tPU:synthetic\tLB:synthetic\tDT:2016-7-1\tSM:${samplename}" 
	bwa mem -M -t 4 -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ 4 -bS > ${align_res}/$samplename.bam
	exit_code=$?
	if [ ! $exit_code -eq 0 ]; then
		echo 'Alignment did NOT work'
		exit
	fi

	if [ ! -s ${align_res}/$samplename.bam ]; then
		echo 'Alignment file empty!'
		exit
	fi

	numalign=$(samtools view -c ${align_res}/$samplename.bam)
	if [  $numalign -eq 0 ]; then
		echo 'Empty bam file'
		exit
	fi

	samtools flagstat ${align_res}/$samplename.bam >> $reports/${samplename}.summary.txt

	samtools sort -o $align_res/$samplename.sorted.bam -@ 2    ${align_res}/$samplename.bam

	################################################################### Now, do marking duplicates
	java -jar $picarddir/picard.jar MarkDuplicates \
	      I=$align_res/$samplename.sorted.bam\
	      O=$align_res/$samplename.dedup.bam \
	      M=$reports/$samplename.dedup.txt

	exit_code=$?
	if [ ! $exit_code -eq 0 ]; then
		echo 'Marking duplicates did NOT work'
		exit
	fi

	if [ ! -s ${align_res}/$samplename.dedup.bam ]; then
		echo 'Marking duplicates file empty!'
		exit
	fi

	numalign=$(samtools view -c ${align_res}/$samplename.dedup.bam)
	if [  $numalign -eq 0 ]; then
		echo 'Empty bam file'
		exit
	fi

	samtools index ${align_res}/$samplename.dedup.bam

	############################################################################ Base recalibration stage:
	java -jar $gatkdir/GenomeAnalysisTK.jar \
		-T BaseRecalibrator\
		-R $reference \
		-I ${align_res}/$samplename.dedup.bam \
		-L $targeted\
		-knownSites $dbsnp\
	 	-knownSites $TG_1000G\
		-o $reports/${samplename}.recal.table

	exit_code=$?
	if [ ! $exit_code -eq 0 ]; then
		echo 'Base Recalibrator did NOT work'
		exit
	fi

	java -jar $gatkdir/GenomeAnalysisTK.jar\
		-T PrintReads\
		-R $reference\
		-I  ${align_res}/$samplename.dedup.bam \
		-L $targeted \
		-BQSR $reports/${samplename}.recal.table\
		-o $align_res/${samplename}.recal.bam
	exit_code=$?
	if [ ! $exit_code -eq 0 ]; then
		echo 'Base Recalibrator did NOT work'
		exit
	fi
	if [ ! -s ${align_res}/$samplename.recal.bam ]; then
		echo 'Base recalibrated file empty!'
		exit
	fi

	numalign=$(samtools view -c ${align_res}/$samplename.recal.bam)
	if [  $numalign -eq 0 ]; then
		echo 'Empty bam file'
		exit
	fi

	########################################################################## Variant calling stage:
	java -jar $gatkdir/GenomeAnalysisTK.jar \
                -T HaplotypeCaller\
                -R $reference\
		-I $align_res/$samplename.recal.bam\
		--emitRefConfidence GVCF \
		--dbsnp $dbsnp \
		-o $vars/$samplename.raw.snps.indels.g.vcf
       exit_code=$?
       if [ ! $exit_code -eq 0 ]; then
               echo 'HaplotypeCaller did NOT work'
               exit
       fi
       if [ ! -s $vars/$samplename.raw.snps.indels.g.vcf ]; then
               echo 'Raw vcf file empty!'
               exit
       fi
	

echo " -V $vars/$samplename.raw.snps.indels.g.vcf " >> $result/raw_variants.txt

done < sampleinfo

java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs\
        -R $reference\
	$(cat $result/raw_variants.txt)	\
	-o $result/joint_called.vcf

exit_code=$?
if [ ! $exit_code -eq 0 ]; then
	echo 'GenotypeGVCFs did NOT work'
        exit
fi
if [ ! -s $result/joint_called.vcf ]; then
	echo 'joint calling file empty!'
        exit
fi

