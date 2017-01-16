#!/bin/bash

# This file shoule be for bqsr and vc using:
# 1. GATK
# 2. Freebayes

echo -e "################## Parsing command line arguments ##############"
set -x
if [ $# -lt 6 ]; then
   echo -e "$0: error in calling the script, revise the arguments!" |  mail -s "accreditation pipeline" azzaea@gmail.com
   exit
fi
inputbam=$1
align_res=$2
samplename=$3
reports=$4
email=$5
analysis=$6
tool=$7
reference=
targeted=
dbsnp138=
Mills=
TG_1000Gindels=

echo -e "##################### Loading needed modules ###################"
set +x

module load gatk/3.6
Module load freebayes/1.1.0 
gatkdir="/usr/src/gatk/gatk-3.6/"
###################################### Base recalibration stage: 
java -Xmx10g -XX:-UseGCOverheadLimit -jar $gatkdir/GenomeAnalysisTK.jar \
	-T BaseRecalibrator\
	-R $reference \
	-I $inputbam \
	-L $targeted\
	-knownSites $dbsnp138\
	-knownSites $Mills\
	-knownSites $TG_1000Gindels\
	-o $reports/${samplename}.recal.table\
	-nct 4

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


############################################## Variant calling stage:
java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T HaplotypeCaller\
	-R $reference\
	-I  ${align_res}/$samplename.dedup.bam \
	-BQSR $reports/${samplename}.recal.table\
	--dbsnp $dbsnp138 \
	-o $vars/$samplename.raw.calls.haplotypecaller.vcf\
	-nct 4


 exit_code=$? 
       if [ ! $exit_code -eq 0 ]; then 
               echo 'HaplotypeCaller did NOT work' 
               exit 
       fi 
       if [ ! -s $vars/$samplename.raw.calls.haplotypecaller.vcf ]; then 
               echo 'Raw vcf file empty!' 
               exit 
       fi 


################################## Freebayes
freebayes -f $reference ${align_res}/$samplename.dedup.bam > $vars/$samplename.raw.calls.freebayes.vcf

