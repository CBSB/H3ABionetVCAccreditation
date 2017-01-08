#!/bin/bash

# This file shoule be for bqsr and vc using:
# 1. GATK
# 2. Freebayes

module load gatk/3.6
gatkdir="/usr/src/gatk/gatk-3.6/"
###################################### Base recalibration stage: 
#-knownSites are verified here - from https://software.broadinstitute.org/gatk/guide/article?id=1247
java -Xmx10g -XX:-UseGCOverheadLimit -jar $gatkdir/GenomeAnalysisTK.jar \
	-T BaseRecalibrator\
	-R $reference \
	-I ${align_res}/$samplename.dedup.bam \
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
 
: <<'comment_PrintReads' 
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
comment_PrintReads 


############################################## Variant calling stage:
        # -dbsnp here is correct as per: https://software.broadinstitute.org/gatk/guide/article?id=1247
        java -jar $gatkdir/GenomeAnalysisTK.jar \
                -T HaplotypeCaller\
                -R $reference\
                -I  ${align_res}/$samplename.dedup.bam \
                -BQSR $reports/${samplename}.recal.table\
                --emitRefConfidence GVCF \
                --dbsnp $dbsnp138 \
                -o $vars/$samplename.raw.snps.indels.g.vcf\
                -nct 4

#REMEBER:  you don't need to gvcf the output and then convert

 exit_code=$? 
       if [ ! $exit_code -eq 0 ]; then 
               echo 'HaplotypeCaller did NOT work' 
               exit 
       fi 
       if [ ! -s $vars/$samplename.raw.snps.indels.g.vcf ]; then 
               echo 'Raw vcf file empty!' 
               exit 
       fi 

