#!/bin/bash

# This file shoule be for bqsr and vc using:
# 1. GATK
# 2. Freebayes

set -x
echo -e "################## Parsing command line arguments ##############"
set -x
if [ $# -lt 11 ]; then
   echo -e "$0: error in calling the script, revise the arguments!" |  mail -s "accreditation pipeline" azzaea@gmail.com
   exit
fi
inputbam=$1
align_res=$2
samplename=$3
reports=$4
email=$5
analysis=$6
reference=$7
targeted=$8
dbsnp138=$9
Mills=${10}
TG_1000Gindels=${11}
vars_res=${12}
dbsnp129=${13}
tool=${14}

echo -e "##################### Loading needed modules ###################"
set +x
module load gatk/3.6
module load freebayes/1.1.0 
gatkdir="/usr/src/gatk/gatk-3.6/"
set -x 

###################################### Base recalibration stage: 
start=`date `
java -Xmx10g -XX:-UseGCOverheadLimit -jar $gatkdir/GenomeAnalysisTK.jar \
	-T BaseRecalibrator\
	-R $reference \
	-I $inputbam \
	-L $targeted\
	-knownSites $dbsnp138\
	-knownSites $Mills\
	-knownSites $TG_1000Gindels\
	-o $reports/${samplename}.targeted.recal.table \
	 -nct 4
end=`date `
echo -e "BaseRecalibrator\t$start\t$end" >> $reports/timings.$analysis

start=`date `
java -jar $gatkdir/GenomeAnalysisTK.jar\
	-T PrintReads\
	-R $reference\
	-I $inputbam\
	-BQSR $reports/${samplename}.targeted.recal.table\
	-o $align_res/${samplename}.targeted.recal.bam
 end=`date `
 echo -e "PrintReads\t$start\t$end" >> $reports/timings.$analysis

./check_bam.sh $align_res/${samplename}.targeted.recal.bam bqsr_gatk $reports $samplename $email

############################################## Variant calling stage:
start=`date `
java -jar $gatkdir/GenomeAnalysisTK.jar \
	-T HaplotypeCaller\
	-R $reference\
	-I $inputbam\
	-BQSR $reports/${samplename}.targeted.recal.table\
	--dbsnp $dbsnp138 \
	-o $vars_res/$samplename.raw.calls.haplotypecaller.targeted.vcf\
	-L $targeted \
	-nct 4
end=`date `
echo -e "HaplotypeCaller\t$start\t$end" >> $reports/timings.$analysis

./check_vcf.sh $vars_res/$samplename.raw.calls.haplotypecaller.targeted.vcf gatk $reports $gatkdir $reference $dbsnp129 $email

./filter_vcf.sh $vars_res/$samplename.raw.calls.haplotypecaller.targeted.vcf $gatkdir $reference $vars_res gatk $samplename
################################## Freebayes
start=`date `
freebayes -t $targeted -= -f $reference $inputbam > $vars_res/$samplename.raw.calls.freebayes.targeted.vcf
end=`date `
echo -e "freebayes\t$start\t$end" >> $reports/timings.$analysis

./check_vcf.sh $vars_res/$samplename.raw.calls.freebayes.targeted.vcf freebayes $reports $gatkdir $reference $dbsnp129 $email

./filter_vcf.sh $vars_res/$samplename.raw.calls.freebayes.targeted.vcf $gatkdir $reference $vars_res freebayes $samplename


#################################### Validate both vcf files:
./compare_vcfs.sh  $vars_res/$samplename.raw.calls.haplotypecaller.targeted.vcf $vars_res/$samplename.raw.calls.freebayes.targeted.vcf
