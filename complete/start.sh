#!/bin/bash

############################################# Analysis variables:
projectdir=/home/aeahmed/assessment
result=${projectdir}/results
samplename=
read1=
read2=
rgheader="@RG\tID:Set1\tSM:CBSB_Khartoum\tPL:illumina\tPU:synthetic\tLB:synthetic\tDT:2016-12-12"
email=

referencedir=${projectdir}/genome/hg19/
bwa_index=${referencedir}/ucsc.hg19.fasta
novoalign_index=
reference=${referencedir}/ucsc.hg19.fasta
dbsnp129=${referencedir}/dbsnp_138.hg19.excluding_sites_after_129.vcf
dbsnp138=${referencedir}/dbsnp_138.hg19.vcf
Mills=${referencedir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
TG_1000Gindels=${referencedir}/1000G_phase1.indels.hg19.sites.vcf
TG_1000Gsnps=${referencedir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf

targeted=chr1

################################# These are output directories. Let's create them!

############################################# Load necessary modules
module load bwa/0.7.15
module load samtools
module load gatk/3.6
module load picard-tools/2.6.0

############################################# Specifying path for java tools:
picarddir="/usr/src/picard-tools/picard-tools-2.6.0" #similar in all
gatkdir="/usr/src/gatk/gatk-3.6/"

set -x

############################################## These are output directories. Let's create them!
rm $result/raw_variants.txt

while read line ; do

	align_res=${result}/${samplename}/align
	vars=${result}/${samplename}/variants
	reports=$result/reports

	mkdir -p ${align_res}
	mkdir -p $reports
	mkdir -p $vars

	###################################################################### Actual start of the pipeline
	
	./aligh.sh $read1 $read2 $index $rgheader

	# bwa mem -M -t 4 -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ 4 -bS > ${align_res}/$samplename.bam
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

	samtools sort -o $align_res/$samplename.sorted.bam -@ 4 ${align_res}/$samplename.bam

	################################################################### Now, do marking duplicates
	java -Xmx10g -XX:-UseGCOverheadLimit -jar $picarddir/picard.jar MarkDuplicates \
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

	###################################### Base recalibration stage: 
	#-knownSites are verified here - from https://software.broadinstitute.org/gatk/guide/article?id=1247
	java -Xmx10g -XX:-UseGCOverheadLimit -jar $gatkdir/GenomeAnalysisTK.jar \
	############################################################################ Base recalibration stage:
	java -jar $gatkdir/GenomeAnalysisTK.jar \
		-T BaseRecalibrator\
		-R $reference \
		-I ${align_res}/$samplename.dedup.bam \
		-L $targeted\
		-knownSites $dbsnp138\
	 	-knownSites $Mills\
		-knownSites $TG_1000Gindels\
		-o $reports/${samplename}.recal.table\
		-nct 4
		-knownSites $dbsnp\
	 	-knownSites $TG_1000G\
		-o $reports/${samplename}.recal.table

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


## You need to add gatk VariantEval with known sites: $dbsnp129

java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T VariantEval\
        -R $reference\
	-o $result/output.eval.grp\
	--eval:gatk  $result/joint_called.vcf\
	-D $dbsnp129\
	-noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary


echo "Pipeline (vanialla gatk) run finished! yeeeey!"| mail -s "accreditation pipeline" "azzaea@gmail.com"
