#!/bin/bash

############################################# Analysis variables:
projectdir=/home/aeahmed/assessment
referencedir=${projectdir}/genome/hg19/
bwa_index=${referencedir}/ucsc.hg19.fasta
reference=${referencedir}/ucsc.hg19.fasta
dbsnp129=${referencedir}/dbsnp_138.hg19.excluding_sites_after_129.vcf
dbsnp138=${referencedir}/dbsnp_138.hg19.vcf
Mills=${referencedir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
TG_1000Gindels=${referencedir}/1000G_phase1.indels.hg19.sites.vcf
TG_1000Gsnps=${referencedir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf

result=${projectdir}/results/run2

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
######################################## Preparing for alignment
rm $result/raw_variants.txt

reports=$result/reports
mkdir -p $reports

dstat -t -c --disk-util --top-cpu --top-io --top-mem -s --output $reports/profile.complete_pipeline.log 30 600 &
echo PROCESS START_TIME END_TIME >$reports/timings

while read line ; do
	samplename=$(echo "$line" | cut -d ' ' -f1 )
	read1=$(echo "$line" | cut -d ' ' -f2)
	read2=$(echo "$line" | cut -d ' ' -f3)

	align_res=${result}/${samplename}/align
	vars=${result}/${samplename}/variants

	mkdir -p ${align_res}
	mkdir -p $vars

	######################################### Actual start of the pipeline
	rgheader="@RG\tID:Set1\tSM:CBSB_Khartoum\tPL:illumina\tPU:synthetic\tLB:synthetic\tDT:2016-12-12"
	
	start=`date`	
	bwa mem -M -t 4 -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ 4 -bS > ${align_res}/$samplename.bam
	exit_code=$?
	end=`date`
	
	echo 'BWA' $start $end >> $reports/timings

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
	
	start=`date`
	samtools sort -o $align_res/$samplename.sorted.bam -@ 4 ${align_res}/$samplename.bam
	end=`date`
	echo 'sorting'  $start $end >> $reports/timings

	################################################################### Now, do marking duplicates
	start=`date`
	java -Xmx10g -XX:-UseGCOverheadLimit -jar $picarddir/picard.jar MarkDuplicates \
	      I=$align_res/$samplename.sorted.bam\
	      O=$align_res/$samplename.dedup.bam \
	      M=$reports/$samplename.dedup.txt
	exit_code=$?
	end=`date`
	
	echo 'deduplication'  $start $end >> $reports/timings

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

	start=`date`
	samtools index ${align_res}/$samplename.dedup.bam
	end=`date`
	echo 'indexing'  $start $end >> $reports/timings

	###################################### Base recalibration stage: 
	#-knownSites are verified here - from https://software.broadinstitute.org/gatk/guide/article?id=1247
	start=`date`
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
	end=`date`
	echo 'bqsr'  $start $end >> $reports/timings

	if [ ! $exit_code -eq 0 ]; then
		echo 'Base Recalibrator did NOT work'
		exit
	fi

: <<'comment_PrintReads'
	start=`date`
	java -jar $gatkdir/GenomeAnalysisTK.jar\
		-T PrintReads\
		-R $reference\
		-I  ${align_res}/$samplename.dedup.bam \
		-L $targeted \
		-BQSR $reports/${samplename}.recal.table\
		-o $align_res/${samplename}.recal.bam
	exit_code=$?
	end=`date`
	echo 'PrintReads' $start $end >> $reports/timings 
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

	start=`date`
	java -jar $gatkdir/GenomeAnalysisTK.jar \
                -T HaplotypeCaller\
                -R $reference\
		-I  ${align_res}/$samplename.dedup.bam \
		-BQSR $reports/${samplename}.recal.table\
		--emitRefConfidence GVCF \
		--dbsnp $dbsnp138 \
		-o $vars/$samplename.raw.snps.indels.g.vcf\
		-nct 4
       exit_code=$?
       end=`date`
       echo 'hc'  $start $end >> $reports/timings

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

start=`date`
java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs\
        -R $reference\
	$(cat $result/raw_variants.txt)	\
	-o $result/joint_called.vcf
exit_code=$?
end=`date`
echo 'genotypeGVCFs'  $start $end >> $reports/timings

if [ ! $exit_code -eq 0 ]; then
	echo 'GenotypeGVCFs did NOT work'
9        exit
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
