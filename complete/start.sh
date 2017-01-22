#!/bin/bash

############################################# Analysis variables:
set -x
echo -e "This analysis is starting at \t `date`"
projectdir=/home/CBSB_UofK.Chr1_50X.set1
email=azzaea@gmail.com

result=${projectdir}/results_2
samplename=CBSB_Khartoum
read1=${projectdir}/data/reads/H3A_NextGen_assessment.Chr1_50X.set1_read1.fq
read2=${projectdir}/data/reads/H3A_NextGen_assessment.Chr1_50X.set1_read2.fq
rgheader="@RG\tID:Set1\tSM:CBSB_Khartoum\tPL:illumina\tPU:synthetic\tLB:synthetic\tDT:2016-12-12"

referencedir=${projectdir}/data/genome/hg19/
bwa_index=${referencedir}/ucsc.hg19.fasta
novoalign_index=${referencedir}/ucsc.hg19.nix
reference=${referencedir}/ucsc.hg19.fasta
dbsnp129=${referencedir}/dbsnp_138.hg19.excluding_sites_after_129.vcf
dbsnp138=${referencedir}/dbsnp_138.hg19.vcf
Mills=${referencedir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
TG_1000Gindels=${referencedir}/1000G_phase1.indels.hg19.sites.vcf
TG_1000Gsnps=${referencedir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf

targeted=$projectdir/data/TruSeq_exome_targeted_regions.hg19.bed
picarddir="/usr/src/picard-tools/picard-tools-2.6.0"
gatkdir="/usr/src/gatk/gatk-3.6/"

################################### Specific analysis options and tools
processing=normal #{normal | ubam}
analysis=vcall #{align | sort | dedup | index | nothing provokes the entire pipeline :)}
align_tool=bwa #{bwa | novoalign}
sort_tool=samtools #{samtools | sambamba| picard | novosort}
dedup_tool=picard #{picard | samtools | sambamba | novosort}
index_tool=samtools #{samtools | sambamba | picard | novosort}
vcall_tool= #{gatk | freebayes}

################################### Output directories
if [ $processing == "ubam" ]; then
	samplename=$samplename.ubam
fi
qc_res=$result/${samplename}/qc
align_res=${result}/$samplename/align
vars_res=${result}/$samplename/variants
reports=$result/tools_reports

mkdir -p $qc_res ${align_res} $vars_res $reports

############################################# Actual start of the pipeline
dstat -t -c --disk-util --top-cpu --top-io --top-mem -s --output $reports/profile.complete_pipeline_${analysis}.log 1 18000 &
dsprocess=$(ps -aux | grep dstat | awk '{ print $2}' |head -n 1)
echo -e "PROCESS\tSTART_TIME\tEND_TIME" >$reports/timings.$analysis

./fastqc.sh $read1 $read2 $qc_res

if [ $processing == "normal" ]; then
	./align.sh $read1 $read2 $bwa_index $novoalign_index $rgheader $align_res $samplename $reports $email $analysis $align_tool
	if [ $analysis == "align" ];then
		echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" | mail -s "accreditation pipeline" $email
		kill -9 $dsprocess
		exit
	fi
	
	./sort.sh ${align_res}/${samplename}.aligned.${align_tool}.bam $align_res $samplename $reports $email $analysis $sort_tool
	if [ $analysis == "sort" ];then
		echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" | mail -s "accreditation pipeline" $email
		kill -9 $dsprocess
		exit
	fi
elif [ $processing == "ubam" ]; then
	./ubam_generation.sh $picarddir $read1 $read2 $align_res $samplename $result $reports $bwa_index $reference 
fi

#############################

./coverage.sh $align_res/$samplename.sorted.bam $targeted $reports $samplename

./markdup.sh $align_res/$samplename.sorted.bam $align_res $samplename $reports $email $analysis $align_tool $dedup_tool

if [ $analysis == "dedup" ];then
	echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" | mail -s "accreditation pipeline" $email
	kill -9 $dsprocess
	exit
fi

./index.sh ${align_res}/$samplename.dedup.$dedup_tool.bam $align_res $samplename $reports $email $analysis $index_tool
if [ $analysis == "index" ];then
	echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" | mail -s "accreditation pipeline" $email
	kill -9 $dsprocess
	exit
fi

./bqvc.sh ${align_res}/$samplename.dedup.$dedup_tool.bam $align_res $samplename $reports $email $analysis $reference $targeted $dbsnp138 $Mills $TG_1000Gindels $vars_res $dbsnp129 $vcall_tool


kill -9 $dsprocess

echo "Pipeline (vanialla gatk) run finished! yeeeey!"| mail -s "accreditation pipeline" $email
