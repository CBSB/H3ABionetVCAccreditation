#!/bin/bash

############################################# Analysis variables:
set -x
echo -e "This analysis is starting at \t `date`"
projectdir=/home/assessment
email=azzaea@gmail.com

result=${projectdir}/results
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

targeted=chr1

################################### Specific analysis options and tools
analysis=sort
align_tool=bwa
sort_tool=
dedup_tool=

################################### Output directories
qc_res=$result/${samplename}/qc
align_res=${result}/$samplename/align
vars_res=${result}/$samplename/variants
delivery=$result/delivery 
reports=$result/reports

mkdir -p $qc_res ${align_res} $vars_res $delivery $reports

############################################# Actual start of the pipeline
dstat -t -c --disk-util --top-cpu --top-io --top-mem -s --output $reports/profile.complete_pipeline_${analysis}.log 1 18000 &
echo -e "PROCESS\tSTART_TIME\tEND_TIME" >$reports/timings.$analysis

./align.sh $read1 $read2 $bwa_index $novoalign_index $rgheader $align_res $samplename $reports $email $align_tool

if [ $analysis == "align" ];then
	echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" | mail -s "accreditation pipeline" $email
	exit
fi

./sort.sh ${align_res}/${samplename}.aligned.${align_tool}.bam $align_res $samplename $reports $email $analysis $sort_tool

if [ $analysis == "sort" ];then
	echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" | mail -s "accreditation pipeline" $email
	exit
fi

################################################################### Now, do marking duplicates
./markdup.sh $align_res/$samplename.sorted.$sort_tool.bam $align_res $samplename $reports $email $analysis $dedup_tool

./index.sh ${align_res}/$samplename.dedup.bam

./bqvc.sh


## You need to add gatk VariantEval with known sites: $dbsnp129

java -jar $gatkdir/GenomeAnalysisTK.jar \
        -T VariantEval\
        -R $reference\
	-o $result/output.eval.grp\
	--eval:gatk  $result/joint_called.vcf\
	-D $dbsnp129\
	-noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary


echo "Pipeline (vanialla gatk) run finished! yeeeey!"| mail -s "accreditation pipeline" $email
