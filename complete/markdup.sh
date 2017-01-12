#!/bin/bash

#Checking deduplication options:
# 1. picard
# 2. sambamba
# 3. samblaster
# 4. samtools

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

echo -e "##################### Loading needed modules ###################"
set +x
module load picard-tools/2.6.0
module load samtools/1.3.1
module load sambamba/0.6.5
module load novocraft/3.06.04
picarddir="/usr/src/picard-tools/picard-tools-2.6.0"
set -x

if [ `expr ${#tool}` -lt 1  ]; then
	echo -e "##################################################"
        echo -e "No tool is selected, doing samtools, picard, novosort & sambamba"
        echo -e "##################################################"
        if [ ! -d "$align_res/benchmarking" ]; then
                 mkdir -p $align_res/benchmarking
        fi
        for i in $(seq 1 1 8); do
                echo -e "##################################################"
                echo -e "Doing sorting for #cores= \t $i"
                echo -e "##################################################"

                start=`date `
		java -Xmx10g -XX:-UseGCOverheadLimit -jar $picarddir/picard.jar MarkDuplicates  I=$inputbam  O=$align_res/benchmarking/$samplename.dedup.picard.$i.bam  M=$reports/$samplename.dedup.picard.$i.txt
		end=`date `
                echo -e "picard_$i\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}/benchmarking/$samplename.dedup.picard.$i.bam dedup_picard_$i $reports $samplename $email
                echo -e "################ picard markduplicates done for $i cores #############"

                start=`date `
                samtools rmdup  $inputbam  $align_res/benchmarking/$samplename.dedup.samtools.$i.bam 
                end=`date `
                echo -e "samtools_$i\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}/benchmarking/$samplename.dedup.samtools.$i.bam dedup_samtools_$i $reports $samplename $email
                echo -e "################ samtools markduplicates done for $i cores #############"

                start=`date `
                sambamba markdup -t $i  $inputbam  $align_res/benchmarking/$samplename.dedup.sambamba.$i.bam 
                end=`date `
                echo -e "sambamba_$i\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}/benchmarking/$samplename.dedup.sambamba.$i.bam dedup_sambamba_$i $reports $samplename $email
                echo -e "################ sambamba markduplicates done for $i cores #############"

		############# This is how far I got: 11/1/2017######## :)
	        start=`date `
                novosort markdup -t $i  $inputbam  $align_res/benchmarking/$samplename.dedup.sambamba.$i.bam 
                end=`date `
                echo -e "sambamba_$i\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}/benchmarking/$samplename.dedup.sambamba.$i.bam dedup_sambamba_$i $reports $samplename $email
                echo -e "################ sambamba markduplicates done for $i cores #############"


		#samtools: "Remove potential PCR duplicates: if multiple read pairs have identical external coordinates, only retain the pair with highest mapping quality. In the paired-end mode, this command ONLY works with FR orientation and requires ISIZE is correctly set. IT does not work for unpaired reads (e.g. two ends mapped to different chromosomes or orphan reads)."
	       # picard can give you the option to remove or mark
       	       # sambabmab: "Marks (by default) or removes duplicate reads. For determining whether a read is a duplicate or not, the same criteria as in Picard are used." It also allows specifying the number of threads to use. It also gives control to the compression level of the resulting file. However, "External sort is not implemented. Thus, memory consumption grows by 2Gb per each 100M reads. Check that you have enough RAM before running the tool." 
