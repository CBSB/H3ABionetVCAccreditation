#!/bin/bash

#This should contains options for indexing. 
# 1. samtools
# 2. picard
# 3. sambamba
# 4. novosort

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
module load samtools/1.3.1
module load picard-tools/2.6.0
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
		samtools sort -o $align_res/benchmarking/$samplename.sorted.samtools.$i.bam -@ $i $inputbam #-l=6
		end=`date `
		echo -e "samtools_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.sorted.samtools.$i.bam sorting_samtools_$i $reports $samplename $email
                echo -e "################ samtools sort done for $i cores #############"
		
		start=`date `
		sambamba sort -o $align_res/benchmarking/$samplename.sorted.sambamba.$i.bam -t $i $inputbam #-l=6
		end=`date `
		echo -e "sambamba_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.sorted.sambamba.$i.bam sorting_sambamba_$i $reports $samplename $email
		echo -e "################ sambamba sort done for $i cores #############"
		
		start=`date `
		java -jar $picarddir/picard.jar SortSam I=$inputbam O=$align_res/benchmarking/$samplename.sorted.picard.$i.bam SORT_ORDER=coordinate 
		end=`date `
		echo -e "picard_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.sorted.picard.$i.bam sorting_picard_$i $reports $samplename $email
		echo -e "################ picard sort done for $i cores #############"

		start=`date `
		novosort -c $i $inputbam > $align_res/benchmarking/$samplename.sorted.novosort.$i.bam #-6
		end=`date `
		echo -e "novosort_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.sorted.novosort.$i.bam sorting_novosort_$i $reports $samplename $email
		# novosort does a lot: sorting, merging, indexing, and marking duplicates at once
		# there is no multi-threading option in picard
		# picard does not allow control of compression levels. novosort compresses to 6 by default, while the others default to z-lib's default level 6. The default is put as a comment above
	done
else
	case $tool in
	samtools)
		start=`date `
		samtools sort -o $align_res/$samplename.sorted.samtools.bam -@ 4 $inputbam
		end=`date `
		echo -e "samtools\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}//$samplename.sorted.samtools.bam sorting_samtools $reports $samplename $email
		;;
	sambamba)
		start=`date `
		sambamba sort -o $align_res//$samplename.sorted.sambamba.bam -t 4 $inputbam
end=`date `
                echo -e "sambamba\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}//$samplename.sorted.sambamba.bam sorting_sambamba $reports $samplename $email
		;;
	picard)
		start=`date `
		java -jar $picarddir/picard.jar SortSam I=$inputbam O=$align_res//$samplename.sorted.picard.bam SORT_ORDER=coordinate
                end=`date `
                echo -e "picard\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}//$samplename.sorted.picard.bam sorting_picard $reports $samplename $email
		;;
	novosort)
		start=`date `
		novosort -c 4 $inputbam > $align_res//$samplename.sorted.novosort.bam
		end=`date `
                echo -e "novosort\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}//$samplename.sorted.novosort.bam sorting_novosort $reports $samplename $email

	esac
fi
