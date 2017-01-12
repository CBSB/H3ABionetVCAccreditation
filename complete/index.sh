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
		echo -e "Doing indexting for #cores= \t $i"
		echo -e "##################################################"
		
		start=`date `
		samtools index $inputbam > $align_res/benchmarking/$samplename.indexed.samtools.$i.bam  
		end=`date `
		echo -e "samtools_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.indexed.samtools.$i.bam indexing_samtools_$i $reports $samplename $email
                echo -e "################ samtools index done for $i cores #############"
		
		start=`date `
		sambamba index -t $i $inputbam > $align_res/benchmarking/$samplename.indexed.sambamba.$i.bam  #-l=6
		end=`date `
		echo -e "sambamba_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.indexed.sambamba.$i.bam indexing_sambamba_$i $reports $samplename $email
		echo -e "################ sambamba index done for $i cores #############"
		
		start=`date `
		java -jar $picarddir/picard.jar BuildBamIndex I=$inputbam O=$align_res/benchmarking/$samplename.indexed.picard.$i.bam 
		end=`date `
		echo -e "picard_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.indexed.picard.$i.bam indexing_picard_$i $reports $samplename $email
		echo -e "################ picard index done for $i cores #############"

		start=`date `
		novosort -c $i $inputbam -i -o $align_res/benchmarking/$samplename.indexed.novosort.$i.bam #-6
		end=`date `
		echo -e "novosort_$i\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.indexed.novosort.$i.bam sorting_novosort_$i $reports $samplename $email
		# samtools index is not  multi-threaded tool!
		# novosort does a lot: sorting, merging, indexing, and marking duplicates at once
		# there is no multi-threading option in picard
	done
else
	case $tool in
	samtools)
		start=`date `
		samtools index $inputbam > $align_res/$samplename.indexed.samtools.bam
		end=`date `
		echo -e "samtools\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/benchmarking/$samplename.indexed.samtools.bam indexing_samtools $reports $samplename $email
		;;
	sambamba)
		start=`date `
		sambamba index -t 4 $inputbam > $align_res/benchmarking/$samplename.indexed.sambamba.bam	
		end=`date `
                echo -e "sambamba\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}/benchmarking/$samplename.indexed.sambamba.bam indexing_sambamba $reports $samplename $email
		;;
	picard)
		start=`date `
		java -jar $picarddir/picard.jar BuildBamIndex I=$inputbam O=$align_res/$samplename.indexed.picard.bam 
                end=`date `
                echo -e "picard\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}/$samplename.indexed.picard.bam indexing_picard $reports $samplename $email
		;;
	novosort)
		start=`date `
		novosort -c 4 $inputbam -i -o $align_res/$samplename.indexed.novosort.bam #-6
		end=`date `
                echo -e "novosort\t$start\t$end" >> $reports/timings.$analysis
                ./check_bam.sh ${align_res}/$samplename.indexed.novosort.bam sorting_novosort $reports $samplename $email

	esac
fi
