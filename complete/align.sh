#!/bin/bash

## this script should contain options for doing alignment. Either bwa or novoalign

# the output should be saved in align_res, and specify the tool used for alginment and resources!
set -x
if [ $# -lt 10 ]; then
   echo -e "$0: error in calling the script, revise the arguments!" |  mail -s "accreditation pipeline" azzaea@gmail.com
   exit
fi
read1=$1
read2=$2
bwa_index=$3
novoalign_index=$4
rgheader=$5
align_res=$6
samplename=$7
reports=$8
email=$9
analysis=${10}
tool=${11}

set +x
module load bwa/0.7.15
module load novocraft/3.06.04
module load samtools/1.3.1
set -x

if [ `expr ${#tool}` -lt 1  ]; then
	echo -e "##################################################"
	echo -e "No tool is selected, doing both bwa and novoalign"
	echo -e "##################################################"
	mkdir -p $align_res/benchmarking
	for i in $(seq 1 1 8); do
		echo -e "##################################################"
		echo -e "Doing bwa and novoalign for #cores= \t $i"
		echo -e "##################################################"
		start=`date `
		bwa mem -M -t $i -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ $i -bS > ${align_res}/benchmarking/$samplename.aligned.bwa.$i.bam
		end=`date `
		echo -e "BWA_$i\t$start\t$end" >> $reports/timings.align
		./check_bam.sh ${align_res}/benchmarking/$samplename.aligned.bwa.$i.bam Alignment_bwa_$i $reports $samplename $email
		echo -e "################ BWA done for $i cores #############"
		start=`date `
		novoalign -c $i -d $novoalign_index -f $read1 $read2 -o SAM $rgheader | samtools view -@ $i -bS  > ${align_res}/benchmarking/$samplename.aligned.novoalign.$i.bam
		end=`date `
		echo -e "Novoalign_$i\t$start\t$end" >> $reports/timings.align
		./check_bam.sh ${align_res}/benchmarking/$samplename.aligned.novoalign.$i.bam Alignment_novoalign_$i $reports $samplename $email
		echo -e "############## NOVOALIGN done for $i cores ##########"
	done
else
	if [ $tool == 'bwa' ]; then
		echo -e "##################################################"
		echo -e "Choosing bwa for alignment for #cores= \t 4"
		echo -e "##################################################"
		start=`date `
		bwa mem -M -t 4 -R "$rgheader" $bwa_index $read1 $read2  | samtools view -@ 4 -bS > ${align_res}/$samplename.aligned.bwa.bam
		end=`date `
		echo -e "BWA\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/$samplename.aligned.bwa.bam Alignment_bwa $reports $samplename $email

	elif [ $tool == 'novoalign' ]; then
		echo -e "##################################################"
		echo -e "Choosing novoalign for alignment for #cores= \t 4"
		echo -e "##################################################"
		start=`date `
		novoalign -c 4 -d $novoalign_index -f $read1 $read2 -o SAM $rgheader | samtools view -@ 4 -bS  > ${align_res}/$samplename.aligned.novoalign.bam
		end=`date `
		echo -e "Novoalign\t$start\t$end" >> $reports/timings.$analysis
		./check_bam.sh ${align_res}/$samplename.aligned.novoalign.bam Alignment_novoalign $reports $samplename $email

	fi
fi
echo -e "##################################################"
echo -e "################# Alignment stage processing done #################"
echo -e "##################################################"

