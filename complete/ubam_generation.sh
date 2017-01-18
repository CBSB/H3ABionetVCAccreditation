#!/bin/bash

# Removing mapping artifacts (clipping end of reference, missing mates, reads splitting into adapters...), sorting reads by mapping position, adding read group info if needed and adding cigar information

set -x
if [ $# -lt 8 ]; then
   echo -e "$0: error in calling the script, revise the arguments!" |  mail -s "accreditation pipeline" azzaea@gmail.com
   exit
fi

picarddir=$1
read1=$2
read2=$3
align_res=$4
samplename=$5
results=$6
reports=$7
bwa_index=$8

mkdir $results/tmp

set +x
module load picard-tools/2.6.0
module load bwa/0.7.15
set -x

java -Xmx8G -jar $picarddir/picard.jar FastqToSam FASTQ=$read1 FASTQ2=$read2 OUTPUT=$align_res/unmappedbam.sam\
	READ_GROUP_NAME="Set1" SAMPLE_NAME=$samplename PLATFORM="illumina" PLATFORM_UNIT="synthetic" LIBRARY_NAME="synthetic" RUN_DATE="2016-12-12"

# According to the broad's website the inputs need to be queryname sorted first, and this is doen automatically by default by the tool

java -Xmx8G -jar $picarddir/picard.jar MarkIlluminaAdapters I=$align_res/unmappedbam.sam O=$align_res/markilluminaadapters.bam M=$reports/markilluminaadapters_metrics.txt TMP_DIR=$results/tmp #O and M are the outputs from this command!
# Marking adapter sequences to minimize thier effect on the alignment and map unmappable reads. This bascially marks read pairs with XT tag to indicate if it is adapter sequence. A tmp_dir was used just to highlight that it is always a key option when working with large files

set -o pipefail # stop the pipeline if any step reported an error

java -Xmx8G -jar $picarddir/picard.jar SamToFastq I=$align_res/markilluminaadapters.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT  CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | \
bwa mem -M -t 4 -p $bwa_index /dev/stdin| \
java -Xmx16G -jar $picarddir/picard.jar MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=$align_res/unmappedbam.sam O=$align_res/$samplename.sorted.bam R=$reference/ucsc.hg19.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=true CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=BestMapq ATTRIBUTES_TO_RETAIN=XS ALIGNER_PROPER_PAIR_FLAGS=false UNMAP_CONTAMINANT_READS=true TMP_DIR=$results/tmp 


# A piped option would be preferred for the 3 steps that follow should more samples be involved (saves memory & speed) --> Something about different performance though, so better if one checks!

#SamToFastq options:
# take read identifier, read sequences and base scores to create a sanger fastq file. setting CLIPPING_ATTRIBUTE=XT and CLIPPING_ACTION=2, SamToFastq changes the quality scores of bases marked by XT to two (which is low quality on phred scale). This  effectively removes the adapter portion of sequences from contributing to downstream read alignment and alignment scoring metrics >> the resulting fastq file hass all meta data: read group, alignment, flags and tags in addition to the read query name, read sequences and read base quality scores

#MergeBamAlignment options:
#ALIGNER_PROPER_PAIR_FLAGS=false\ # default: reassess and reassign the aligners' proper pair flags
#CREATE_INDEX=true\ # by default the output is coordinate sorted, and this option creates the index (bai) file
#ADD_MATE_CIGAR=true\ # by default to add MC tag
#CLIP_ADAPTERS=true\ # done automatically in the data from uiuc/igb, but should not hurt to verify
#CLIP_OVERLAPPING_READS=true\ # by default soft-clips ends so mates do not overlap
#INCLUDE_SECONDARY_ALIGNMENTS=true\ # secondary aligning reads may contain interesting info
#MAX_INSERTIONS_OR_DELETIONS=-1\ #changed to allow any number of insertions or deletions
#PRIMARY_ALIGNMENT_STRATEGY=BestMapq\
#ATTRIBUTES_TO_RETAIN=XS\ #carryover the XS tag from the alignment, which reports BWA-MEM's suboptimal alignment scores. The XS tag in unnecessary GATK bp, 
#UNMAP_CONTAMINANT_READS=true\ #detects reads from foreign organisms if found on the read
#TMP_DIR=$results/tmp  

