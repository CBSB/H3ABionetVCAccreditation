#!/bin/bash

######### Removing mapping artifacts (clipping end of reference, missing mates, reads splitting into adapters...), sorting reads by mapping position, adding read group info if needed and adding cigar information

set -x
if [ $# -lt 9 ]; then
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
tool=${10}


set +x
module load picard-tools/2.6.0
picarddir="/usr/src/picard-tools/picard-tools-2.6.0"
set -x


mkdir -p $results/$samplename_ubam
cd $results/$samplename_ubam

java -Xmx8G -jar $picarddir/picard.jar FastqToSam FASTQ=$read1 FASTQ2=$read2 OUTPUT=unmappedbam.sam READ_GROUP_NAME=$rgheader SAMPLE_NAME=$samplename
# According to the broad's website the inputs need to be queryname sorted first, and this is doen automatically by default by the tool

java -Xmx8G -jar $picarddir/picard.jar MarkIlluminaAdapters I=unmappedbam.sam O=markilluminaadapters.bam M=markilluminaadapters_metrics.txt TMP_DIR=$results/tmp #O and M are the outputs from this command!
# Marking adapter sequences to minimize thier effect on the alignment and map unmappable reads. This bascially marks read pairs with XT tag to indicate if it is adapter sequence. A tmp_dir was used just to highlight that it is always a key option when working with large files


# A piped option would be preferred for the 3 steps that follow should more samples be involved (saves memory & speed) --> Something about different performance though, so better if one checks!

##### Unpiped processing:
java -Xmx8G -jar $picarddir/picard.jar SamToFastq I=markilluminaadapters.bam FASTQ=read1_cleaned.fq SECOND_END_FASTQ=read2_cleaned.fq CLIPPING_ATTRIBUTE=XT  CLIPPING_ACTION=2 NON_PF=true 
# take read identifier, read sequences and base scores to create a sanger fastq file. setting CLIPPING_ATTRIBUTE=XT and CLIPPING_ACTION=2, SamToFastq changes the quality scores of bases marked by XT to two (which is low quality on phred scale). This  effectively removes the adapter portion of sequences from contributing to downstream read alignment and alignment scoring metrics >> the resulting fastq file hass all meta data: read group, alignment, flags and tags in addition to the read query name, read sequences and read base quality scores

module load bwa/0.7.10
bwa mem -M -t 12 $reference/human read1_cleaned.fq read2_cleaned.fq > a.alignedbam.sam
#the alignment file has automatically generated read group info


java -Xmx16G -jar $picarddir/picard.jar MergeBamAlignment R=$reference/ucsc.hg19.fasta UNMAPPED_BAM=unmappedbam.sam ALIGNED_BAM=a.alignedbam.sam O=cleaned_bam.bam ALIGNER_PROPER_PAIR_FLAGS=false CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=true CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=BestMapq ATTRIBUTES_TO_RETAIN=XS UNMAP_CONTAMINANT_READS=true TMP_DIR=$results/tmp


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

#### Piped process:

set -o pipefail # stop the pipeline if any step reported an error

java -Xmx8G -jar $picarddir/picard.jar SamToFastq I=markilluminaadapters.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT  CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | \
bwa mem -M -t 12 -p $reference/human /dev/stdin| \
java -Xmx16G -jar $picarddir/picard.jar MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=unmappedbam.sam O=cleaned_bam_piped.bam R=$reference/ucsc.hg19.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=true CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=BestMapq ATTRIBUTES_TO_RETAIN=XS ALIGNER_PROPER_PAIR_FLAGS=false UNMAP_CONTAMINANT_READS=true TMP_DIR=$results/tmp 

##### Was it useful
# To finalize this stage, rememeber to sort and index your default bam file from the p1 stage, so that you can see how effective these strategies hase been:

cp  .
module load novocraft/3.02

novosort -c 12 --compression 5 --forcesort --index ../p2_alignment/default/a.default.0.bam -o a.default.bam

# operation distributed among 12 cores, forcing sorting, indexing the output while setting the compression level to 5 to match picard's default




echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
