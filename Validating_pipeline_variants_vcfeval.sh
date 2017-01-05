#!/bin/bash
#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/src/log.vcfeval.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/src/log.vcfeval.e 
#PBS -M aeahmed@illinois.edu
#PBS -m abe

######### Paths defintions:
set -x
output=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/results/run8
goldenFile=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/reads/synthetic.Oct2016.WES30x/sample.1/Merged/synthetic.Oct2016.WES30x_merged_golden.vcf
workflowFile=$output/delivery/jointVCFs/jointVCFcalled.vcf
targeted_region=/home/groups/hpcbio_shared/azza/TargetedRegions-Azza-has-permission/ZachUses_Illumina_truseq_exome_targeted_regions.hg19.chr.bed

########################### Preparatory stages:
set +x
module load tabix

set -x
bgzip -c $goldenFile > $goldenFile.gz
tabix -p vcf $goldenFile.gz

bgzip -c $workflowFile > $workflowFile.gz
tabix -p vcf $workflowFile.gz

########################### Preparing rtg files:
rtg=/home/groups/hpcbio_shared/azza/apps/rtg-tools-3.7/rtg
referencedir=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/
$rtg format -o $referencedir/HG19_GATKbundle2.8_noDecoys $referencedir/HG19_GATKbundle2.8_noDecoys.fa

########################### Comparison stage:

set -x
$rtg vcfeval \
 -t $referencedir/HG19_GATKbundle2.8_noDecoys \
 -b $goldenFile.gz -c $workflowFile.gz -o $output/variant_comparison_vcfeval\
 --squash-ploidy \
 --sample ALT,synthetic.Oct2016.WES30x
 --evaluation-regions=$targeted_region

$rtg rocplot $output/variant_comparison_vcfeval/weighted_roc.tsv.gz \
  --png=$output/variant_comparison_vcfeval/roc.png \
  --title="Variant Calling Pipeline Variants  VS  Synthetic reads variants" \
  --precision-sensitivity

$rtg vcfstats $output/variant_comparison_vcfeval/*.vcf.gz $goldenFile.gz $workflowFile.gz > $output/variant_comparison_vcfeval/comp_stats.txt
