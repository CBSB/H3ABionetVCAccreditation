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

result=${projectdir}/results

targeted=chr1

input=$result/joint_called.vcf

############################################# Load necessary modules
module load bcftools/1.3.1

set -x

which=both

bcftools query -f -c $which "%CHROM\t%ID\t%QUAL\t%INFO/QD" $input -o $result/annotations/QD
bcftools query -f -c $which "%CHROM\t%ID\t%QUAL\t%INFO/MQ" $input -o $result/annotations/MQ
bcftools query -f -c $which "%CHROM\t%ID\t%QUAL\t%INFO/FS" $input -o $result/annotations/FS
bcftools query -f -c $which "%CHROM\t%ID\t%QUAL\t%INFO/SOR" $input -o $result/annotations/SOR
bcftools query -f -c $which "%CHROM\t%ID\t%QUAL\t%INFO/MQRankSum" $input -o $result/annotations/MQRankSum
bcftools query -f -c $which "%CHROM\t%ID\t%QUAL\t%INFO/ReadPosRankSum" $input -o $result/annotations/ReadPosRankSum
bcftools query -f -c $which "%CHROM\t%ID\t%QUAL\t%INFO/GQ" $input -o $result/annotations/ReadPosRankSum

echo "Pipeline (vcftools) run finished! yeeeey!"| mail -s "accreditation pipeline" "azzaea@gmail.com"
