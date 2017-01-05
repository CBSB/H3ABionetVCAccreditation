#!/bin/bash
module load bamtstats/1.25
input=~/assessment/results/CBSB_Khartoum/align/CBSB_Khartoum.dedup.bam

java -jar /usr/src/bamstats/bamstats-1.25/BAMStats-1.25.jar -i $input > ~/assessment/results/coverage
