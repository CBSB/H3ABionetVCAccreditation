#!/bin/bash

read1=$1
read2=$2
qc_folder=$3

module load fastqc/0.11.5 

fastqc $read1 $read2 --outdir=$qc_folder
