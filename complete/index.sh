#!/bin/bash

# try indexing:
#1. samtools
#2. sambamba

samtools index ${align_res}/$samplename.dedup.bam
