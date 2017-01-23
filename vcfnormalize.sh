#!/bin/bash

module load bcftools/1.3.1 

bcftools isec   $path_to_gatkvcf $path_to_freebayesvcf  -p -c both
