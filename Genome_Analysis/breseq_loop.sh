#!/bin/bash

# This script takes a folder with the short Illumina reads from evolved clones and maps them vs the reference genome of the strain to get genetic variants using breseq

input_dir=$1 # A folder with all the short reads of the evolved clones from a given strain
reference=$2 # A gbk file of the reference strain
output_dir=$3 # The folder in which the output is desired

for reads1 in $input_dir/*_R1_001_val_1.fq.gz
do
	reads2=${reads1%%_R1_001_val_1.fq.gz}"_R2_001_val_2.fq.gz"
	sample=$( basename $( echo ${reads1%%_R1_001_val_1.fq.gz}))
	breseq -r $reference -j 25 -n $sample -o $output_dir/$sample $reads1 $reads2 -p # Write -p flag in case of polyclonal samples, if not, remove
done

