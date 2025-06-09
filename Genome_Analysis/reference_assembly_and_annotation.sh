#!/bin/bash

# This script takes as input the directories in which long and short reads are stored for a given strain, and performs the hybrid assembly of the reference genome

long_data_dir=$1
short_data_dir=$2
out_dir=$3
assembly_dir=$4

for dir in $long_data_dir* # Directory with each folder containing the ID of each strain (folder name) and long reads in a single FASTQ file
do
	sample=$( basename $dir )
	reads1=$short_data_dir/$sample"_R1_001.fastq.gz"
	reads2=$short_data_dir/$sample"_R2_001.fastq.gz"
	long_reads=$long_data_dir/$sample/*.fastq.gz # a single fastq containing the long reads of the given strain
	unicycler -t 17 -1 $reads1 -2 $reads2 -l $long_reads --racon_path /home/pbe/Documents/programs/racon/build/bin/racon -o $out_dir/$sample # Modify racon path to installation path in local machine
done

data_dir=/home/pbe/Downloads/ereA2_interaction/genomes

source activate bakta

for genome in $out_dir/*/*.fasta
do
	sample=$( basename $dir )
	samplename=${sample%%.fasta}
	bakta --db /home/pbe/anaconda3/envs/bakta/databases/db --threads 22 --complete --strain $samplename --output assembly_dir/$samplename $genome # Modify db path to installation path in local machine
done

source deactivate

