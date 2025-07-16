#!/bin/bash

# Loop to download BV-BRC genomes, from the list of ID downloaded from the filtered (complete and good quality) metadata table from (https://www.bv-brc.org/view/Bacteria/2)
# stored in a local file called IDs_complete_genomes_bacteria.csv

for i in `cat /home/pbe/Documents/database_IS1.2/db_data/IDs_complete_genomes_bacteria.csv`
do
	wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.fna"
done

# Point towards downloaded genomes folder

data_dir=$1
out_dir=$2

# Get each contig length

echo -e "File\tID\tLength" > "$output_file"
for fasta_file in $data_dir/*.fna; do
    awk -v file="$fasta_file" '
    /^>/ {
        if (seqlen > 0) {
            print file "\t" id "\t" seqlen
        }
        id = substr($1, 2)
        seqlen = 0  # Reset sequence length for new sequence
    }
    !/^>/ {
        seqlen += length($0)
    }
    END {
        if (seqlen > 0) {
            print file "\t" id "\t" seqlen
        }
    }' "$fasta_file" >> "$output_file"
done

# This next code takes the FASTA files of the BV-BRC genomes downloaded (in data_dir) and annotates the presence of partial ARGs, full IS elements and plasmid replicons with plasmidfinder
# Then summarizes the abricate files into individual tables (per database)

# Loop over folder and run ABRicate

for genome in /home/pbe/Documents/database_IS1.2/db_data_non_enterobacteria/genomes/*.fna
do
	barc=$( basename $genome )
	name=${barc%%.fna}
	echo $barc $name
	/home/pbe/Documents/programs/abricate/bin/abricate --minid 80 --mincov 30 --threads 28 --db megares_w_targets $genome > $out_dir/megares/$name.tab
	/home/pbe/Documents/programs/abricate/bin/abricate --threads 28 --db ISfinder_db $genome > $out_dir/ISfinder_db/$name.tab
	/home/pbe/Documents/programs/abricate/bin/abricate --threads 28 --db plasmidfinder $genome > $out_dir/plasmidfinder/$name.tab
done

# cat abricate files and summarize into a table (with the same columns as the ABRicate file, not the ABRicate --summary option)

# point towards your output directory (depending on the database used in ABRicate)

abricate_out=$3

for file in $abricate_dir/*.tab
do
	filename=$( basename $file )
	gID=${filename%%.tab}
	tail -n +2 $file >> $abricate_out/db_sumABRicate_bacteria.tab
done
