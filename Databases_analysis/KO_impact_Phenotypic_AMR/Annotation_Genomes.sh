#Code for extracting the ISs and AMR-related genes for each genome. We set the parameters as minimum coverage 5% to better identify all the genes 
# as we have non-closed assemblies. 

genomes="~/IS/Databases/NCBI/Betalactams/genomes"

for i in "$genomes"/*;
do
	strain=$(basename $i | cut -d '.' -f 1,2 | cut -d '_' -f 1,2)
	abricate --db megares_targets --threads 30 --mincov 5 --minid 50 --csv $i >> "~/IS/Databases/NCBI/Betalactams/Megares/$strain.csv"
	abricate --db ISfinder --threads 30 --mincov 5 --minid 50 --csv $i >> "~/IS/Databases/NCBI/Betalactams/ISfinder/$strain.csv"
done


#Unify the annotations of the databases in a single dataframe
for file in ~/IS/Databases/NCBI/Betalactams/ISfinder/*.csv;
do
	filename=$( basename $file)
	tail -n +2 $file >> ~/IS/Databases/NCBI/Betalactams/plamids_Polymyxins_NCBI_data.csv
done

for file in ~/IS/Databases/NCBI/Betalactams/Megares/*.csv;
do
	filename=$( basename $file)
	tail -n +2 $file >> ~/IS/Databases/NCBI/Betalactams/plamids_Polymyxins_NCBI_data.csv
done
