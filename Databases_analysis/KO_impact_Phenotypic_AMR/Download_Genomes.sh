#Code for downloading the genomes from the NCBI database using the metadata 
cd ~/IS/Databases/NCBI/Betalactams/
mkdir genomes
cd genomes
email="enteremail@mail.com"

#Set the metadata and creating a new csv with the information of the SAM number 
# and the ID from the database.
csv="/IS/Databases/NCBI/Betalactams/betalactams.csv"
downloads ="/IS/Databases/NCBI/Betalactams/NCBI_betalactams_downloads.csv"

echo "BioSample_ID,Assembly_ID" > "$downloads"
ids=$(awk -F',' 'NR > 1 {print $1}' "$csv")
for biosample_id in $ids; 
do
  echo "Processing BioSample ID: $biosample_id"
  assembly_id=$(esearch -db biosample -query "$biosample_id" -email "$email" | \
                elink -target assembly | \
                esummary | \
                xtract -pattern DocumentSummary -element AssemblyAccession)

  #IF the assembly is not found, continue
  if [ -z "$assembly_id" ]; then
    echo "No assembly found for BioSample $biosample_id"
    continue
  fi
  
  echo "Found Assembly ID: $assembly_id"
  echo "$biosample_id,$assembly_id" >> "$downloads"
  
  # Get the FTP path for the assembly,get the filename from the path, download the FASTA and extract it
 ftp_path=$(esearch -db assembly -query "$assembly_id" -email "$email" | \
             esummary | \
             xtract -pattern DocumentSummary -element FtpPath_GenBank)
             
  #If it could not be found, continue
  if [ -z "$ftp_path" ]; then
    echo "No FTP path found for Assembly $assembly_id"
    continue
  fi
  
  if [ -f ${assembly_id}_genomic.fna ] ; then
 	echo "Already downloaded"
  	continue 
  else 
  	filename=$(basename "$ftp_path")
  	wget -O "${assembly_id}_genomic.fna.gz" "${ftp_path}/${filename}_genomic.fna.gz"
  	gunzip "${assembly_id}_genomic.fna.gz"
  fi 
  
  echo "Downloaded genome for BioSample $biosample_id as ${assembly_id}_genomic.fna"
 
  
done 
