
# Set working directory where your ABRicate summarized tables are stored

setwd("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results")

# Load libraries

library(dplyr)
library(tidyr)
library(xlsx)
library(ggplot2)
library(ggsci)
library(forcats)
library(ggpubfigs)

### Load BVBRC complete+good quality genomes (CSV table filtered and downloaded from the official website)

complete_genomes_metadata <- read.csv("BVBRC_complete_good_genome.csv", colClasses = c("Genome.ID" = "character"))

# Filter by Enterobacteriaceae family, now with those non-enterobacteria

enterobacteria_metadata <- complete_genomes_metadata %>% filter(Family == "Enterobacteriaceae")
non_enterobacteria_metadata <- complete_genomes_metadata %>% filter(Family != "Enterobacteriaceae")

rm(complete_genomes_metadata) # Remove to save mem

# Go to directory with ABRicate results in original format (files obtained from merging by cat command stripping headers from the BVBRC_parsing.sh script)
# And load data

setwd("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results/original_format_non_enterobacteria")

plasmidfinder_tab <- read.csv("plasmidfinder_sumABRicate_non_enterobacteria.tsv", sep = "\t", header = FALSE)
plasmidfinder_tab <- plasmidfinder_tab %>% select(-c(15))
colnames(plasmidfinder_tab) <- c("Genome.ID","SEQUENCE","START","END","STRAND","GENE","COVERAGE","COVERAGE_MAP","GAPS","X.COVERAGE","X.IDENTITY","DATABASE","ACCESSION","PRODUCT")
isfinder_tab <- read.csv("isfinder_sumABRicate_non_enterobacteria.tsv", sep = "\t", header = FALSE)
isfinder_tab <- isfinder_tab %>% select(-c(15))
colnames(isfinder_tab) <- c("Genome.ID","SEQUENCE","START","END","STRAND","GENE","COVERAGE","COVERAGE_MAP","GAPS","X.COVERAGE","X.IDENTITY","DATABASE","ACCESSION","PRODUCT")
megares_tab <- read.csv("megares_w_targets_sumABRicate_non_enterobacteria.tsv", sep = "\t", header = FALSE)
megares_tab <- megares_tab %>% select(-c(15))
colnames(megares_tab) <- c("Genome.ID","SEQUENCE","START","END","STRAND","GENE","COVERAGE","COVERAGE_MAP","GAPS","X.COVERAGE","X.IDENTITY","DATABASE","ACCESSION","PRODUCT")

headers_info <- read.csv("headers_non_enterobacteria_info_table.tsv", sep = "\t", header = TRUE)
colnames(headers_info) <- c("Genome.ID", "SEQUENCE", "Length")


# Reformat Genome.ID column from MEGARes table (strip path and extension)
plasmidfinder_tab$Genome.ID <- gsub("/home/pbe/Documents/database_IS1.2/db_data_non_enterobacteria/genomes/", "", plasmidfinder_tab$Genome.ID)
plasmidfinder_tab$Genome.ID <- gsub(".fna", "", plasmidfinder_tab$Genome.ID)
isfinder_tab$Genome.ID <- gsub("/home/pbe/Documents/database_IS1.2/db_data_non_enterobacteria/genomes/", "", isfinder_tab$Genome.ID)
isfinder_tab$Genome.ID <- gsub(".fna", "", isfinder_tab$Genome.ID)
megares_tab$Genome.ID <- gsub("/home/pbe/Documents/database_IS1.2/db_data_non_enterobacteria/genomes/", "", megares_tab$Genome.ID)
megares_tab$Genome.ID <- gsub(".fna", "", megares_tab$Genome.ID)
headers_info$Genome.ID <- gsub("/home/pbe/Documents/database_IS1.2/db_data_non_enterobacteria/genomes/", "", headers_info$Genome.ID)
headers_info$Genome.ID <- gsub(".fna", "", headers_info$Genome.ID)

# Split gene metadata into useful columns

megares_tab <- megares_tab %>% separate(GENE, c("MEGARES.ID", "DRUG", "DRUG.CLASS", "RESISTANCE.CLASS", "GENE", "QUAL"), sep = "\\|")

# Split IS metadata into useful columns

isfinder_tab <- isfinder_tab %>% separate(GENE, c("IS.element", "IS.family", "IS.group"), sep = "\\_")

# As many replicons contain more than an Inc, handle multiplicity of contigs in the table by merging the Inc genes in the same row (going from 23350 to 15732 contigs identified)
# Keeping the info in all the columns only for the first Inc appearing in the table
# Merging the Inc genes in the "GENE" column split by "/"

plasmidfinder_tab <- plasmidfinder_tab %>%
  group_by(SEQUENCE) %>%
  summarise(
    GENE = paste(sort(GENE, decreasing = TRUE), collapse = "/"), # Order the Inc groups to avoid repetitions
    across(everything(), first)
  ) %>%
  ungroup()

# Add specific ID to both headers_info table and plasmidfinder table which combines both genome and contig to avoid downstream crashing

headers_info$replicon.ID <- paste(headers_info$Genome.ID, headers_info$SEQUENCE, sep = ";")

plasmidfinder_tab$replicon.ID <- paste(plasmidfinder_tab$Genome.ID, plasmidfinder_tab$SEQUENCE, sep = ";")

# Merge plasmidfinder simplified table with header info table to classify replicons in chromosomes and plasmids (by replicon contig ID extracted from the .fna files)

headers_info <- headers_info %>% select(c(3,4)) # Get only length and replicon ID

replicon_tab <- full_join(headers_info, plasmidfinder_tab, by = "replicon.ID")

# Reshape replicon tab to add SEQUENCE information to chromosome or unknown plasmids
# Remove preexisting SEQUENCE column

replicon_tab <- replicon_tab %>% select(-c(3))

replicon_tab <- replicon_tab %>% separate(replicon.ID, c("Genome.ID", "SEQUENCE"), sep = ";")

# Filter out genomes which are not completely closed, i.e. with more than 39 contigs which is the maximum number of plasmids described for a bacterium
# Mitić, N.S., Malkov, S.N., Kovačević, J.J. et al. Structural disorder of plasmid-encoded proteins in Bacteria and Archaea. BMC Bioinformatics 19, 158 (2018). https://doi.org/10.1186/s12859-018-2158-6

n_occurrences <- replicon_tab %>% count(Genome.ID, name = "n_replicons") %>%
  filter(n_replicons < 39)

replicon_tab <- replicon_tab %>% filter(Genome.ID %in% n_occurrences$Genome.ID)

# And include info about replicon (classify into chromosome or plasmid)
# Length of threshold to classify betweem unknown plasmid or chromosome stablished at 350 kbp, based on the observations of
# in diCenzo GC, Finan TM.2017.The Divided Bacterial Genome: Structure, Function, and Evolution. Microbiol Mol Biol Rev81:10.1128/mmbr.00019-17.https://doi.org/10.1128/mmbr.00019-17
# Correct those chromosomes in which plasmidfinder found replicons

replicon_tab <- replicon_tab %>%
  mutate(GENE = ifelse(is.na(GENE) & Length < 350000, "Unclassified_plasmid",
                       ifelse(Length >= 350000, "Chromosome", GENE)))

################ Analyze KOs in AMR-related genes caused by IS elements ################

### Filter disrupted AMR-related genes

# First, we look for potential KO genes by filtering those AMR-related genes present twice in the same contig (i.e. both extremes of the gene disrupted)
# and with a coverage lower than the threshold used to find partial alignments (30%)

target_genes <- megares_tab %>%
  group_by(SEQUENCE, MEGARES.ID) %>%
  filter(X.COVERAGE < 70, n() >= 2) # Consider multiplicity of a gene in the same contig

# Get potential regions in which IS elements might be inserted by getting the extremes of the gene disrupted.
# Setting 60kb as maximum length between start and end disruption points (maximum size described for a transposon
# Bennett, P. M. (1991). Transposable elements and transposition in bacteria. Modern microbial genetics, 324-364.)
# 700 bp as minimum length between start and end disruption points (minimum size for an IS element)
# and coverage maps as indicative of matching extremes... i.e: left end: ==========...... right end: .......=======
# where the potential insertion site would be: ==========........... -> [IS/transposon (> 60kbp)] <- ............=======

# Create empty data frame in which to store potential AMR-related targeted genes

potential_KOs <- data.frame()

# Iterate through rows of AMR-related genes dataframe
for (row in 2:length(rownames(target_genes))) {
  if (target_genes[row-1,]$MEGARES.ID == target_genes[row,]$MEGARES.ID) { # Compare row with previous row, the extremes of the genes are reflected in 2 columns with matching gene names, and filter out genes which confer resistance through SNPs
    breakpoint_begin <- target_genes[row-1,]$END
    breakpoint_end <- target_genes[row,]$START
    breakpoint_range <- breakpoint_end-breakpoint_begin
    #print(breakpoint_range)
    #print(target_genes[row,])
    if (breakpoint_range < 3000 & breakpoint_range > 500) { # Filter by coords matching the criteria
      if (substr(target_genes[row-1,]$COVERAGE_MAP, 15, 15) == "." & substr(target_genes[row,]$COVERAGE_MAP, 1, 1) == "." | substr(target_genes[row-1,]$COVERAGE_MAP, 15, 15) == "=" & substr(target_genes[row,]$COVERAGE_MAP, 1, 1) == "=") { # Filter by coverage map matching both gene ends; account for possible reverse coverage map
        potential_KOs <- rbind(potential_KOs, c(breakpoint_begin, breakpoint_end, target_genes[row,]$Genome.ID, target_genes[row,]$SEQUENCE, target_genes[row,]$DRUG, target_genes[row,]$DRUG.CLASS, target_genes[row,]$RESISTANCE.CLASS, target_genes[row,]$GENE, target_genes[row,]$MEGARES.ID))
      }
    }
  }
}

# Reshape df to correct colnames
colnames(potential_KOs) <- c("Breakpoint_start", "Breakpoint_end", "Genome.ID", "SEQUENCE", "DRUG", "DRUG.CLASS", "RESISTANCE.CLASS", "GENE", "MEGARES.ID")

# Now, having the potential KO genes, the sequence and coordinates in which IS elements might have inserted
# we'll use the potential_KOs dataframe as filter to get those KOs mediated by IS elements (summarized in isfinder_tab).

# Here, we'll merge the info from potential_KOs dataframe and isfinder_tab, adding the IS elements mediating the KOs.

# First, filter out the IS elements of strains which do not contain potential KOs for computational eficiency and save mem

isfinder_tab_KOs <- isfinder_tab %>% filter(Genome.ID %in% potential_KOs$Genome.ID)

# Now, append to potential_KOs IS elements which are in the same sequence and inside the breakpoint range (i.e. disrupting the gene)
# add +- 100 bp of margin to avoid breakpoint errors

KOs_by_ISs <- data.frame(Breakpoint_start = c(NA), Breakpoint_end = c(NA), Genome.ID = c(NA), SEQUENCE = c(NA), DRUG = c(NA), DRUG.CLASS = c(NA), RESISTANCE.CLASS = c(NA), GENE = c(NA), MEGARES.ID = c(NA), IS.element = c(NA), IS.family = c(NA), IS.group = c(NA)) # Init with dummy NA row

#colnames(KOs_by_ISs) <- c("Breakpoint_start", "Breakpoint_end", "Genome.ID", "SEQUENCE", "DRUG", "DRUG.CLASS", "RESISTANCE.CLASS", "GENE", "IS.element", "IS.family", "IS.group")

for (amrrow in 1:length(rownames(potential_KOs))){
  # For each KO gene, define the isertion range
  begin_point <-  as.numeric(potential_KOs[amrrow, ]$Breakpoint_start)-100
  end_point <- as.numeric(potential_KOs[amrrow, ]$Breakpoint_end)+100
  # And append IS elements which match the Genome.ID, Sequence.ID and insertion range
  #print(begin_point)
  #print(end_point)
  ISs_genome <- isfinder_tab_KOs %>% filter(SEQUENCE == potential_KOs[amrrow, ]$SEQUENCE) # Get the IS elements encoded in the target contig
  if (!length(rownames(ISs_genome)) == 0) { # If the contig with the potential KO contains IS elements, begin the loop
    for (isrow in 1:length(rownames(ISs_genome))) { # Parse the IS elements in the contig
      if (ISs_genome[isrow,]$SEQUENCE == potential_KOs[amrrow,]$SEQUENCE & ISs_genome[isrow,]$Genome.ID == potential_KOs[amrrow,]$Genome.ID & ISs_genome[isrow,]$START > begin_point & ISs_genome[isrow,]$END < end_point) { # Ensure the Genome ID and contig coincide, then check that the IS element is between the ends of the disrupted gene
        #print(potential_KOs[amrrow, ]$Genome.ID)
        #print(potential_KOs[amrrow, ]$SEQUENCE)
        #print(potential_KOs[amrrow, ]$GENE)
        #print(ISs_genome[isrow,]$IS.family)
        KOs_by_ISs <- rbind(KOs_by_ISs, c(potential_KOs[amrrow,]$Breakpoint_start, potential_KOs[amrrow,]$Breakpoint_end, potential_KOs[amrrow,]$Genome.ID, potential_KOs[amrrow,]$SEQUENCE, potential_KOs[amrrow,]$DRUG, potential_KOs[amrrow,]$DRUG.CLASS, potential_KOs[amrrow,]$RESISTANCE.CLASS, potential_KOs[amrrow,]$GENE, potential_KOs[amrrow,]$MEGARES.ID ,ISs_genome[isrow,]$IS.element, ISs_genome[isrow,]$IS.family, ISs_genome[isrow,]$IS.group)) # Append individual elements
        #colnames(KOs_by_ISs) <- c("Breakpoint_start", "Breakpoint_end", "Genome.ID", "SEQUENCE", "DRUG", "DRUG.CLASS", "RESISTANCE.CLASS", "GENE", "IS.element", "IS.family", "IS.group")
      }
    }
  }
}

KOs_by_ISs <- na.omit(KOs_by_ISs) # Remove dummy row

KOs_by_ISs %>%
  select(Genome.ID) %>%
  unique() %>%
  count() # From 1609 disrupted genes, 832 are disrupted by IS elements (9.7% of the database)

KOs_by_ISs$TE_length <- as.numeric(KOs_by_ISs$Breakpoint_end) - as.numeric(KOs_by_ISs$Breakpoint_start) # Add length of the TE in between the ends of the gene

####################### Save output #######################
write.xlsx(KOs_by_ISs, "KOs_w_targets_by_ISs_non_enterobacteria_definitive.xlsx", sheetName = "KOs_by_ISs")

# Load to avoid running the loop each execution

KOs_by_ISs <- read.xlsx("KOs_w_targets_by_ISs_non_enterobacteria_definitive.xlsx", sheetIndex = 1)
KOs_by_ISs <- KOs_by_ISs %>% select(-c(1))
#KOs_by_ISs_chr <- read.xlsx("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results/original_format_non_enterobacteria/KOs_by_ISs_w_Tns_non_enterobacteria.xlsx", sheetIndex = 1)

# Rename gene with mismatched fasta header "Chloramphenicol_phosphotransferase" -> "ARAF_like"

KOs_by_ISs$GENE <- gsub("Chloramphenicol_phosphotransferase", "ARAF_like", KOs_by_ISs$GENE)

KOs_by_ISs$GENE.ID <- paste(KOs_by_ISs$GENE, KOs_by_ISs$MEGARES.ID, sep = "_")

# First, classify the KOs caused by a single IS element and those caused by potential transposons (i.e. more than 1 IS in between the ends of the gene)

nKOs <- KOs_by_ISs %>%
  count(Breakpoint_start, Breakpoint_end, Genome.ID, SEQUENCE) # Merge and count by disruption event

nKOs$intgroup <- paste0(nKOs$Breakpoint_start, nKOs$Breakpoint_end, nKOs$Genome.ID, nKOs$SEQUENCE)
nKOs <- nKOs %>% select(n, intgroup)

KOs_by_ISs$intgroup <- paste0(KOs_by_ISs$Breakpoint_start, KOs_by_ISs$Breakpoint_end, KOs_by_ISs$Genome.ID, KOs_by_ISs$SEQUENCE)

KOs_by_ISs <- merge(KOs_by_ISs, nKOs, by = "intgroup")

rm(nKOs) # Remove temporal df

# Classify events

KOs_by_ISs <- KOs_by_ISs %>% mutate(TE_class = case_when(n == 1 ~ "IS", n > 1 ~ "Transposon"))

# Now, analyze those mediated by IS elements (keep out transposons)
# First, merge metadata

KOs_by_ISs <- merge(KOs_by_ISs, non_enterobacteria_metadata, by ="Genome.ID")

# Filter out incomplete assemblies

KOs_by_ISs <- KOs_by_ISs %>% filter(Genome.ID %in% n_occurrences$Genome.ID)

# Remove duplicates in Tn category due to repeated target but different IS element (more than 1 class between begin and end)

KOs_by_ISs <- KOs_by_ISs %>% select(-c(3)) # Remove empty column

# Split, get only Tns

KOs_by_Tns <- KOs_by_ISs %>% filter(TE_class == "Transposon")

# Unify by KO event

KOs_by_Tns <- KOs_by_Tns %>%
  group_by(across(-c(IS.element, IS.family, IS.group))) %>%
  summarise(IS.element = paste(unique(IS.element), collapse = "; "),
            IS.family = paste(unique(IS.family), collapse = "; "),
            IS.group = paste(unique(IS.group), collapse = "; ")) %>%
  ungroup()

KOs_by_Tns$IS.family <- "Tn"

# Now bind to ISs

KOs_by_ISs <- KOs_by_ISs %>% filter(TE_class == "IS")

# And regroup

KOs_by_ISs <- rbind(KOs_by_ISs, KOs_by_Tns)

# Save again

write.xlsx(KOs_by_ISs, "KOs_by_ISs_w_targets_w_Tns_non_enterobacteria_definitive.xlsx")

### Run from here to avoid rerunning the code

KOs_by_ISs <- read.xlsx("KOs_by_ISs_w_targets_w_Tns_non_enterobacteria_definitive.xlsx", sheetIndex = 1)

KOs_by_ISs <- KOs_by_ISs %>% select(-c(1))

KOs_by_ISs <- KOs_by_ISs %>%
  filter(Genus != "", Genus != "Raoultella")

#write.xlsx(KOs_by_ISs, "KOs_by_ISs_w_Tns_non_enterobacteria.xlsx")


# Get those KOs in chromosomal ARGs

replicon_tab <- replicon_tab %>% select(c("SEQUENCE", "GENE"))

colnames(replicon_tab) <- c("SEQUENCE", "Replicon")

KOs_by_ISs <- merge(KOs_by_ISs, replicon_tab, by = "SEQUENCE")

KOs_by_ISs_chr <- KOs_by_ISs %>% filter(Replicon == "Chromosome")

# Save

write.xlsx(KOs_by_ISs_chr, "KOs_by_ISs_w_Tns_chr_non_enterobacteria_definitive.xlsx")
