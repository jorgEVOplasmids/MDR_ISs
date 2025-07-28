
# Set working directory in the folder with the reshaped tables resulting from the reshape_BVBRC_results.sh and reshape_BVBRC_results_non_enterobacteria.sh

setwd("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results")

# Load libraries

library(dplyr)
library(tidyr)
library(xlsx)
library(ggplot2)
library(ggsci)
library(forcats)
library(ggpubr)
library(stringr)
library(tidyverse)

######################################## Load the data ########################################

### Load BVBRC complete+good quality genomes

complete_genomes_metadata <- read.csv("BVBRC_complete_good_genome.csv", colClasses = c("Genome.ID" = "character"))

# Filter by family

enterobacteria_metadata <- complete_genomes_metadata %>% filter(Family == "Enterobacteriaceae")
non_enterobacteria_metadata <- complete_genomes_metadata %>% filter(Family != "Enterobacteriaceae")

rm(complete_genomes_metadata) # Remove to save mem

###### ###### ###### ###### ###### ###### Load non-enterobacteria results ###### ###### ###### ###### ###### 

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

### Reshape metadata to remove repeat Genome.ID

non_enterobacteria_metadata <- non_enterobacteria_metadata %>% filter(Genome.ID %in% n_occurrences$Genome.ID)

# Load Chromosomal KOs

KOs_by_ISs_non_ent <- read.xlsx("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results/original_format_non_enterobacteria/KOs_by_ISs_w_Tns_chr_non_enterobacteria_definitive.xlsx", sheetIndex = 1)

KOs_by_ISs_non_ent <- KOs_by_ISs_non_ent %>% select(-c(1))

KOs_by_ISs_non_ent <- KOs_by_ISs_non_ent %>%
  filter(Genus != "")

# Get those KOs in chromosomal ARGs

replicon_tab <- replicon_tab %>% select(c("SEQUENCE", "GENE"))

colnames(replicon_tab) <- c("SEQUENCE", "Replicon")

KOs_by_ISs_non_ent <- merge(KOs_by_ISs_non_ent, replicon_tab, by = "SEQUENCE")

KOs_by_ISs_chr_non_ent <- KOs_by_ISs_non_ent %>% filter(Replicon.x== "Chromosome")

# Filter out false positives (missannotated CTX genes)

KOs_by_ISs_chr_non_ent <- KOs_by_ISs_chr_non_ent %>% filter(MEGARES.ID != "MEG_2378", MEGARES.ID != "MEG_2430", MEGARES.ID != "MEG_2435", MEGARES.ID != "MEG_8646")

###### ###### ###### ###### ###### ###### Load enterobacteria results ###### ###### ###### ###### ###### 

setwd("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results")

### Load BVBRC complete+good quality genomes

complete_genomes_metadata <- read.csv("BVBRC_complete_good_genome.csv")

# Filter by Enterobacteriaceae family

enterobacteria_metadata <- complete_genomes_metadata %>% filter(Family == "Enterobacteriaceae")

rm(complete_genomes_metadata) # Remove to save mem

# Go to directory with ABRicate results in original format (files obtained from merging by cat command stripping headers)
# And load data

#setwd("D:/database_IS1.2/abricate_summary_results/original_format")
setwd("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/abricate_summary_results/original_format")

plasmidfinder_tab <- read.csv("plasmidfinder_sumABRicate.tsv", sep = "\t")
isfinder_tab <- read.csv("isfinder_sumABRicate.tsv", sep = "\t")
megares_tab <- read.csv("megares_w_targets_sumABRicate.tsv", sep = "\t", header = FALSE)
megares_tab <- megares_tab %>% select(-c(15))
headers_info <- read.csv("headers_info_table.tsv", sep = "\t")

colnames(megares_tab) <- colnames(plasmidfinder_tab)

# Reformat Genome.ID column from MEGARes table (strip path and extension)

megares_tab$Genome.ID <- gsub("/home/pbe/Documents/database_IS1.2/db_data/genomes/", "", megares_tab$Genome.ID)
megares_tab$Genome.ID <- gsub(".fna", "", megares_tab$Genome.ID)

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

# Add metadata to IS and plasmid table
isfinder_tab <- merge(isfinder_tab, enterobacteria_metadata, by = "Genome.ID")
plasmidfinder_tab <- merge(plasmidfinder_tab, enterobacteria_metadata, by = "Genome.ID")

### Reshape metadata to remove repeat Genome.ID

enterobacteria_metadata <- enterobacteria_metadata %>% filter(Genome.ID %in% n_occurrences$Genome.ID)

# Load Chromosomal KOs

KOs_by_ISs_ent <- read.xlsx("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results/original_format/KOs_by_ISs_w_metadata_enterobacteria_definitive.xlsx", sheetIndex = 1)

KOs_by_ISs_ent <- KOs_by_ISs_ent %>% select(-c(1))

KOs_by_ISs_ent <- KOs_by_ISs_ent %>%
  filter(Genus != "", Genus != "Raoultella")

# Get those KOs in chromosomal ARGs

replicon_tab <- replicon_tab %>% select(c("SEQUENCE", "GENE"))

colnames(replicon_tab) <- c("SEQUENCE", "Replicon")

KOs_by_ISs_ent <- merge(KOs_by_ISs_ent, replicon_tab, by = "SEQUENCE")

KOs_by_ISs_chr_ent <- KOs_by_ISs_ent %>% filter(Replicon == "Chromosome")

# Filter out false positives (missannotated CTX genes)

KOs_by_ISs_chr_ent <- KOs_by_ISs_chr_ent %>% filter(MEGARES.ID != "MEG_2378", MEGARES.ID != "MEG_2430", MEGARES.ID != "MEG_2435", MEGARES.ID != "MEG_8646")

# Merge results fixing column names

KOs_by_ISs_chr_non_ent <- KOs_by_ISs_chr_non_ent %>% select(-c(Replicon.y))
colnames(KOs_by_ISs_chr_non_ent) <- colnames(KOs_by_ISs_chr_ent)

KOs_by_ISs_chr <- rbind(KOs_by_ISs_chr_non_ent, KOs_by_ISs_chr_ent)

# Avoid repeated row KOs

KOs_by_ISs_chr <- KOs_by_ISs_chr %>% distinct(Genome.ID, MEGARES.ID, IS.element, .keep_all = TRUE)

### Inspect the results altogether

complete_genomes_metadata <- rbind(non_enterobacteria_metadata, enterobacteria_metadata)

# Classify depending on the KO presence or not

complete_genomes_metadata$KO <- ifelse(complete_genomes_metadata$Genome.ID %in% KOs_by_ISs_chr$Genome.ID, "Y", "N")

# And plot the results together for all the BV-BRC genomes

# First plot the KOs by Family
# And plot fig for genomes in database

complete_genomes_metadata_freq <- complete_genomes_metadata %>%
  filter(Family != "") %>%
  group_by(Family) %>%
  summarise(n = n()) %>%
  filter(n > 300) %>%
  inner_join(complete_genomes_metadata, by = "Family")

top_families <- complete_genomes_metadata_freq$Family

complete_genomes_kod <- complete_genomes_metadata %>%
  filter(Family != "", KO == "Y") %>%
  group_by(Family) %>%
  summarise(n = n())

complete_genomes_metadata_freq <- complete_genomes_metadata_freq %>% mutate(KO = factor(KO, levels = c("Y", "N")))

order_families <- fct_rev(fct_infreq(complete_genomes_metadata_freq$Family)) %>% unique()

# Same plot as S. Fig. 4 A but by family
ggplot(complete_genomes_metadata_freq, aes(x = fct_rev(fct_infreq(Family)), fill = factor(KO, levels = c("N", "Y")))) +
  geom_bar(position = "stack", alpha = 0.8, size = 0.2, col = "black") +
  xlab("Family") +
  ylab("Number of genomes") +
  scale_fill_manual(values = c("Y" = "black", "N" = "grey")) +  # Custom colors
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank())


# Compute relative frequencies
relative_data <- complete_genomes_metadata %>%
  count(Family, KO) %>%  
  group_by(Family) %>%
  mutate(freq = n / sum(n)) %>%  # Calculate proportion
  ungroup() %>%
  mutate(KO = factor(KO, levels = c("Y", "N")))  # Ensure "Y" stacks first

# Plot
# Same as S. Fig. 4 but the frequencies are relative to 
relative_data %>%
  filter(Family %in% top_families) %>%
  ggplot(aes(x = fct_rev(factor(Family, levels = c("Enterobacteriaceae", "Bacillaceae", "Staphylococcaceae",
                                                   "Pseudomonadaceae", "Streptococcaceae", "Lactobacillaceae",
                                                   "Enterococcaceae", "Moraxellaceae", "Mycobacteriaceae",
                                                   "Alcaligenaceae", "Burkholderiaceae", "Xanthomonadaceae",
                                                   "Streptomycetaceae", "Vibrionaceae", "Campylobacteraceae",
                                                   "Pasteurellaceae", "Corynebacteriaceae", "Yersiniaceae",
                                                   "Helicobacteraceae", "Flavobacteriaceae", "Morganellaceae",
                                                   "Rhizobiaceae", "Listeriaceae", "Neisseriaceae", "Erwiniaceae",
                                                   "Clostridiaceae", "Aeromonadaceae", "Brucellaceae", "Paenibacillaceae",
                                                   "Bifidobacteriaceae"))), y = freq, fill = factor(KO, levels = c("N", "Y")))) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, size = 0.2, color = "black") +
  xlab("Family") +
  ylab("Relative Frequency") +
  scale_fill_manual(values = c("Y" = "black", "N" = "grey")) +  # Custom colors
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank())


# Create a table with the frequency

relative_dataframe <- complete_genomes_metadata %>%
  filter(Family != "") %>%
  count(Family, KO) %>%  
  group_by(Family) %>%
  mutate(freq = n / sum(n)) %>%  # Calculate proportion
  ungroup() %>%
  mutate(KO = factor(KO, levels = c("Y", "N")))

relative_dataframe <- relative_dataframe %>% filter(Family %in% top_families, KO == "Y")

###########

# And plot fig for genomes in database

complete_genomes_metadata <- complete_genomes_metadata %>%
  distinct(Genome.ID, .keep_all = TRUE)

complete_genomes_metadata_freq <- complete_genomes_metadata %>%
  filter(Species != "") %>%
  group_by(Species) %>%
  summarise(n = n()) %>%
  filter(n > 200) %>%
  inner_join(complete_genomes_metadata, by = "Species")

top_species <- complete_genomes_metadata_freq$Species

complete_genomes_kod <- complete_genomes_metadata %>%
  filter(Species != "", KO == "Y") %>%
  group_by(Species) %>%
  summarise(n = n())

order_species <- fct_rev(fct_infreq(complete_genomes_metadata_freq$Species)) %>% unique()


complete_genomes_metadata_filtered <- complete_genomes_metadata %>%
  add_count(Species, name = "species_count") %>%
  filter(species_count > 200)

complete_genomes_metadata_filtered_KO <- complete_genomes_metadata_filtered %>% filter(Genome.ID %in% KOs_by_ISs_chr$Genome.ID)

# Plot S. Fig. 4A
ggplot(complete_genomes_metadata_freq, aes(x = fct_rev(fct_infreq(Species)), fill = factor(KO, levels = c("N", "Y")))) +
  geom_bar(aes(fill = Family), position = "stack", alpha = 0.8, size = 0.2, col = "black") +
  xlab("Species") +
  ylab("Number of genomes") +
  scale_fill_manual(values = c("Y" = "black", "N" = "grey")) +  # Custom colors
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank())

# Compute relative frequencies
relative_data <- complete_genomes_metadata %>%
  count(Species, KO) %>%  
  group_by(Species) %>%
  mutate(freq = n / sum(n)) %>%  # Calculate proportion
  ungroup() %>%
  mutate(KO = factor(KO, levels = c("Y", "N")))  # Ensure "Y" stacks first

# Plot
relative_data %>%
  filter(Species %in% top_species) %>%
  ggplot(aes(x = fct_rev(factor(Species, levels = c("Escherichia coli", "Klebsiella pneumoniae", "Salmonella enterica",
                                                    "Staphylococcus aureus", "Pseudomonas aeruginosa", "Acinetobacter baumannii",
                                                    "Mycobacterium tuberculosis", "Bordetella pertussis", "Enterococcus faecalis",
                                                    "Campylobacter jejuni", "Helicobacter pylori", "Bacillus subtilis",
                                                    "Bacillus velezensis", "Enterococcus faecium", "Listeria monocytogenes",
                                                    "Streptococcus pneumoniae", "Streptococcus pyogenes", "Lactiplantibacillus plantarum",
                                                    "Enterobacter hormaechei", "Clostridioides difficile", "Citrobacter freundii",
                                                    "Streptococcus agalactiae"))), y = freq, fill = factor(KO, levels = c("N", "Y")))) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, size = 0.2, color = "black") +
  xlab("Species") +
  ylab("Relative Frequency") +
  scale_fill_manual(values = c("Y" = "black", "N" = "grey")) +  # Custom colors
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank())


# Create a table with the frequency

relative_dataframe <- complete_genomes_metadata %>%
  filter(Species != "") %>%
  count(Species, KO) %>%  
  group_by(Species) %>%
  mutate(freq = n / sum(n)) %>%  # Calculate proportion
  ungroup() %>%
  mutate(KO = factor(KO, levels = c("Y", "N")))

relative_dataframe <- relative_dataframe %>% filter(Species %in% top_species, KO == "Y")

# At species level
complete_genomes_metadata_freq <- complete_genomes_metadata %>%
  filter(Species != "") %>%
  group_by(Species) %>%
  summarise(n = n()) %>%
  filter(n > 200) %>%
  inner_join(complete_genomes_metadata, by = "Species")

top_species <- complete_genomes_metadata_freq$Species

complete_genomes_kod <- complete_genomes_metadata %>%
  filter(Species != "", KO == "Y") %>%
  group_by(Species) %>%
  summarise(n = n())

complete_genomes_metadata_freq <- complete_genomes_metadata_freq %>% mutate(KO = factor(KO, levels = c("Y", "N")))

order_species <- fct_rev(fct_infreq(complete_genomes_metadata_freq$Species)) %>% unique()
# Same as S. Fig. 4 A but relative to the n of each species
ggplot(complete_genomes_metadata_freq, aes(x = fct_rev(fct_infreq(Species)), fill = factor(KO, levels = c("N", "Y")))) +
  geom_bar(position = "stack", alpha = 0.8, size = 0.2, col = "black") +
  xlab("Species") +
  ylab("Number of genomes") +
  scale_fill_manual(values = c("Y" = "black", "N" = "grey")) +  # Custom colors
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank())


# Compute relative frequencies
relative_data <- complete_genomes_metadata %>%
  count(Species, KO) %>%  
  group_by(Species) %>%
  mutate(freq = n / sum(n)) %>%  # Calculate proportion
  ungroup() %>%
  mutate(KO = factor(KO, levels = c("Y", "N")))  # Ensure "Y" stacks first

# Create a table with the frequency

relative_dataframe <- complete_genomes_metadata %>%
  filter(Species != "") %>%
  count(Species, KO) %>%  
  group_by(Species) %>%
  mutate(freq = n / sum(n)) %>%  # Calculate proportion
  ungroup() %>%
  mutate(KO = factor(KO, levels = c("Y", "N")))

relative_dataframe <- relative_dataframe %>% filter(Species %in% top_species, KO == "Y")

######################################################### PLOT THE TARGETS ######################################################### 
# Avoid mismatched rows
KOs_by_ISs_chr$TE_class[KOs_by_ISs_chr$IS.element == "ISPsp7"] <- "IS"
KOs_by_ISs_chr$IS.family[KOs_by_ISs_chr$IS.element == "ISPsp7"] <- "IS30"
KOs_by_ISs_chr$IS.family[KOs_by_ISs_chr$IS.element == "ISAba1"] <- "IS4"

# Define palette for IS families
is_colors <- c(
  "IS1" = "#6a3d9a",
  "IS3" = "#1f78b4",
  "IS4" = "#e6ab02",
  "IS5" = "#33a02c",
  "IS6" = "#fb9a99",
  "IS110" = "#d95f02",
  "IS1380" = "#7570b3",
  "IS1595" = "#e7298a",
  "IS200/IS605" = "#66a61e",
  "IS21" = "#b2df8a",
  "IS256" = "#a6761d",
  "IS30" = "#666666",
  "IS607" = "#fdbf6f",
  "IS630" = "#cab2d6",
  "IS66" = "#ff7f00",
  "IS982" = "#1b9e77",
  "ISL3" = "#b15928",
  "ISNCY" = "#8dd3c7"
)

# S. Fig. 4B
KOs_by_ISs_chr %>%
  filter(Genus != "") %>%
  mutate(
    IS.family = factor(IS.family, levels = names(is_colors))
  ) %>%
  ggplot(aes(x = fct_rev(fct_infreq(IS.family)), fill = IS.family, color = IS.family)) +
  geom_bar(position = "stack", alpha = 0.7) +
  xlab("IS family")+
  ylab("# KOs")+
  #geom_hline(yintercept = 10) +
  #facet_wrap(~Genus, scales = "free_y")+
  #scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Set breaks and labels
  #annotation_logticks(sides = "l")+
  scale_fill_manual(values = is_colors) +
  scale_color_manual(values = is_colors) +
  coord_flip()+
  theme_bw(base_size = 22)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.text.y = element_text(size = 16),
        #strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank())

################################################################################################################

# Define the 4 families to highlight
highlight_families <- c("IS1", "IS3", "IS4", "IS5")

# Define new color palette: keep only the 4, reassign IS3 to IS30's color, and add 'Other' as grey
highlight_colors <- c(
  "IS1" = "#6a3d9a",
  "IS3" = "#1f78b4",  # using IS30's color instead of grey for IS3
  "IS4" = "#e6ab02",
  "IS5" = "#33a02c",
  "Other" = "grey70"
)

# Desired order for stacking and legend
highlight_order <- names(highlight_colors)

# S. Fig. 4C
KOs_by_ISs_chr %>%
  filter(Genus != "") %>%
  add_count(RESISTANCE.CLASS, name = "class_count") %>%
  filter(class_count > 15) %>%
  mutate(
    RESISTANCE.CLASS = str_replace_all(RESISTANCE.CLASS, "_", " "),
    IS.family.cat = ifelse(IS.family %in% highlight_order, IS.family, "Other"),
    IS.family.cat = factor(IS.family.cat, levels = highlight_order)  # normal order for legend
  ) %>%
  ggplot(aes(
    x = fct_rev(fct_infreq(RESISTANCE.CLASS)),
    fill = IS.family.cat,
    color = IS.family.cat
  )) +
  geom_bar(position = position_stack(reverse = TRUE), alpha = 0.7, size = 0.3) +  # Reverse the fill stack only
  xlab("Functional category") +
  ylab("# KOs") +
  scale_fill_manual(values = highlight_colors, name = "IS family") +
  scale_color_manual(values = highlight_colors, name = "IS family") +
  coord_flip() +
  theme_bw(base_size = 22)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank())
# S. Fig. 4D
KOs_by_ISs_chr %>%
  filter(Genus != "") %>%
  add_count(GENE, name = "class_count") %>%
  filter(class_count > 9) %>%
  mutate(
    IS.family.cat = ifelse(IS.family %in% highlight_order, IS.family, "Other"),
    IS.family.cat = factor(IS.family.cat, levels = highlight_order)  # normal order for legend
  ) %>%
  ggplot(aes(
    x = fct_rev(fct_infreq(GENE)),
    fill = IS.family.cat,
    color = IS.family.cat
  )) +
  geom_bar(position = position_stack(reverse = TRUE), alpha = 0.7, size = 0.3) +  # Reverse the fill stack only
  xlab("Target gene") +
  ylab("# KOs") +
  scale_fill_manual(values = highlight_colors, name = "IS family") +
  scale_color_manual(values = highlight_colors, name = "IS family") +
  theme_bw(base_size = 22)+
  coord_flip()+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 22),
        #strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank())

### ### ### ### ### ### ### ### ### Run the loop for all the strains together ### ### ### ### ### ### ### ### ### 

##### Reload replicon tab #####

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

replicon_tab <- replicon_tab %>%
  select(c("Length", "Genome.ID", "SEQUENCE", "GENE"))

# Make grouping variable (Genome.ID + SEQUENCE) and select only replicon identification (GENE)

replicon_tab$Group <- paste(replicon_tab$Genome.ID, replicon_tab$SEQUENCE, sep = "_")

replicon_tab <- replicon_tab %>%
  select(c("Group", "GENE", "Length"))

colnames(replicon_tab) <- c("Group", "Replicon", "Length")

# Now, introduce replicon information into IS table

isfinder_tab$Group <- paste(isfinder_tab$Genome.ID, isfinder_tab$SEQUENCE, sep = "_")

isfinder_tab <- merge(isfinder_tab, replicon_tab, by = "Group")

# And summarize number of IS element per chromosome or plasmid in each genome

n_IS_elements <- isfinder_tab %>%
  group_by(IS.element, Group) %>% ### Summaryze by IS.element to go down 1 level in following analyses
  summarise(n_IS_element = n())

n_IS_elements <- as.data.frame(n_IS_elements)

# Save the info of number of each IS.family per replicon in the replicon table

replicon_tab <- merge(replicon_tab, n_IS_elements)

replicon_tab <- replicon_tab %>% separate(Group, c("Genome.ID", "SEQUENCE"), sep = "_")
replicon_tab$Group <- paste(replicon_tab$Genome.ID, replicon_tab$SEQUENCE, sep = "_")


### Repeat for non-enterobacteria

# Reload all the data from the first lines and get replicon_tab for non enterobacteria and then run the following code


replicon_tab_non_ent <- full_join(headers_info, plasmidfinder_tab, by = "replicon.ID")

# Reshape replicon tab to add SEQUENCE information to chromosome or unknown plasmids
# Remove preexisting SEQUENCE column

replicon_tab_non_ent <- replicon_tab_non_ent %>% select(-c(3))

replicon_tab_non_ent <- replicon_tab_non_ent %>% separate(replicon.ID, c("Genome.ID", "SEQUENCE"), sep = ";")

# Filter out genomes which are not completely closed, i.e. with more than 39 contigs which is the maximum number of plasmids described for a bacterium
# Mitić, N.S., Malkov, S.N., Kovačević, J.J. et al. Structural disorder of plasmid-encoded proteins in Bacteria and Archaea. BMC Bioinformatics 19, 158 (2018). https://doi.org/10.1186/s12859-018-2158-6

n_occurrences <- replicon_tab_non_ent %>% count(Genome.ID, name = "n_replicons") %>%
  filter(n_replicons < 39)

replicon_tab_non_ent <- replicon_tab_non_ent %>% filter(Genome.ID %in% n_occurrences$Genome.ID)

# And include info about replicon (classify into chromosome or plasmid)
# Length of threshold to classify betweem unknown plasmid or chromosome stablished at 350 kbp, based on the observations of
# in diCenzo GC, Finan TM.2017.The Divided Bacterial Genome: Structure, Function, and Evolution. Microbiol Mol Biol Rev81:10.1128/mmbr.00019-17.https://doi.org/10.1128/mmbr.00019-17
# Correct those chromosomes in which plasmidfinder found replicons

replicon_tab_non_ent <- replicon_tab_non_ent %>%
  mutate(GENE = ifelse(is.na(GENE) & Length < 350000, "Unclassified_plasmid",
                       ifelse(Length >= 350000, "Chromosome", GENE)))

replicon_tab_non_ent <- replicon_tab_non_ent %>%
  select(c("Length", "Genome.ID", "SEQUENCE", "GENE"))

# Make grouping variable (Genome.ID + SEQUENCE) and select only replicon identification (GENE)

replicon_tab_non_ent$Group <- paste(replicon_tab_non_ent$Genome.ID, replicon_tab_non_ent$SEQUENCE, sep = "_")

replicon_tab_non_ent <- replicon_tab_non_ent %>%
  select(c("Group", "GENE", "Length"))

colnames(replicon_tab_non_ent) <- c("Group", "Replicon", "Length")

# Now, introduce replicon information into IS table

isfinder_tab$Group <- paste(isfinder_tab$Genome.ID, isfinder_tab$SEQUENCE, sep = "_")

isfinder_tab <- merge(isfinder_tab, replicon_tab_non_ent, by = "Group")

# And summarize number of IS element per chromosome or plasmid in each genome

n_IS_elements <- isfinder_tab %>%
  group_by(IS.element, Group) %>% ### Summaryze by IS.element to go down 1 level in following analyses
  summarise(n_IS_element = n())

n_IS_elements <- as.data.frame(n_IS_elements)

# Save the info of number of each IS.family per replicon in the replicon table

replicon_tab_non_ent <- merge(replicon_tab_non_ent, n_IS_elements)

replicon_tab_non_ent <- replicon_tab_non_ent %>% separate(Group, c("Genome.ID", "SEQUENCE"), sep = "_")
replicon_tab_non_ent$Group <- paste(replicon_tab_non_ent$Genome.ID, replicon_tab_non_ent$SEQUENCE, sep = "_")


replicon_tab_full <- rbind(replicon_tab, replicon_tab_non_ent)


# And now include metadata information to have all samples info for replicon + IS load

metadata_w_replicon_IS_info <- merge(complete_genomes_metadata, replicon_tab_full, by = "Genome.ID")

# First, classify the genomes of the main table in KO by IS or not

metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>%
  mutate(KO = ifelse(Genome.ID %in% KOs_by_ISs_chr$Genome.ID, "yes", "no"))

# Classify replicon in chromosome or plasmid and add to a faceting column (Replicon.class)

metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>%
  mutate(Replicon.class = ifelse(Replicon == "Chromosome", "Chromosome", "Plasmid"))

metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>% select(c("Genome.ID", "Genome.Name", "Genus", "Species", "Plasmids", "Contigs", "Size",
                                                                        "Isolation.Source", "Host.Name", "SEQUENCE", "Replicon", "Length",
                                                                        "IS.element", "n_IS_element", "Group", "KO", "Replicon.class"))

################################################# STATISTICAL TESTS (Results shown in Fig. 3E) ################################################# 

### Contingency table will be:
########################
#  No KO    #    KO    #
################################################################################
# Genomes w/ plasmids w/IS element causing chr KO       #           #          #
################################################################################
# Genomes w/o plasmids w/IS element causing chr KO      #           #          #
################################################################################

# Create a loop to iterate through the elements causing AMR determinants KOs
# And perform statistical test for each IS element

table_pvalues <- data.frame(IS=character(0), pvalue=numeric(0), chi_stat=numeric(0), n_KO_encoding=numeric(0), n_KO_not_encoding=numeric(0), n_no_KO_not_encoding=numeric(0), n_no_KO_encoding=numeric(0))

for(IS_causing_KO in unique(KOs_by_ISs_chr$IS.element)) {
  # Get genomes of Klebsiella showing KOs by the IS element
  genomes_KOd <- KOs_by_ISs_chr %>%
    filter(IS.element == IS_causing_KO)
  
  # Count the number of Genomes
  n_total_KOs <- genomes_KOd %>%
    select(Genome.ID) %>%
    n_distinct()
  
  # Get the replicons and divide the genomes by those which encode the IS element in plasmids and those which do not
  
  n_KOd_encoding <- metadata_w_replicon_IS_info %>%
    filter(Genome.ID %in% genomes_KOd$Genome.ID) %>%
    filter(Replicon.class == "Plasmid", IS.element == IS_causing_KO) %>%
    select(Genome.ID) %>%
    n_distinct()
  
  n_KOd_not_encoding <- n_total_KOs - n_KOd_encoding
  
  # Get the total number of Genomes of Klebsiella which do not show KOs by the IS
  
  n_total_no_KOs <- metadata_w_replicon_IS_info %>%
    filter(!Genome.ID %in% genomes_KOd$Genome.ID) %>%
    select(Genome.ID) %>%
    n_distinct()
  
  # Get the replicons and divide the genomes by those which encode the IS in plasmids and those which do not
  
  n_not_KOd_encoding <- metadata_w_replicon_IS_info %>%
    filter(!Genome.ID %in% genomes_KOd$Genome.ID) %>%
    filter(Replicon.class == "Plasmid", IS.element == IS_causing_KO) %>%
    select(Genome.ID) %>%
    n_distinct()
  
  # And simply substract to get the genomes which do not encode the IS in plasmids
  
  n_no_KOd_not_encoding <- n_total_no_KOs - n_not_KOd_encoding
  
  # Build contingency table
  
  cont_table <- data.frame(
    "NO_KOd" = c(n_no_KOd_not_encoding, n_not_KOd_encoding),
    "KOd" = c(n_KOd_not_encoding, n_KOd_encoding),
    row.names = c("Not_encoding", "Encoding"),
    stringsAsFactors = FALSE
  )
  colnames(cont_table) <- c("No_KOd", "KOd")
  
  test <- fisher.test(cont_table)
  
  table_pvalues <- table_pvalues %>% add_row (IS= IS_causing_KO, pvalue=test$p.value, n_KO_encoding=n_KOd_encoding, n_KO_not_encoding=n_KOd_not_encoding, n_no_KO_not_encoding=n_no_KOd_not_encoding, n_no_KO_encoding=n_not_KOd_encoding)
  
}

# Adjust pvalues

table_pvalues <- table_pvalues %>% filter(!is.na(pvalue))
table_pvalues$padj <- p.adjust(table_pvalues$pvalue, method = "BH")
table_pvalues$category <- ifelse(table_pvalues$padj > 0.05, "ns", "sign")

# Reorder the dataframe by the padj column
table_pvalues <- table_pvalues %>%
  arrange(padj) %>%
  mutate(IS = factor(IS, levels = IS))  # Set factor levels to the ordered IS column
# Define target families
target_families <- c("IS1", "IS3", "IS4", "IS5")

# Replace Family with grouped labels
table_pvalues_repl <- table_pvalues_repl %>%
  mutate(Family_grouped = ifelse(Family %in% target_families, Family, "Other"))

# And plot Fig. 3E

table_pvalues_repl %>%
  ggplot(aes(x = IS, y = -log10(padj), 
             color = Family_grouped, shape = category)) +
  geom_point(alpha = 0.7, size = 3, stroke = 1) +
  xlab("IS") +
  ylab("-log10(p-value)") +
  scale_shape_manual(values = c("ns" = 1, "sign" = 19), name = "Significance") +
  scale_color_manual(values = c(
    "IS1" = "#810f7cff",
    "IS3" = "#3399ccff",
    "IS4" = "#5ca85cff",
    "IS5" = "#f5b841ff",
    "Other" = "#575757ff"
  ), name = "IS Family") +
  theme_bw(base_size = 13) +
  ylim(0,300)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    strip.background = element_blank()
  )

# Filtering ISKpn26 to fix y-axis

table_pvalues_repl %>%
  #filter(IS != "ISKpn26")%>%
  ggplot(aes(x = IS, y = -log10(padj), 
             color = Family_grouped, shape = category)) +
  geom_point(alpha = 0.7, size = 3, stroke = 1) +
  xlab("IS") +
  ylab("-log10(p-value)") +
  scale_shape_manual(values = c("ns" = 1, "sign" = 19), name = "Significance") +
  scale_color_manual(values = c(
    "IS1" = "#810f7cff",
    "IS3" = "#3399ccff",
    "IS4" = "#5ca85cff",
    "IS5" = "#f5b841ff",
    "Other" = "#575757ff"
  ), name = "IS Family") +
  theme_bw(base_size = 13) +
  ylim(0,20)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    strip.background = element_blank()
  )

# Save results to avoid running again the heavy loop
#write.xlsx(table_pvalues, "final_pvalues_plasmidencodedISs.xlsx")

table_pvalues <- read.xlsx("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results/original_format_non_enterobacteria/table_pvalues_complete_version.xlsx", sheetIndex = 1)

########################### NOW ANALYZE THE REPLICON PART (Fig. 3D and S. Fig 5) ###########################

cheatsheet_ISs <- isfinder_tab %>% select(c(IS.element, IS.family, IS.group)) %>% unique()

# Add IS family and group info to the metadata w IS+replicon table

metadata_w_replicon_IS_info_complete <- left_join(metadata_w_replicon_IS_info, cheatsheet_ISs, by = "IS.element")

# Calculate the frequency of a given IS family per species and replicon
highlight_families <- c("IS1", "IS3", "IS4", "IS5")
highlight_species <- c("Klebsiella pneumoniae", "Escherichia coli", "Acinetobacter baumannii")

metadata_w_replicon_IS_info_complete <- metadata_w_replicon_IS_info_complete %>% mutate(IS.main = ifelse(IS.family %in% highlight_families, IS.family, "Other"))
metadata_w_replicon_IS_info_complete <- metadata_w_replicon_IS_info_complete %>% mutate(species.main = ifelse(Species %in% highlight_species, Species, "Other"))

###### Regression #######

### Correlate the number of plasmids within a genome with the number of genomes showing KOs in each category
# show the total n of genomes and the n of genomes showing KOs for each n of plasmids

# Calculate the minimum #Plasmids to keep 95% of the database distribution
filtered_data <- metadata_w_replicon_IS_info %>%
  distinct(Genome.ID, .keep_all = TRUE)
quantile(filtered_data$Plasmids, probs = 0.95) # 5 plasmids

# Fill NA values in the Plasmids column with 0
metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>%
  mutate(
    Plasmids = ifelse(is.na(Plasmids), 0, Plasmids),
    Plasmid_Group = ifelse(Plasmids > 4, ">=5", as.character(Plasmids))
  )

metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>%
  mutate(
    Plasmids = ifelse(is.na(Plasmids), 0, Plasmids))

metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>% filter(Plasmids != ">5")

mean_plasmids <- mean(filtered_data$Plasmids, na.rm = TRUE)

metadata_w_replicon_IS_info %>%
  distinct(Genome.ID, .keep_all = TRUE) %>%
  ggplot(aes(x = Plasmids)) +
  geom_bar(fill = "#A7C6ED", color = "#6A9CCF", size = 1) +
  geom_vline(xintercept = mean_plasmids, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(
    x = "Number of Plasmids",
    y = "Number of Genomes"
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

# Group genomes by n plasmids, keeping those with <=5
metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>%
  mutate(
    Plasmids = ifelse(is.na(Plasmids), 0, Plasmids),
    Plasmid_Group = ifelse(Plasmids > 4, ">=5", as.character(Plasmids))
  )

metadata_w_replicon_IS_info <- metadata_w_replicon_IS_info %>% filter(Plasmids != ">5")

# Group by unique genomes
unique_genomes <- metadata_w_replicon_IS_info %>%
  group_by(Genome.ID) %>%
  summarise(
    Plasmid= first(Plasmids),       # Take the unique plasmid group for each genome
    KO = if_else(any(KO == "yes"), "yes", "no") # Check if any row for the genome has KO as "yes"
  )

# Count genomes by number of plasmids and KO status
count_data <- unique_genomes %>%
  group_by(Plasmid, KO) %>%
  summarise(Count = n(), .groups = "drop")

# Merge "5" and ">=6" into a new group ">=5"
count_data_new <- count_data %>%
  #mutate(Plasmid_Group = ifelse(Plasmid_Group %in% c("5", ">=6"), ">=5", Plasmid_Group)) %>%
  group_by(Plasmid, KO) %>%
  summarise(Count = sum(Count), .groups = "drop")

# Plot the barplot
ggplot(count_data_new, aes(x = factor(Plasmid, levels = c("0", "1", "2", "3", "4", "5")), y = Count, fill = KO, col = KO)) +
  geom_bar(stat = "identity", position = "stack", size = 1) +
  labs(
    x = "Number of Plasmids",
    y = "Number of genomes carrying n plasmids",
    fill = "KO Status"
  ) +
  scale_fill_manual(values = c("#A7C6ED", "#FF8C4A")) +
  scale_color_manual(values = c("#6A9CCF", "#FF7043")) +
  theme_bw(base_size = 20)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank())

# Count genomes by plasmid group and KO status
count_data <- unique_genomes %>%
  group_by(Plasmid, KO) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Plasmid) %>%
  mutate(Proportion = Count / sum(Count))  # Calculate proportions

# Merge "5" and ">=6" into ">=5", and recalculate proportions
count_data_new <- count_data %>%
  #mutate(Plasmid_Group = ifelse(Plasmid_Group %in% c("5", ">=6"), ">=5", Plasmid_Group)) %>%
  group_by(Plasmid, KO) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Plasmid) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Plot the barplot as proportions
ggplot(count_data_new, aes(x = factor(Plasmid, levels = c("0", "1", "2", "3", "4", ">=5")), y = Proportion, fill = KO, col = KO)) +
  geom_bar(stat = "identity", position = "stack", size = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "# Plasmids",
    y = "Genomes with KOs (%)",
  ) +
  scale_fill_manual(values = c("lightgrey", "#FF8C4A")) +
  scale_color_manual(values = c("darkgrey", "#FF7043")) +
  theme_bw(base_size = 20)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        #legend.position = "none",
        aspect.ratio = 1,
        strip.background = element_blank())

# Filter for KO = "yes"
yes_data <- count_data %>%
  filter(KO == "yes") %>%
  mutate(Plasmid_Group = ifelse(Plasmid_Group == ">=6", 6, as.numeric(Plasmid_Group))) # Convert to numeric for plotting

no_data <- count_data %>%
  filter(KO == "no") %>%
  mutate(Plasmid_Group = ifelse(Plasmid_Group == ">=6", 6, as.numeric(Plasmid_Group))) # Convert to numeric for plotting

# Filter for KO = "yes"
yes_data <- count_data_new %>%
  filter(KO == "yes") %>%
  mutate(Plasmid_Group = ifelse(Plasmid_Group == ">=5", 5, as.numeric(Plasmid_Group))) # Convert to numeric for plotting

no_data <- count_data_new %>%
  filter(KO == "no") %>%
  mutate(Plasmid_Group = ifelse(Plasmid_Group == ">=5", 5, as.numeric(Plasmid_Group))) # Convert to numeric for plotting

yes_data <- as.data.frame(yes_data)

# Plot points with correlation stats using stat_cor (Fig. 3D)
yes_data %>%
  filter(Plasmid < 6) %>%
  ggplot(aes(x = Plasmid, y = Proportion * 100)) +
  geom_smooth(method = "lm", se = TRUE, col = "black") +
  stat_cor(method = "pearson", label.x = 0, label.y = max(yes_data$Proportion) * 0.9 * 100, size = 6) + # Note: multiplied by 100 to match y-axis scale
  geom_point(size = 5, color = "#613f75ff", alpha = 0.8) +
  labs(
    x = "# Plasmids",
    y = "Genomes with KOs (%)"
  ) +
  theme_bw(base_size = 22) +
  ylim(-2, 25) +  # Add white space below 0
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    strip.background = element_blank()
  )

######### PLOT PLASMID PRESENCE DEPENDING ON THE KO CATEGORY OF THE GENOMES, RELATIVE TO THE TOTAL DB (%) #########

metadata_w_replicon_IS_info

# Count genomes by number of plasmids and KO status
count_data_replicon_KO <- metadata_w_replicon_IS_info %>%
  filter(Replicon != "Chromosome") %>%
  group_by(Replicon, KO) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(KO) %>%
  mutate(Proportion = Count / sum(Count))  # Calculate proportions

replicon_order <- count_data_replicon_KO %>%
  filter(KO == "yes", Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  arrange(desc(Proportion)) %>%
  pull(Replicon) %>%
  unique()

count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  mutate(Replicon = factor(Replicon, levels = replicon_order)) %>%
  ggplot(aes(y = Replicon, x = Proportion))+
  geom_bar(stat = "identity", position = "stack", size = 1) +
  facet_wrap(~factor(KO, levels = c("yes", "no"))) +
  theme_bw(base_size = 15)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        #legend.position = "none",
        #aspect.ratio = 1,
        strip.background = element_blank())

count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  mutate(Replicon = factor(Replicon, levels = replicon_order)) %>%
  ggplot(aes(y = Replicon, x = Proportion)) +
  geom_bar(stat = "identity", position = "stack", size = 1) +
  facet_wrap(~factor(KO, levels = c("yes", "no"))) +
  theme_bw(base_size = 15) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.background = element_blank())


# Create ordered list based on "yes"
replicon_order <- count_data_replicon_KO %>%
  filter(KO == "yes", Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  arrange(desc(Proportion)) %>%
  pull(Replicon) %>%
  unique()

# Include all replicons (from both yes and no) and reverse
all_replicons <- count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  pull(Replicon) %>%
  unique()

final_order <- rev(c(replicon_order, setdiff(all_replicons, replicon_order)))

# Plot with reversed factor levels
count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  mutate(Replicon = factor(Replicon, levels = final_order)) %>%
  ggplot(aes(y = Replicon, x = Proportion, fill = KO, col = KO)) +
  geom_bar(stat = "identity", position = "stack", size = 1, alpha = 0.7) +
  #facet_wrap(~factor(KO, levels = c("yes", "no"))) +
  theme_bw(base_size = 15) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.background = element_blank())



# Reorder by total (summed) proportion across KO categories
replicon_order <- count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  group_by(Replicon) %>%
  summarise(total_prop = sum(Proportion, na.rm = TRUE)) %>%
  arrange(desc(total_prop)) %>%
  pull(Replicon)

# Apply to plot (S. Fig. 5A)
count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  mutate(Replicon = factor(Replicon, levels = rev(replicon_order))) %>%  # reverse to have biggest at top
  ggplot(aes(y = Replicon, x = Proportion*100, fill = KO, col = KO)) +
  geom_bar(stat = "identity", position = "stack", size = 1, alpha = 0.8) +
  theme_bw(base_size = 18) +
  xlim(0,20)+
  ylab("Plasmid replicons") +
  xlab("Proportion of genomes (%)")+
  scale_fill_manual(values = c("#b3b3b3ff", "#613d73cc"))+
  scale_color_manual(values = c("#b3b3b3ff", "#613d73cc"))+
  facet_wrap(~factor(KO, levels = c("yes", "no")), nrow = 1) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.background = element_blank())

### Perform stats to see if some replicons are enriched in KOd genomes (S. Fig. 5B)

replicon_n_plus10 <- count_data_replicon_KO %>%
  group_by(Replicon) %>%
  summarise(count_n = sum(Count)) %>%
  filter(count_n > 10)

stats_replicon_table <- count_data_replicon_KO %>% filter(Replicon %in% replicon_n_plus10$Replicon)

# Loop over each replicon

### Contingency table will be:
########################
#  No KO    #    KO    #
################################################################################
# Genomes w/ plasmid      #           #          #
################################################################################
# Genomes w/o plasmid     #           #          #
################################################################################

# Create a loop to iterate through the elements causing AMR determinants KOs
# And perform chisq for each IS element

table_pvalues_rep <- data.frame(Replicon=character(0), pvalue=numeric(0), chi_stat=numeric(0))

for(replicon in unique(stats_replicon_table$Replicon)) {
  # Get genomes of Klebsiella showing KOs by the IS element
  genomes_KOd <- metadata_w_replicon_IS_info %>%
    filter(KO == "yes")
  
  # Count the number of Genomes
  n_total_KOs <- genomes_KOd %>%
    select(Genome.ID) %>%
    n_distinct()
  
  # Get the replicons and divide the genomes by those which encode the IS element in plasmids and those which do not
  
  n_KOd_bearing <- metadata_w_replicon_IS_info %>%
    filter(Genome.ID %in% genomes_KOd$Genome.ID) %>%
    filter(Replicon == replicon) %>%
    select(Genome.ID) %>%
    n_distinct()
  
  n_KOd_not_bearing <- n_total_KOs - n_KOd_bearing
  
  n_total_no_KOs <- metadata_w_replicon_IS_info %>%
    filter(!Genome.ID %in% genomes_KOd$Genome.ID) %>%
    select(Genome.ID) %>%
    n_distinct()
  
  # Get the replicons and divide the genomes by those which encode IS26 in plasmids and those which do not
  
  n_not_KOd_bearing <- metadata_w_replicon_IS_info %>%
    filter(!Genome.ID %in% genomes_KOd$Genome.ID) %>%
    filter(Replicon == replicon) %>%
    select(Genome.ID) %>%
    n_distinct()
  
  n_no_KOd_not_bearing <- n_total_no_KOs - n_not_KOd_bearing
  
  # Build contingency table
  
  cont_table_rep <- data.frame(
    "NO_KOd" = c(n_no_KOd_not_bearing, n_not_KOd_bearing),
    "KOd" = c(n_KOd_not_bearing, n_KOd_bearing),
    row.names = c("Not_bearing", "Bearing"),
    stringsAsFactors = FALSE
  )
  colnames(cont_table_rep) <- c("No_KOd", "KOd")
  
  test <- fisher.test(cont_table_rep)
  
  table_pvalues_rep <- table_pvalues_rep %>% add_row (Replicon = replicon, pvalue=test$p.value)
  
}

# Adjust pvalues

table_pvalues_rep <- table_pvalues_rep %>% filter(!is.na(pvalue))

table_pvalues_rep$padj <- p.adjust(table_pvalues_rep$pvalue, method = "BH")

table_pvalues_rep$category <- ifelse(table_pvalues_rep$padj > 0.05, "ns", "sign")

################################################### PLOT info together (Plasmid+IS) ################################################### 

# We can plot the plasmid-IS abundance in a heatmap which could also include the statistics

# We've got the IS info in the ISfinder_tab, and the plasmid replicon info in the replicon_tab
# Hence we can merge both tables, and count the occurrences of each IS element in each replicon

# IMPORTANT!!!! MERGE PREVIOUSLY LOADED REPLICON TABS

# Summarize the n of each IS.element in each sequence ID
# To avoid bias in the sampling of the database, we can get the relative abundance of each IS.element in each replicon and genome

replicon_is_relative <- replicon_tab_full %>%
  filter(Replicon != "Chromosome", Replicon != "Unclassified_plasmid") %>% # filter out chromosome and plasmids not classified
  group_by(Genome.ID, Replicon) %>% # group by plasmid replicon type and is element
  mutate(relative_abundance = n_IS_element/sum(n_IS_element)) %>%
  ungroup()

long_replicon_is <- replicon_tab_full %>%
  filter(Replicon != "Chromosome", Replicon != "Unclassified_plasmid") %>%
  group_by(Genome.ID, Replicon) %>%
  mutate(rel_abundance = n_IS_element / sum(n_IS_element)) %>%
  ungroup() %>%
  select(Replicon, IS.element, rel_abundance) %>%
  group_by(Replicon, IS.element) %>%
  summarise(mean_abundance = mean(rel_abundance), .groups = "drop") %>%
  filter(IS.element %in% KOs_by_ISs_chr$IS.element) %>%
  replace(is.na(.), 0)

# Filter and only get the top replicons

top_replicon <- count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  mutate(Replicon = factor(Replicon, levels = rev(replicon_order)))

top_replicons <- unique(top_replicon$Replicon)

long_replicon_is %>%
  filter(Replicon %in% top_replicons) %>%
  ggplot(aes(x = IS.element, y = factor(Replicon, levels = levels(top_replicons)), fill = mean_abundance)) +
  geom_tile() +
  scale_fill_viridis_c() + # For continuous color scale
  labs(x = "IS Element", y = "Replicon", fill = "Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) # Rotate x-axis labels for better readability

### Plot in order of IS elements, based on their family and order by prevalence in the knock-outs

long_replicon_is %>%
  filter(Replicon %in% top_replicons) %>%
  ggplot(aes(x = factor(IS.element, levels = c(table_pvalues$IS)), y = factor(Replicon, levels = levels(top_replicons)), fill = mean_abundance)) +
  geom_tile() +
  scale_fill_viridis_c() + # For continuous color scale
  labs(x = "IS Element", y = "Replicon", fill = "Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) # Rotate x-axis labels for better readability


# We can do the same but normalizing by genome, as the replicon content varies a lot


replicon_is_relative <- replicon_tab_full %>%
  filter(Replicon != "Chromosome", Replicon != "Unclassified_plasmid") %>% # filter out chromosome and plasmids not classified
  group_by(Genome.ID) %>% # group by plasmid replicon type and is element
  mutate(norm_abundance = n_IS_element/sum(n_IS_element)) %>%
  ungroup()

long_replicon_is <- replicon_tab_full %>%
  filter(Replicon != "Chromosome", Replicon != "Unclassified_plasmid") %>% # filter out chromosome and plasmids not classified
  group_by(Genome.ID) %>% # group by plasmid replicon type and is element
  mutate(norm_abundance = n_IS_element/sum(n_IS_element)) %>%
  ungroup() %>%
  select(Replicon, IS.element, norm_abundance) %>%
  group_by(Replicon, IS.element) %>%
  summarise(mean_abundance = mean(norm_abundance), .groups = "drop") %>%
  filter(IS.element %in% KOs_by_ISs_chr$IS.element) %>%
  replace(is.na(.), 0)

# Filter and only get the top replicons

top_replicon <- count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  mutate(Replicon = factor(Replicon, levels = rev(replicon_order)))

count_data_replicon_KO %>%
  filter(Proportion > 0.01, Replicon != "Unclassified_plasmid") %>%
  mutate(Replicon = factor(Replicon, levels = rev(replicon_order))) %>% filter(KO == "yes")

top_replicons <- unique(count_data_replicon_KO$Replicon)

table_pvalues <- read.xlsx("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/database_IS1.2/last_version/database_IS1.2/abricate_summary_results/original_format_non_enterobacteria/table_pvalues_complete_version.xlsx", sheetIndex = 1)

# Add IS family to long_replicon_is

cheatsheet_ISs <- read.csv("/home/jorge/Downloads/sequences_IS_headers.csv", sep = "\t")

colnames(cheatsheet_ISs) <- c("IS.element", "IS.family", "IS.group")
cheatsheet_ISs <- cheatsheet_ISs %>% distinct(IS.element, IS.family, IS.group)

cheatsheet_ISs_sorted <- cheatsheet_ISs[order(cheatsheet_ISs$IS.family),]
cheatsheet_ISs_sorted$IS.element <- as.factor(cheatsheet_ISs_sorted$IS.element)

cheatsheet_ISs_sorted <- cheatsheet_ISs_sorted[!duplicated(cheatsheet_ISs_sorted$IS.element),]

top_replicons <- count_data_replicon_KO %>%
  filter(Proportion > 0.01) %>% select(Replicon)

# Heatmap shown in S. Fig. 5B
long_replicon_is %>%
  filter(Replicon %in% top_replicons) %>%
  ggplot(aes(x = factor(IS.element, levels = c(cheatsheet_ISs_sorted$IS.element)), y = factor(Replicon, levels = levels(top_replicons)), fill = as.numeric(mean_abundance))) +
  geom_tile() +
  #scale_fill_distiller(palette = "Purples", direction = 1) +
  scale_fill_gradient(low = "grey97", high = "steelblue")+
  labs(x = "IS Element", y = "Replicon", fill = "Mean Relative Abundance") +
  theme_bw(base_size = 15) +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank())

long_replicon_is_filtered <- long_replicon_is %>%
  filter(Replicon %in% top_replicons)

#### CHECK WHETHER A PLASMID-IS COMBINATION IS ENRICHED IN THE KO SUBGROUP (Results shown in S. Fig. 5B) ####

replicon_tab_enrichment <- replicon_tab_full %>% mutate(KO = ifelse(Genome.ID %in% KOs_by_ISs_chr$Genome.ID, "yes", "no"))

# Filter out chromosomal and non classified plasmids

replicon_tab_enrichment_filtered <- replicon_tab_enrichment %>%
  filter(Replicon != "Chromosome", Replicon != "Unclassified_plasmid") %>% # filter out chromosome and plasmids not classified
  distinct(Genome.ID, Replicon, IS.element, KO) %>%
  mutate(presence = 1)

# Get 0's

enrichment_combinations_full <- expand.grid(Genome.ID = unique(replicon_tab_enrichment_filtered$Genome.ID),
                                            Replicon = unique(replicon_tab_enrichment_filtered$Replicon),
                                            IS.element = unique(replicon_tab_enrichment_filtered$IS.element), stringsAsFactors = FALSE)


# More efficient solution

replicon_tab_enrichment <- replicon_tab_full %>% mutate(KO = ifelse(Genome.ID %in% KOs_by_ISs_chr$Genome.ID, "yes", "no"))

filtered_tab <- replicon_tab_enrichment %>%
  filter(Replicon != "Chromosome", Replicon != "Unclassified_plasmid") %>%
  distinct(Genome.ID, Replicon, IS.element, KO)


combos <- filtered_tab %>%
  distinct(Replicon, IS.element)

genomes_by_group <- filtered_tab %>%
  distinct(Genome.ID, KO) %>%
  count(KO) %>%
  deframe()

# Main loop with safety checks
fisher_results <- map_dfr(1:nrow(combos), function(i) {
  this_combo <- combos[i, ]
  
  # Genomes where this combo is present
  presence <- filtered_tab %>%
    filter(Replicon == this_combo$Replicon,
           IS.element == this_combo$IS.element) %>%
    distinct(Genome.ID, KO)
  
  # Count of presence per group
  presence_counts <- presence %>%
    count(KO) %>%
    complete(KO = c("yes", "no"), fill = list(n = 0)) %>%
    deframe()
  
  # Calculate presence/absence
  a <- presence_counts["yes"] %||% 0  # Present in KO = yes
  c <- presence_counts["no"] %||% 0   # Present in KO = no
  b <- (genomes_by_group["yes"] %||% 0) - a  # Absent in KO = yes
  d <- (genomes_by_group["no"] %||% 0) - c   # Absent in KO = no
  
  # Only run test if all values are nonnegative and finite
  if (any(c(a, b, c, d) < 0) || any(!is.finite(c(a, b, c, d)))) {
    return(NULL)
  }
  
  test <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  
  tibble(
    Replicon = this_combo$Replicon,
    IS.element = this_combo$IS.element,
    KO_yes_present = a,
    KO_yes_absent = b,
    KO_no_present = c,
    KO_no_absent = d,
    p_value = test$p.value,
    odds_ratio = unname(test$estimate)
  )
})

# Adjust p-values
fisher_results <- fisher_results %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_adj)

#write.xlsx(fisher_results, "plasmid_is_association_results.xlsx")

fisher_results_filtered <- fisher_results %>% filter(Replicon %in% long_replicon_is$Replicon, IS.element %in% long_replicon_is$IS.element)

#write.xlsx(fisher_results_filtered, "fisher_results_filtered.xlsx")
