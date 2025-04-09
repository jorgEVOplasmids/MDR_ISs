setwd("~/IS/Databases/NCBI/Betalactams/")
library(tidyr)
library(dplyr)
library(stringr)
#library(openxlsx)
library(readxl)
library(ggplot2)
library(readr)

### Load our databases: AMR related genes identified with Megares database and Insertion sequences identified with ISfinder
## We also upload the general database of NCBI for betalactams and the sequences that have downloaded correctly.

Meg <- read.csv("Megares_Betalactams_NCBI_data.csv", header = FALSE)
IS <- read.csv("IS_Betalactams_NCBI_data.csv", header = FALSE)
NCBI_original <- read.csv("betalactams.csv", header = TRUE, sep = ",")
NCBI_downloads <- read.csv("NCBI_betalactams_downloads.csv", header = TRUE, sep = ",")

NCBI_downloads_unique <- NCBI_downloads[!duplicated(NCBI_downloads$BioSample_ID), ]

###################
## We filtered the original database with the genomes that have downloaded, merging  both dataframes
# and ruling out those without the Assembly_ID
NCBI <- merge(
  NCBI_downloads_unique,NCBI_original,
  by.x = "BioSample_ID",  
  by.y = "X.BioSample", 
  all.y = TRUE           
)

NCBI <- NCBI[!is.na(NCBI$Assembly_ID), ]

NCBI$File_number <- NA
NCBI$File_number <- sapply(strsplit(NCBI$Assembly_ID, "_"), 
                           function(x) paste(x[1], x[2], sep = "_"))
write.csv(NCBI, "final_lactams_NCBI.csv", row.names =  FALSE)

#------------------------------------------------------------------------------
# We set the same format for all the dataframes
colnames(IS)<-c("File", "Sequence", "Start", "End", "Strand", "Gene", "Coverage",
                "Coverage_map", "Gaps", "%Coverage", "%Identity", "Database", 
                "Accession", "Product", "Resistance")

IS$File <- sub(".*/(.*)\\.fna$", "\\1", IS$File)
IS <- IS[-c(12,13,14,15)]
IS <- IS %>%
  separate(Gene, into=c("Element", "Family", "Gene"), sep = "_")


colnames(Meg)<-c("File", "Sequence", "Start", "End", "Strand", "Gene", "Coverage",
                 "Coverage_map", "Gaps", "%Coverage", "%Identity", "Database", 
                 "Accession", "Product", "Resistance")
Meg$File <- sub(".*/(.*)\\.fna$", "\\1", Meg$File)
Meg <- Meg[-c(12,13,14,15)]
Meg <- Meg %>%
  separate(Gene, into=c("Code", "Type", "Family", "Element", "Gene", "Other"), 
           sep = "\\|")
Meg <- Meg[-c(11)]
 
Meg$File_number <- sapply(strsplit(Meg$File, "_"), 
                          function(x) paste(x[1], x[2], sep = "_"))


# Then we discard all the coincidences in Meg with ISs inserted to avoid false positives.
clean_CTX <- read.csv("abricate_mincov5_ISs_in_megares.csv", header = FALSE)
clean_CTX <- clean_CTX[-c(3,4,5,6,7,8,9,10,11,12,13,14)]
colnames(clean_CTX) <- c("File", "Gene")
clean_CTX <- clean_CTX %>%
  separate(Gene, into=c("Code", "Type", "Family", "Element", "Gene", "Other"), 
           sep = "\\|")

Meg <- Meg[!(Meg$Code %in% clean_CTX$Code), ]

###############################################################################
                          
# In order to detect AMR-related genes disrupted by ISs, we set a loop. 
# If a same genome has data in the Megares and IS dataframe, within the same
# sequence (contig) and the IS is adyacent to the gen, it is added in the
# genral results.
betalactam_general <- data.frame()

for (i in 1:nrow(Meg)) {
  file_gen <- Meg$File[i]
  secuencia_gen <- Meg$Sequence[i]
  start_gen <- Meg$Start[i]
  end_gen <- Meg$End[i]
  
  IS_adjacent <- IS %>%
    filter(File == file_gen &
             Sequence == secuencia_gen & 
             (
               (Start == end_gen) | (Start == end_gen + 1) | 
               (End == start_gen) | (End == start_gen - 1)
             ))
  
  if (nrow(IS_adjacent) > 0) {
    merged_data<- merge(Meg[i, ], IS_adjacent, by = NULL)
    merged_data$Type_of_insertion <- "IS_adjacent"
    betalactam_general <- rbind (betalactam_general, merged_data)
  }
}


betalactam_results <- betalactam_general[-c(5,11,13,16,17,18,25,26,27)]
betalactam_results <- betalactam_results %>% filter(X.Coverage < 100, )
betalactam_results$File_number <- sapply(strsplit(betalactam_results$File, "_"), 
                                         function(x) paste(x[1], x[2], sep = "_"))

# We added the Species to the results and  match the genomes to the different Antibiotics
# and Resistance  phenotype of betalactams. For betalactamics and polymyxin data, it is set
# as a relationship "many-to-many", while for the rest of antibiotics the data is acquired 
# using the same loop used for Species. 
                                         
betalactam_results$Species <- NA
for (i in 1:nrow(betalactam_results)) {
  file_meg <- betalactam_results$File_number[i]
  match_BV <- NCBI %>%
    filter(File_number == file_meg)
  
  if (nrow(match_BV) > 0) {
    betalactam_results$Species[i] <- match_BV$Scientific.name[1] 
  }
}

betalactam_results <- betalactam_results %>% relocate (Species, .before = File.x)
                                         
betalactam_result$Assembly_ID <- sapply(strsplit(betalactam_result$File, "_"), function(x) paste(x[1], x[2], sep = "_"))

betalactam_result <- betalactam_result %>%
  left_join(NCBI %>% select(Assembly_ID, Antibiotic, Resistance.phenotype),
            by = "Assembly_ID", relationship = "many-to-many")
#write_csv(betalactam_result, "betalactam_results.csv", col_names = TRUE)


#Lastly, we identified and count the KO genes and elements to analyse all the atb together
ko_genes <- betalactam_result %>%
  group_by(Antibiotic,Resistance.phenotype, Element)%>%
  count(Gene)

colnames(ko_genes)[5] <- "N"
colnames(ko_genes)[2] <- "Phenotype"
write_csv(ko_genes, "target_betalactams.csv", col_names = TRUE)


ko_elements <- betalactam_results %>%
  group_by(Phenotype)%>%
  count(Element.x)
write_csv(ko_elements, "target_elements.csv", col_names = TRUE)




# Plotting the resistant vs. susceptible samples of general betalactamics. 
resistance_counts <- NCBI %>%
  distinct(BioSample_ID, Resistance.phenotype) %>%  
  count(Resistance.phenotype)                      

resistance_counts <- resistance_counts %>%
  mutate(
    highlight = ifelse(Resistance.phenotype == "sensible" & n == 12, "highlight", "normal")
  )

ggplot(resistance_counts, aes(x = Resistance.phenotype, y = n, fill = Resistance.phenotype)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(
    title = "Count of Resistance Phenotypes",
    x = "Resistance Phenotype",
    y = "Count"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("resistant" = "orange", "sensible" = "grey")) +
  #scale_y_continuous(breaks = seq(0, max(resistance_counts$n) + 5, by = 5)) +
  theme(legend.position = "none")


#separating different ATBs from Betalactamic family 
resistance_betalactams <- NCBI %>%
  distinct(X.BioSample, Antibiotic, Resistance.phenotype) %>%  
  group_by(Antibiotic, Resistance.phenotype) %>%
  summarise(count = n(), .groups = 'drop')

ggplot(resistance_betalactams, aes(x = Resistance.phenotype, y = count, fill = Resistance.phenotype)) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_wrap( ~ Antibiotic, scales = "free_y")
labs(
  title = "Count of Resistance Phenotypes",
  x = "Resistance Phenotype",
  y = "Count"
) +
  theme_classic() +
  scale_fill_manual(values = c("resistant" = "orange", "sensible" = "grey")) +
  #scale_y_continuous(breaks = seq(0, max(resistance_counts$n) + 5, by = 5)) +
  theme(legend.position = "none")
