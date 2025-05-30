setwd("~/IS/Databases/NCBI/")
#library(xlsx)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)
library(readr)
library(ggplot2)
library(wesanderson)
library(scales)
library(patchwork)
library(ggh4x)
library(ggbreak)


#Upload the databases for each ATB and combine them 
chloram_NCBI <- read.csv("Chloram/final_chloram_NCBI.csv")
cipro_NCBI <- read.csv("Cipro/final_cipro_NCBI.csv")
polymyxin_NCBI <- read.csv("Polymyxins/final_polymyxins_NCBI.csv")
kana_NCBI <- read.csv("Kana/final_kana_NCBI.csv")
lactams_NCBI <- read.csv("Betalactams/final_lactams_NCBI.csv")


#NCBI_DATABASE <- rbind(chloram_NCBI, kana_NCBI, cipro_NCBI, chloram_NCBI,
                       #lactams_NCBI, polymyxin_NCBI)
#write.csv(NCBI_DATABASE, "NCBI_complete_database.csv", row.names = FALSE)



#Function to calculate the number o resistant/susceptibles per antibiotic
#(in the case of betalactams and polymyxins broken down)

calculate_resistance_distribution <- function(df) {
  resistance_counts <- df %>%
    group_by(Antibiotic, Resistance.phenotype) %>%
    summarise(Count = n(), .groups = 'drop')
  
  antibiotic_sum <- resistance_counts %>%
    pivot_wider(
      names_from = Resistance.phenotype, 
      values_from = Count, 
      values_fill = 0
    ) %>%
    mutate(Total = resistant + susceptible)
  
  totals <- data.frame(
    Total_Resistant = sum(antibiotic_sum$resistant),
    Total_Susceptible = sum(antibiotic_sum$susceptible),
    Total_Strains = sum(antibiotic_sum$Total)
  )
  
  return(list(
    antibiotic_sum = antibiotic_sum,
    totals = totals
  ))
}

chloram_distribution<- calculate_resistance_distribution(chloram_NCBI)
cipro_distribution <- calculate_resistance_distribution(cipro_NCBI)
kana_distribution <- calculate_resistance_distribution(kana_NCBI)
polymyxin_distribution <- calculate_resistance_distribution(polymyxin_NCBI)
betalactams_distribution <- calculate_resistance_distribution(lactams_NCBI)


#general_NCBI <- rbind(chloram_general_distribution, cipro_general_distribution, 
#                      kana_general_distribution, betalactams_general_distribution, 
#                      polymyxin_general_distribution)


general_distribution_NCBI <-  rbind(chloram_distribution$antibiotic_sum, 
                                    cipro_distribution$antibiotic_sum, 
                                    kana_distribution$antibiotic_sum, 
                                    betalactams_distribution$antibiotic_sum, 
                                    polymyxin_distribution$antibiotic_sum)

#write.csv(general_distribution_NCBI, "distribution_NCBI_database.csv", 
           #row.names = FALSE)

rm(chloram_distribution, cipro_distribution, kana_distribution, 
   betalactams_distribution, polymyxin_distribution)


#Set relevant targets for the antibiotics, except colistin and polymyxin B
elements <- c("Drug_and_biocide_ABC_efflux_pumps",
              "Drug_and_biocide_and_metal_RND_efflux_pumps",
              "Drug_and_biocide_and_metal_RND_efflux_regulator",
              "Drug_and_biocide_MATE_efflux_pumps",
              "Drug_and_biocide_MFS_efflux_pumps",                     
              "Drug_and_biocide_MFS_efflux_regulator",
              "Drug_and_biocide_RND_efflux_pumps",                     
              "Drug_and_biocide_RND_efflux_regulator",
              "Drug_and_biocide_SMR_efflux_pumps",                     
              "Drug_and_metal_MFS_efflux_pumps",
              "MDR_mutant_porin_proteins",
              "MLS_resistance_MFS_efflux_pumps",
              "Multi-drug_RND_efflux_pumps",                           
              "Multi-drug_RND_efflux_regulator",
              "Mutant_porin_proteins",
              "Phenicol_resistance_MFS_efflux_pumps",
              "Tetracycline_resistance_MFS_efflux_pumps",
              "Tetracycline_resistance_MFS_efflux_regulator",          
              "Tetracycline_resistance_ribosomal_protection_proteins",
              "Tetracycline_transcriptional_repressor")


#Now, we upload the results from our previous analysis (both the raw results and 
# the elements counts) and filter them using the relevant targets

targets_chloram <- read.csv("Chloram/target_chloram.csv")
targets_chloram$atb <- "chloramphenicol"
targets_chloram <- targets_chloram[targets_chloram$Element %in% elements, ]



targets_cipro <- read_xlsx("Cipro/target_cipro.xlsx", sheet = 2)
targets_cipro$atb <- "ciprofloxacin"
targets_cipro <- targets_cipro[targets_cipro$Element %in% elements,]

  
targets_colistin <- read.csv("Polymyxins/target_poly.csv")
targets_colistin <- targets_colistin[targets_colistin$Element == "Colistin-resistant_mutant", ]
  
targets_kana <- read.csv("Kana/targets_kanamycin.csv")
targets_kana <- targets_kana[targets_kana$Element %in% elements, ]
targets_kana$atb<- "kanamycin"


targets_betalac <- read.csv("Betalactams/target_betalactams.csv")
targets_betalac <- targets_betalac[targets_betalac$Element %in% elements, ]


combined_atb <- rbind(targets_cipro, targets_chloram, targets_kana, 
                      targets_colistin, targets_betalac)

rm(targets_cipro, targets_chloram, targets_colistin, targets_kana, targets_betalac)
#write.csv(combined_atb, "all_relevant_targets.csv")


################################################################################

#Now, we calculate the frequency of having one of these genes disrupted by an IS 
# considering the data samples (resistant vs susceptible)
#We discarded those atbs with less than 100 samples overall and those with less
#than 5 samples in one of the phenotypic condiitons.

colnames(combined_atb)[5] <- "Antibiotic"

total_frequencies <- combined_atb %>%   
  mutate(Phenotype = tolower(Phenotype)) %>%   
  left_join(
    general_distribution_NCBI %>%
      select(Antibiotic, susceptible, resistant),
    by = "Antibiotic"
  ) %>%   
  filter(susceptible + resistant >= 100) %>%
  filter (susceptible > 5 & resistant > 5) %>%
  group_by(Antibiotic, Phenotype, Element) %>%   
  mutate(
    frequency = case_when(
      Phenotype == "susceptible" ~ N / susceptible,
      Phenotype == "resistant" ~ N / resistant
    )
  ) %>%   
  ungroup() %>%   
  select(Antibiotic, Phenotype, Element, Gene, N, frequency)

#For better analysing the antibiotics, we add the antibiotic family and group them
total_frequencies$Family <- NA
total_frequencies <- total_frequencies %>%
  mutate(Family = case_when(
    Antibiotic == "ciprofloxacin" ~ "Fluoroquinolones",
    Antibiotic == "chloramphenicol" ~ "Chloramphenicol",
    Antibiotic == "kanamycin" ~ "Aminoglycosides",
    Antibiotic == "colistin" ~ "Polymyxins",
    Antibiotic == "polymyxin B" ~ "Polymyxins",
    Antibiotic == "amoxicillin-clavulanic acid" ~ "Penicillin",
    Antibiotic == "ampicillin" ~ "Penicillin",
    Antibiotic == "ampicillin-sulbactam" ~ "Penicillin",
    Antibiotic == "aztreonam" ~ "Monobactam",
    Antibiotic == "cefazolin" ~ "Cephalosphorin 1 generation",
    Antibiotic == "cefepime" ~ "Cephalosphorin 4 generation",
    Antibiotic == "cefotaxime" ~ "Cephalosphorin 3 generation",
    Antibiotic == "cefotaxime-clavulanic acid" ~ "Cephalosphorin 3 generation",
    Antibiotic == "cefotetan" ~ "Cephalosphorin 2 generation",
    Antibiotic == "cefoxitin" ~ "Cephalosphorin 2 generation",
    Antibiotic == "cefpodoxime" ~ "Cephalosphorin 3 generation",
    Antibiotic == "ceftazidime" ~ "Cephalosphorin 3 generation",
    Antibiotic == "ceftazidime-avibactam" ~ "Cephalosphorin 3 generation",
    Antibiotic == "ceftazidime-clavulanic acid" ~ "Cephalosphorin 3 generation",
    Antibiotic == "ceftiofur" ~ "Cephalosphorin 3 generation",
    Antibiotic == "ceftolozane-tazobactam" ~ "Cephalosphorin 5 generation",
    Antibiotic == "ceftriaxone" ~ "Cephalosphorin 3 generation",
    Antibiotic == "cefuroxime" ~ "Cephalosphorin 2 generation",
    Antibiotic == "cephalothin" ~ "Cephalosphorin 1 generation",
    Antibiotic == "doripenem" ~ "Carbapenem",
    Antibiotic == "ertapenem" ~ "Carbapenem",
    Antibiotic == "imipenem" ~ "Carbapenem",
    Antibiotic == "meropenem" ~ "Carbapenem",
    Antibiotic == "piperacillin-tazobactam" ~ "Penicillin",
    Antibiotic == "ticarcillin-clavulanic acid" ~ "Penicillin"
  ))


# We format the genes to be italic and plot the frequency per gene and condition
format_gene_name <- function(gene) {
  first_three <- substr(gene, 1, 3)
  rest <- substr(gene, 4, nchar(gene))
  return(paste0(tolower(first_three), rest))
}

gene_labels <- setNames(
  paste0("italic('", sapply(unique(total_frequencies$Gene), format_gene_name), "')"),
  unique(total_frequencies$Gene)
)

ggplot(total_frequencies, aes(x = Phenotype, y = frequency, fill = Gene)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Antibiotic, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(labels = parse(text = gene_labels)) + 
  labs(
    title = "Gene Frequency by Phenotype",
    x = "Phenotype",
    y = "Frequency",
    fill = "Gene"
  ) +
  theme_bw(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(), 
  )

#Sum up the frequencies so we dont see a break and set an order
total_frequencies2 <- total_frequencies %>% 
  group_by(Family, Phenotype, Element) %>%
  summarise(frequency = sum(frequency), .groups = 'drop')

total_frequencies2$General_Family <- NA
total_frequencies2 <- total_frequencies2 %>%
  mutate(General_Family = case_when(
    Family == "Cephalosphorin 1 generation" ~ "Cephalosphorin",
    Family == "Cephalosphorin 2 generation" ~ "Cephalosphorin",
    Family == "Cephalosphorin 3 generation" ~ "Cephalosphorin",
    Family == "Cephalosphorin 4 generation" ~ "Cephalosphorin",
    Family == "Cephalosphorin 5 generation" ~ "Cephalosphorin",
    Family == "Aminoglycosides" ~ "Aminoglycosides",
    Family == "Carbapenem" ~ "Carbapenem", 
    Family == "Monobactam" ~ "Monobactam", 
    Family == "Fluoroquinolones" ~ "Fluoroquinolones", 
    Family == "Chloramphenicol" ~ "Chloramphenicol", 
    Family == "Penicillin" ~ "Penicillin", 
    Family == "Polymyxins" ~ "Polymyxins"
  ))


# We plot the frequency per element and condition firstly with the genes in the
#y axis and all atb together (grouping by phenotype) and secondly with the ATB 
#separated and filling the elements

ggplot(total_frequencies2, aes(x = frequency, y = Element, color = General_Family, 
                               fill = General_Family)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
  facet_grid(~ Phenotype, scales = "free_y") +
  scale_x_continuous(labels = scales::percent) +
  labs(
    title = "Broken Elements per Antibiotic Family",
    x = "Frequency",
    y = "Element",
    fill = "General_Family"
  ) +
  theme_bw(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(), 
  )


ggplot(total_frequencies, aes(x =  Phenotype, y = frequency, fill = Element)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ Antibiotic, scales = "free_y") +
  labs(
    title = "Element Frequency by Phenotype",
    x = "Phenotype",
    y = "Frequency",
    fill = "Element"
  ) +
  theme_bw(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.position = "bottom",
    legend.key.size = unit(0.1, "lines"),  
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)  
    +
      guides(fill = guide_legend(ncol = 1)) 
  )


#

## Reorder by antibiotic family
desired_order <- c("Aminoglycosides", "Chloramphenicol", "Fluoroquinolones", 
                   "Polymyxins", "Cephalosphorin 1 generation",
                   "Cephalosphorin 2 generation", "Cephalosphorin 3 generation",
                   "Cephalosphorin 4 generation","Cephalosphorin 5 generation",
                   "Carbapenem", "Monobactam", "Penicillin" 
                   )

total_frequencies3 <- total_frequencies %>% 
  group_by(Family, Antibiotic, Phenotype) %>%
  summarise(frequency = sum(frequency), .groups = 'drop')


total_frequencies3$Family <- factor(total_frequencies3$Family, 
                                    levels = desired_order)

ggplot(total_frequencies3, aes(x = Phenotype, y = frequency, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0.1, 0.))  
  ) +
  facet_nested_wrap(Family ~ Antibiotic, scales = "free_y", ncol = 5) +
  labs(
    title = "Relevant knock-out genes frequency per phenotype",
    x = "Phenotype",
    y = "Frequency (%)",
    fill = "Antibiotic"
  ) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# Now we plot only our experimental atb and the 10 betalactams with the higher 
#number of samples (covering all the sublcases of the betalactams: Cephalosporin,
# Penicillins, Carbapenems and Monobactams).

atb <- c( "chloramphenicol", "ciprofloxacin", "chloramphenicol", "kanamycin",
          "colistin", "polymyxin B", "amoxicillin-clavulanic acid", 
          "ceftriaxone","ampicillin", "meropenem", "ceftiofur", "ceftazidime", 
          "piperacillin-tazobactam","cefepime", "ertapenem", "imipenem", 
          "cefazolin", "cefoxitin", "aztreonam","cefotaxime", 
          "ampicillin-sulbactam")


total_frequencies2 <- total_frequencies[total_frequencies$Antibiotic %in% atb, ]
desired_order2 <- c("Aminoglycosides", "Chloramphenicol", "Fluoroquinolones", 
                   "Polymyxins", "Cephalosphorin 1 generation",
                   "Cephalosphorin 2 generation", "Cephalosphorin 3 generation",
                   "Cephalosphorin 4 generation","Penicillin","Carbapenem", 
                    "Monobactam")

total_frequencies2$Family <- factor(total_frequencies2$Family, 
                                   levels = desired_order2)

ggplot(total_frequencies2,  aes(x =  Phenotype, y = frequency, fill = Phenotype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap( Family ~ Antibiotic, scales = "free_y", ncol = 5) +
  labs(
    title = "Relevant knock-out genes frequency per phenotype",
    x = "Phenotype",
    y = "Frequency (%)",
    fill = "Antibiotic"
  ) +
  theme_bw(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


### STATISTICAL ANALYSIS

# Function to extract data for the contingency table, handling NAs as 0 and 
#performing fisher test for each atb with more than 100 genomes.

perform_fisher_test <- function(antibiotic_name, general_df, combined_df) {
  general_data <- general_df %>%
    filter(Antibiotic == antibiotic_name, Total > 100) %>%
    select(resistant, susceptible)
  
  # Print the raw data for each ATB and check possible NAs
  cat("\nRaw data for", antibiotic_name, ":\n")
  print(general_data)
  
  if(any(is.na(general_data$resistant)) || any(is.na(general_data$susceptible))) {
    cat("  Skipping", antibiotic_name, "- contains NA values\n")
    return(NULL)
  }
  
  # Get total resistant/susceptible counts and extract counts of broken genes
  total_resistant <- general_data$resistant[1]
  total_susceptible <- general_data$susceptible[1]
  total <- total_resistant + total_susceptible
  
  broken_counts <- combined_df %>%
    filter(Antibiotic == antibiotic_name) %>%
    group_by(Phenotype) %>%
    summarize(broken_count = sum(N, na.rm = TRUE))
  
  # Check if we have data for both phenotypes in the combined df and set 0
  # if one phenotype is missing 
  has_resistant <- "resistant" %in% broken_counts$Phenotype
  has_susceptible <- "susceptible" %in% broken_counts$Phenotype
  
  broken_resistant <- 0
  broken_susceptible <- 0
  
  # Update with actual values if htere are
  if(has_resistant) { broken_resistant <- broken_counts$broken_count[broken_counts$Phenotype
                                                                     == "resistant"]}
  if(has_susceptible) {broken_susceptible <- broken_counts$broken_count[broken_counts$Phenotype 
                                                                        == "susceptible"]}
  
  # Calculate non-broken genes counts
  nonbroken_resistant <- max(0, total_resistant - broken_resistant)
  nonbroken_susceptible <- max(0, total_susceptible - broken_susceptible)
  
  # Create contingency table
  contingency_table <- matrix(
    c(broken_resistant, broken_susceptible, nonbroken_resistant, 
      nonbroken_susceptible),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Broken", "NonBroken"),
      c("Resistant", "Susceptible")
    )
  )
  
  # Print the contingency table
  cat("\nContingency table for", antibiotic_name, ":\n")
  print(contingency_table)
  
  # Create summary dataframe (broken/nonbroken + susceptible/resistant) for each antibiotic
  summary_row <- data.frame(
    Antibiotic = antibiotic_name,
    Broken_resistant = broken_resistant,
    Broken_susceptible = broken_susceptible,
    Non_broken_resistant = nonbroken_resistant,
    Non_broken_susceptible = nonbroken_susceptible,
    Total_resistant = total_resistant,
    Total_susceptible = total_susceptible,
    Total = total
  )
  
  # Perform Fisher's test and return the contingency table, summary data and test results
  tryCatch({
    result <- fisher.test(contingency_table)
    return(list(
      antibiotic = antibiotic_name,
      contingency_table = contingency_table,
      p_value = result$p.value,
      odds_ratio = result$estimate,
      confidence_interval = result$conf.int,
      summary_data = summary_row
    ))
  }, error = function(e) {
    cat("  Error in Fisher test for", antibiotic_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Main function to run tests on all antibiotics with total > 100.
run_fisher_tests <- function(general_df, combined_df) {
  antibiotics <- general_df %>%
    filter(Total > 100) %>%
    filter(resistant > 5 & susceptible > 5) %>%
    pull(Antibiotic) %>%
    unique()
  
  results <- list()
  summary_data_list <- list()
  
  for(atb in antibiotics) {
    result <- tryCatch({
      perform_fisher_test(atb, general_df, combined_df)
    }, error = function(e) {
      cat("\nError processing", atb, ":", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(result)) {
      results[[length(results) + 1]] <- result
      summary_data_list[[length(summary_data_list) + 1]] <- result$summary_data
    }
  }
  
  # Dataframe with the results
  if(length(results) > 0) {
    summary <- data.frame(
      Antibiotic = sapply(results, function(x) x$antibiotic),
      P_Value = sapply(results, function(x) x$p_value),
      Odds_Ratio = sapply(results, function(x) x$odds_ratio),
      CI_Lower = sapply(results, function(x) x$confidence_interval[1]),
      CI_Upper = sapply(results, function(x) x$confidence_interval[2]),
      stringsAsFactors = FALSE
    )
    
    # Add significance column (p < 0.05) and BH testing correction (multiple ATB)
    summary$Significant <- summary$P_Value < 0.05
    summary$Adjusted_P_Value <- p.adjust(summary$P_Value, method = "BH")
    summary$Significant_Adjusted <- summary$Adjusted_P_Value < 0.05
    
  
    full_summary <- do.call(rbind, summary_data_list)
    cat("\nFisher's Exact Tests:\n")
    print(summary)
    
    return(list(
      stats_summary = summary,
      distribution_summary = full_summary
    ))
  } 
}

results <- run_fisher_tests(general_distribution_NCBI, combined_atb)
statistics <- results$stats_summary
#write.csv(statistics, "fisher_test_results.csv", row.names = FALSE)
distribution <- results$distribution_summary
write.csv(distribution, "complete_distirbution_KO.csv", row.names = FALSE)


## Calculate Fold Change of KO detected 
total_frequencies4 <- total_frequencies3 %>%
  complete(Antibiotic, 
           Phenotype = c("resistant", "susceptible"), 
           fill = list(frequency = 0))

total_frequencies4 <- total_frequencies4[-c(3)] 

fold_data <-total_frequencies4%>%
  group_by(Antibiotic) %>%
  summarise(
    resistant_freq = frequency[Phenotype == "resistant"],
    susceptible_freq = frequency[Phenotype == "susceptible"],
    fold_change = resistant_freq / susceptible_freq,
    .groups = 'drop'
  ) %>%
  select(Antibiotic, fold_change)

mean(fold_data$fold_change[is.finite(fold_data$fold_change)])
median(fold_data$fold_change[is.finite(fold_data$fold_change)])
sd(fold_data$fold_change[is.finite(fold_data$fold_change)])

