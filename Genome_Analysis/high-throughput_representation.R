setwd("~/High-throughput/")
library(trackViewer)
library(readxl)
library(dplyr)
library(svglite)
library(ggplot2)

#Definitive script for the visualization of the mutational events conferring resistance
#in the sequenced samples of the fluctuation assay (Fig.2) and Supplementary Figure 3.
# Firstly, we plotted the genes with the higher number of events (Fig.2) and then all the genes
# affected (Supp. Fig. 3) per antibiotic. 



######################### Colistin #############################################
# Import our data resulted from analyzing the variant calling (breseq)
# of the sequenced samples of the fluctuation tests and we discard those sample
# without a mutation to explain the resistance phenotype. We also discard mgrB
# mutations because we are plotting it separately as it is a small gene and its 
# visualization can be difficult. For Figure 3, we are also showing mgrB alone.
 
evol_rescue <- read_xlsx("evolutionary_Rescue_all.xlsx", sheet = 5 )

evol_rescue_Colistin <- evol_rescue[!is.na(evol_rescue$Plot_gene) & 
                                      evol_rescue$Plot_gene != "mgrB", ]

#Detect whether the SNPs are non synonymous + non sense, synonymous or intergenic by 
#extracting the aminoacid change
evol_rescue_Colistin <- evol_rescue_Colistin %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))

# Create a single category of the mutational events, adding the information of 
# the type of NJ (whether is mediated by an IS1 or by other element).

evol_rescue_Colistin$events <- ifelse(
  is.na(evol_rescue_Colistin$IS),
  evol_rescue_Colistin$Event,
  paste(evol_rescue_Colistin$Event, evol_rescue_Colistin$IS, sep = "_")
)


# Calculation of the position of the genes for representation according to 
# the reference position.
evol_rescue_Colistin <- evol_rescue_Colistin %>%
  mutate(plot = case_when(
    Plot_gene == "qseC" ~ Position + 0,
    Plot_gene == "mgrB" ~ Position + 1469,
    Plot_gene == "crrB" ~ Position + 1700,
    Plot_gene == "phoP" ~ Position + 4236,
    Plot_gene == "phoQ" ~ Position + 2770,
    Plot_gene == "pmrB" ~ Position + 4915,
  ))

# Setting different heights or scores to differentiate between strains and conditions.
evol_rescue_Colistin$score <- NA
evol_rescue_Colistin <- evol_rescue_Colistin %>%
  mutate(score = case_when(
    STRAIN == "KPN08" ~ 31,
    STRAIN == "KPN16" ~ 1
    
  ))

positions <- evol_rescue_Colistin %>%
  pull(plot)
events <- evol_rescue_Colistin %>%
  pull(events)
strain <- evol_rescue_Colistin %>%
  pull(STRAIN)
condition <- evol_rescue_Colistin %>%
  pull(Condition)
is <- evol_rescue_Colistin %>%
  pull(IS)
score <- evol_rescue_Colistin %>%
  pull(score)

## Set color and side for the events
color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)  

#Initialize the object and apply all the parameters
mut_gr <- GRanges("chr1", IRanges(positions, width = 1))

side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45

## Description of the start and end for each gen in order to represent them.
qsec_start <- 0
qsec_end <- 1349
arrow1 <- GRanges("chr1", ranges = IRanges(start = qsec_start, end = qsec_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")
crrb_start <- 1700
crrb_end <- 2761
arrow2 <- GRanges("chr1", ranges = IRanges(start = crrb_start, end = crrb_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")
phoQ_start <- 2770
phoQ_end <- 4236
arrow3 <- GRanges("chr1", ranges = IRanges(start = phoQ_start, end = phoQ_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

phoP_start <- 4236
phoP_end <- 4907
arrow4 <- GRanges("chr1", ranges = IRanges(start = phoP_start, end = phoP_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

pmrb_start <- 4915
pmrb_end <- 6012
arrow5 <- GRanges("chr1", ranges = IRanges(start = pmrb_start, end = pmrb_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

# Combine all the arrows as features and plotting
features <- c(arrow1,  arrow2,  arrow3, arrow4, arrow5)
svg("ER_COL_except_mgrb.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 6020)), 
          cex = 0.8)
dev.off()




### Next, we plot mgrB independently 
evol_rescue_Colistin_mgrB <- evol_rescue[!is.na(evol_rescue$Plot_gene) & 
                                           evol_rescue$Plot_gene == "mgrB", ]

evol_rescue_Colistin_mgrB <- evol_rescue_Colistin_mgrB %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


evol_rescue_Colistin_mgrB$events <- ifelse(
  is.na(evol_rescue_Colistin_mgrB$IS),
  evol_rescue_Colistin_mgrB$Event,
  paste(evol_rescue_Colistin_mgrB$Event, evol_rescue_Colistin_mgrB$IS, sep = "_")
)

evol_rescue_Colistin_mgrB$score <- NA
evol_rescue_Colistin_mgrB$plot <- evol_rescue_Colistin_mgrB$Position+ 150


evol_rescue_Colistin_mgrB <- evol_rescue_Colistin_mgrB %>%
  mutate(score = case_when(
    STRAIN == "KPN08" ~ 31,
    STRAIN == "KPN16" ~ 1
  ))

positions <- evol_rescue_Colistin_mgrB %>%
  pull(plot)
events <- evol_rescue_Colistin_mgrB %>%
  pull(events)
strain <- evol_rescue_Colistin_mgrB %>%
  pull(STRAIN)
condition <- evol_rescue_Colistin_mgrB %>%
  pull(Condition)
is <- evol_rescue_Colistin_mgrB %>%
  pull(IS)
score <- evol_rescue_Colistin_mgrB %>%
  pull(score)

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)  


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))

side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45

mgr_start <- 150
mgr_end <- 294
arrow_mgr <- GRanges("chr1", ranges = IRanges(start = mgr_start, end = mgr_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")


features <- c(arrow_mgr)
svg("ER_COL_mgrb.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 300)), 
          cex = 0.8)
dev.off()


##### We plot all the mutational events affecting all the strains,
#independently of parallel evolution. 
muts_counter <- read_xlsx("evolutionary_Rescue_all.xlsx", sheet = 5 )
muts_counter <- muts_counter[!is.na(muts_counter$Event), ]

muts_counter <- muts_counter %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


muts_counter$events <- ifelse(
  is.na(muts_counter$IS),
  muts_counter$Event,
  paste(muts_counter$Event, muts_counter$IS, sep = "_")
)


muts_PF <- muts_counter[muts_counter$Condition=="PF", ]
muts_TC <- muts_counter[muts_counter$Condition=="TC", ]



counter_PF <- muts_PF %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_PF$events <- factor(counter_PF$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", "NJ_Y",
                                                          "SNP"))


PF_COL <- ggplot(counter_PF, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") + 
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "DEL" = "#575757ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_COL_PF.svg", plot = PF_COL)



counter_TC <- muts_TC %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_TC$events <- factor(counter_TC$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", 
                                                          "NJ_Y"))



TC_COL <- ggplot(counter_TC, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") +
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "DEL" = "#575757ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_COL_TC.svg", plot = TC_COL)





################################ Ciprofloxacin #########################################
evol_rescue_cipro <- read_xlsx("evolutionary_Rescue_all.xlsx", sheet = 6 )
evol_rescue_cipro <- evol_rescue_cipro %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


evol_rescue_cipro$events <- ifelse(
  is.na(evol_rescue_cipro$IS),
  evol_rescue_cipro$Event,
  paste(evol_rescue_cipro$Event, evol_rescue_cipro$IS, sep = "_")
)

evol_rescue_cipro <- evol_rescue_cipro%>%
  mutate(score = case_when(
    STRAIN == "KPN08" ~ 41, 
    STRAIN == "KPN16" ~ 11
  ))

#### Due to the high number of genes with mutations, we separate them in fragments to visualize easily: 
#First fragment
evol_rescue_1<- evol_rescue_cipro[evol_rescue_cipro$Plot_gene 
                                       %in% c("LysR", "artP", 
                                              "ABC_periplasmic_component",
                                              "ABC_permease_component",
                                              "ABC_binding_protein",
                                              "IM_protein", "Metal_efflux_pump",
                                              "Sugar_permease"), ]

evol_rescue_1 <- evol_rescue_1 %>%
  mutate(Position = as.numeric(Position)) %>%
  mutate(plot = case_when(
    Plot_gene == "LysR" ~ Position + 0,
    Plot_gene == "artP" ~ Position + 946,
    Plot_gene == "ABC_periplasmic_component" ~ Position + 1750,
    Plot_gene == "ABC_permease_component" ~ Position + 2696,
    Plot_gene == "ABC_binding_protein" ~ Position + 4112, 
    Plot_gene == "IM_protein" ~ Position + 4959,
    Plot_gene == "Metal_efflux_pump" ~ Position + 5608,
    Plot_gene == "Sugar_permease" ~ Position + 6536
  ))

positions <- evol_rescue_1 %>%
  pull(plot)
events <- evol_rescue_1 %>%
  pull(events)
strain <- evol_rescue_1 %>%
  pull(STRAIN)
condition <- evol_rescue_1 %>%
  pull(Condition)
is <- evol_rescue_1 %>%
  pull(IS)
score <- evol_rescue_1 %>%
  pull(score)

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45


LysR_start <- 0
LysR_end <- 941
arrow <- GRanges("chr1", ranges = IRanges(start = LysR_start, end = LysR_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

artP_start <- 946
artP_end<- 1745
arrow1 <- GRanges("chr1", ranges = IRanges(start = artP_start, end = artP_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

ABCp_start <- 1750
ABCp_end <- 2691
arrow2 <- GRanges("chr1", ranges = IRanges(start = ABCp_start, end = ABCp_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")
ABCper_start <- 2696
ABCper_end <- 3723
arrow3 <- GRanges("chr1", ranges = IRanges(start = ABCper_start, end = ABCper_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

ABCb_start <- 4112 
ABCb_end <- 4954
arrow4 <- GRanges("chr1", ranges = IRanges(start = ABCb_start, end = ABCb_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")
IM_start <- 4959
IM_end <- 5603
arrow5 <- GRanges("chr1", ranges = IRanges(start = IM_start, end = IM_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

Metaleff_start <- 5608
Metaleff_end <- 6531
arrow6 <- GRanges("chr1", ranges = IRanges(start = Metaleff_start, 
                                           end = Metaleff_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

sugarperm_start <- 6536
sugarperm_end <- 7864
arrow7 <- GRanges("chr1", ranges = IRanges(start = sugarperm_start, 
                                           end = sugarperm_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

features <- c(arrow,  arrow1, arrow2, arrow3, arrow4, arrow5, arrow6, arrow7)
svg("ER_CIP_part1.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 7870)), 
          cex = 0.8)
dev.off()



##### Second fragment
evol_rescue_2<- evol_rescue_cipro[evol_rescue_cipro$Plot_gene 
                                  %in% c("aaeB", 
                                         "ABC_permease_component_2",
                                         "oqxB",
                                         "MFS_efflux_pump",
                                         "membrane_protein"), ]

evol_rescue_2 <- evol_rescue_2 %>%
  mutate(Position = as.numeric(Position)) %>%
  mutate(plot = case_when(
    Plot_gene == "aaeB" ~ Position + 0,
    Plot_gene == "ABC_permease_component_2" ~ Position + 1972,
    Plot_gene == "oqxB" ~Position + 2867,
    Plot_gene == "MFS_efflux_pump" ~ Position + 6024,
    Plot_gene == "membrane_protein" ~ Position + 7225
))

positions <- evol_rescue_2 %>%
  pull(plot)
events <- evol_rescue_2 %>%
  pull(events)
strain <- evol_rescue_2 %>%
  pull(STRAIN)
condition <- evol_rescue_2 %>%
  pull(Condition)
is <- evol_rescue_2 %>%
  pull(IS)
score <- evol_rescue_2 %>%
  pull(score)

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45


aaeB_start <- 0
aaeB_end <- 1967
arrow <- GRanges("chr1", ranges = IRanges(start = aaeB_start, end = aaeB_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABCper2_start <- 1972
ABCper2_end<- 2862
arrow1 <- GRanges("chr1", ranges = IRanges(start = ABCper2_start, end = ABCper2_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

oqxB_start <- 2867
oqxB_end <- 6019
arrow2 <- GRanges("chr1", ranges = IRanges(start = oqxB_start, end = oqxB_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

MFS_eff_start <- 6024
MFS_eff_end <- 7220
arrow3 <- GRanges("chr1", ranges = IRanges(start = MFS_eff_start, end = MFS_eff_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")
MP_start <- 7225
MP_end <- 8274
arrow4 <- GRanges("chr1", ranges = IRanges(start = MP_start, end = MP_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")


features <- c(arrow,  arrow1, arrow2, arrow3, arrow4)
svg("ER_CIP_part2.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 8280)), 
          cex = 0.8)
dev.off()


### Third fragment
evol_rescue_3<- evol_rescue_cipro[evol_rescue_cipro$Plot_gene 
                                  %in% c("ABC_permease_component_3", 
                                         "ABC_periplasmic_component_1",
                                         "ABC_permease_component_4",
                                         "ABC_ATP_component",
                                         "MFS_transporter",
                                         "HpxR", "marR", "LysR_1", "LysR_2",
                                         "ompA"), ]

evol_rescue_3 <- evol_rescue_3 %>%
  mutate(Position = as.numeric(Position)) %>%
  mutate(plot = case_when(
    Plot_gene == "ABC_permease_component_3" ~ Position + 0,
    Plot_gene == "ABC_periplasmic_component_1" ~ Position + 1291,
    Plot_gene == "ABC_permease_component_4" ~ Position + 2552,
    Plot_gene == "ABC_ATP_component" ~ Position + 3506,
    Plot_gene == "MFS_transporter" ~ Position + 4524, 
    Plot_gene == "HpxR" ~ Position + 6369,
    Plot_gene == "marR" ~ Position + 7519,
    Plot_gene == "LysR_1" ~ Position + 7958,
    Plot_gene == "LysR_2" ~ Position + 8865, 
    Plot_gene == "ompA" ~ Position + 9775
  ))

positions <- evol_rescue_3 %>%
  pull(plot)
events <- evol_rescue_3 %>%
  pull(events)
strain <- evol_rescue_3 %>%
  pull(STRAIN)
condition <- evol_rescue_3 %>%
  pull(Condition)
is <- evol_rescue_3 %>%
  pull(IS)
score <- evol_rescue_3 %>%
  pull(score)

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45


ABCper3_start <- 0
ABCper3_end <- 1286
arrow <- GRanges("chr1", ranges = IRanges(start = ABCper3_start, end = ABCper3_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABCp1_start <- 1291
ABCp1_end<- 2547
arrow1 <- GRanges("chr1", ranges = IRanges(start = ABCp1_start, end = ABCp1_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

ABCper4_start <- 2552
ABCper4_end <- 3301
arrow2 <- GRanges("chr1", ranges = IRanges(start = ABCper4_start, end = ABCper4_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

ABC_ATP_start <- 3506
ABC_ATP_end <- 4519
arrow3 <- GRanges("chr1", ranges = IRanges(start = ABC_ATP_start, end = ABC_ATP_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

MFSt_start <- 4524
MFSt_end <- 5930
arrow4 <- GRanges("chr1", ranges = IRanges(start = MFSt_start, end = MFSt_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

hpxR_start <- 6369
hpxR_end <- 7292
arrow5 <- GRanges("chr1", ranges = IRanges(start = hpxR_start, end = hpxR_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

marR_start <- 7519
marR_end <- 7953
arrow6 <- GRanges("chr1", ranges = IRanges(start = marR_start, end = marR_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

LysR1_start <- 7958
LysR1_end <- 8860
arrow7 <- GRanges("chr1", ranges = IRanges(start = LysR1_start, end = LysR1_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

LysR2_start <- 8865
LysR2_end <- 9770
arrow8 <- GRanges("chr1", ranges = IRanges(start = LysR2_start, end = LysR2_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")

ompA_start <- 9775
ompA_end <- 11481
arrow9 <- GRanges("chr1", ranges = IRanges(start = ompA_start, end = ompA_end),
                  height = 0.05, shape = "arrow", fill = "lightgrey")


features <- c(arrow,  arrow1, arrow2, arrow3, arrow4, arrow5, arrow6, arrow7, 
              arrow8, arrow9)

svg("ER_CIP_part3.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 11490)), 
          cex = 0.8)
dev.off()


##### Then we plot the most relevant genes: oqxR y ramR
evol_rescue_cipro <- evol_rescue_cipro[evol_rescue_cipro$Plot_gene 
                                       %in% c("oqxR", "ramR"), ]


evol_rescue_cipro <- evol_rescue_cipro %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


evol_rescue_cipro$events <- ifelse(
  is.na(evol_rescue_cipro$IS),
  evol_rescue_cipro$Event,
  paste(evol_rescue_cipro$Event, evol_rescue_cipro$IS, sep = "_")
)

evol_rescue_cipro <- evol_rescue_cipro%>%
  mutate(score = case_when(
    STRAIN == "KPN08" ~ 41, 
    STRAIN == "KPN16" ~ 11
  ))



evol_rescue_cipro <- evol_rescue_cipro %>%
  mutate(Position = as.numeric(Position)) %>%
  mutate(plot = case_when(
    Plot_gene == "oqxR" ~ Position + 57,
    Plot_gene == "ramR" ~ Position + 567,
  ))

positions <- evol_rescue_cipro %>%
  pull(plot)
events <- evol_rescue_cipro %>%
  pull(events)
strain <- evol_rescue_cipro %>%
  pull(STRAIN)
condition <- evol_rescue_cipro %>%
  pull(Condition)
is <- evol_rescue_cipro %>%
  pull(IS)
score <- evol_rescue_cipro %>%
  pull(score)

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45


oqxR_start <- 27
oqxR_end <- 506
arrow <- GRanges("chr1", ranges = IRanges(start = oqxR_start, end = oqxR_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")
ramR_start <- 547
ramR_end<- 1128
arrow2 <- GRanges("chr1", ranges = IRanges(start = ramR_start, end = ramR_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")


features <- c(arrow,  arrow2)
svg("ER_CIP_main.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 1200)), 
          cex = 0.8)
dev.off()

# We plot all the mutational evetns affecting all the strains,
#independently of parallel evolution. 
evol_rescue <- read_xlsx("evolutionary_Rescue_all.xlsx", sheet = 6 )
muts_counter <- evol_rescue[!is.na(evol_rescue$Event), ]

muts_counter <- muts_counter %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


muts_counter$events <- ifelse(
  is.na(muts_counter$IS),
  muts_counter$Event,
  paste(muts_counter$Event, muts_counter$IS, sep = "_")
)


muts_PF <- muts_counter[muts_counter$Condition=="PF", ]
muts_TC <- muts_counter[muts_counter$Condition=="TC", ]



counter_PF <- muts_PF %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_PF$events <- factor(counter_PF$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "SNP", "NJ_N", 
                                                          "NJ_Y"))



PF <- ggplot(counter_PF, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") + 
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "DEL" = "#575757ff",
                               "SNP" = "#b7c4c8ff", 
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_CIP_PF.svg", plot = PF)



counter_TC <- muts_TC %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_TC$events <- factor(counter_TC$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", "NJ_Y",
                                                          "SNP"))



TC <- ggplot(counter_TC, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") +
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "DEL" = "#575757ff",
                               "SNP" = "#b7c4c8ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_CIP_TC.svg", plot = TC)



############################### Chloramphenicol ######################################## 
evol_rescue_chloram <- read_xlsx("evolutionary_Rescue_all.xlsx", sheet = 7)
evol_rescue_chloram <- evol_rescue_chloram[!evol_rescue_chloram$Plot_gene =="oqxR", ]

evol_rescue_chloram <- evol_rescue_chloram[!is.na(evol_rescue_chloram$Plot_gene), ]

evol_rescue_chloram <- evol_rescue_chloram %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


evol_rescue_chloram$events <- ifelse(
  is.na(evol_rescue_chloram$IS),
  evol_rescue_chloram$Event,
  paste(evol_rescue_chloram$Event, evol_rescue_chloram$IS, sep = "_")
)


evol_rescue_chloram <- evol_rescue_chloram %>%
  mutate(plot = case_when(
    Plot_gene == "ArgT" ~ Position + 534,
    Plot_gene == "AraJ" ~ Position + 1321
  ))


positions <- evol_rescue_chloram %>%
  pull(plot)
events <- evol_rescue_chloram %>%
  pull(events)
strain <- evol_rescue_chloram %>%
  pull(STRAIN)
condition <- evol_rescue_chloram %>%
  pull(Condition)
is <- evol_rescue_chloram %>%
  pull(IS)

mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$label.parameter.rot <- 45


argT_start <- 534
argT_end <- 1316
arrow1 <- GRanges("chr1", ranges = IRanges(start = argT_start, end = argT_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

araJ_start <- 1321
araJ_end <- 3577
arrow3 <- GRanges("chr1", ranges = IRanges(start = araJ_start, end = araJ_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

features <- c(arrow1, arrow2, arrow3)
svg("ERC_CMP.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(530, 3580)), 
          cex = 0.8)
dev.off()

## We plot individually oqxR, as it is the main gene 
evol_rescue_chloram <- read_xlsx("evolutionary_Rescue_all.xlsx", sheet = 7)
evol_rescue_chloram <- evol_rescue_chloram[evol_rescue_chloram$Plot_gene =="oqxR", ]

evol_rescue_chloram <- evol_rescue_chloram %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


evol_rescue_chloram$events <- ifelse(
  is.na(evol_rescue_chloram$IS),
  evol_rescue_chloram$Event,
  paste(evol_rescue_chloram$Event, evol_rescue_chloram$IS, sep = "_")
)


evol_rescue_chloram <- evol_rescue_chloram %>%
  mutate(plot = case_when(
    Plot_gene == "oqxR" ~ Position + 50
  ))


positions <- evol_rescue_chloram %>%
  pull(plot)
events <- evol_rescue_chloram %>%
  pull(events)
strain <- evol_rescue_chloram %>%
  pull(STRAIN)
condition <- evol_rescue_chloram %>%
  pull(Condition)
is <- evol_rescue_chloram %>%
  pull(IS)

mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$label.parameter.rot <- 45
oqxR_start <- 50
oqxR_end <- 529
arrow <- GRanges("chr1", ranges = IRanges(start = oqxR_start, end = oqxR_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")
features <- c(arrow)
svg("ER_CMP_main.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 550)), 
          cex = 0.8)
dev.off()


muts_counter <- read_xlsx("evolutionary_Rescue_all.xlsx", sheet = 7 )
muts_counter <- muts_counter[!is.na(muts_counter$Event), ]

muts_counter <- muts_counter %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


muts_counter$events <- ifelse(
  is.na(muts_counter$IS),
  muts_counter$Event,
  paste(muts_counter$Event, muts_counter$IS, sep = "_")
)


muts_PF <- muts_counter[muts_counter$Condition=="PF", ]
muts_TC <- muts_counter[muts_counter$Condition=="TC", ]


counter_PF <- muts_PF %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_PF$events <- factor(counter_PF$events, levels = c("DEL", "INDEL", "SNP_I",
                                                          "SNP_N", "NJ_N", "NJ_Y",
                                                          "SNP"))


PF <- ggplot(counter_PF, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") + 
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "SNP" = "#b7c4c8ff",
                               "DEL" = "#575757ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_CMP_PF.svg", plot = PF)


counter_TC <- muts_TC %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_TC$events <- factor(counter_TC$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", "NJ_Y",
                                                          "SNP"))



TC <- ggplot(counter_TC, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") + 
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff",
                               "SNP" = "#b7c4c8ff",
                               "DEL" = "#575757ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_CMP_TC.svg", plot = TC)


################################ Kanamycin  ####################################

evol_rescue_kana <- read_xlsx("evolutionary_Rescue_all.xlsx",sheet = 8)
evol_rescue_kana <- evol_rescue_kana[!evol_rescue_kana$Plot_gene == "sbmA", ]

evol_rescue_kana <- evol_rescue_kana %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


evol_rescue_kana$events <- ifelse(
  is.na(evol_rescue_kana$IS),
  evol_rescue_kana$Event,
  paste(evol_rescue_kana$Event, evol_rescue_kana$IS, sep = "_")
)

evol_rescue_kana <- evol_rescue_kana %>%
  mutate(plot = case_when(
    Plot_gene == "ABC_periplasmic" ~ Position + 0,
    Plot_gene == "MD_efflux_RND" ~ Position + 985,
    Plot_gene == "ABC_periplasmic_2" ~ Position + 4097,
    Plot_gene == "ABC_permease_LivM" ~ Position + 5091,
    Plot_gene == "ABC_transmembrane_protein" ~ Position + 6373,
    Plot_gene == "ABC_permease_component" ~ Position + 7055,
    Plot_gene == "ABC_transporter_protein" ~ Position + 8537,
    Plot_gene == "ABC_periplasmic_3" ~ Position + 9126,
    Plot_gene == "ABC_permease_component_2" ~ Position + 10099,
    Plot_gene == "IM_protein_2" ~ Position + 11759,
  ))


positions <- evol_rescue_kana %>%
  pull(plot)
events <- evol_rescue_kana %>%
  pull(events)
strain <- evol_rescue_kana %>%
  pull(STRAIN)
condition <- evol_rescue_kana %>%
  pull(Condition)
is <- evol_rescue_kana %>%
  pull(IS)


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)


mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$label.parameter.rot <- 45


ABC_periplasmic_start <- 0
ABC_periplasmic_end <- 980
arrow <- GRanges("chr1", ranges = IRanges(start = ABC_periplasmic_start, 
                                          end = ABC_periplasmic_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

MD_efflux_start <- 985
MD_efflux_end <- 4092
arrow2 <- GRanges("chr1", ranges = IRanges(start = MD_efflux_start, 
                                          end = MD_efflux_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABC_periplasmic2_start <- 4097
ABC_periplasmic2_end <- 5086
arrow3 <- GRanges("chr1", ranges = IRanges(start = ABC_periplasmic2_start, 
                                          end = ABC_periplasmic2_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

livM_start <- 5091
livM_end <- 6368
arrow4 <- GRanges("chr1", ranges = IRanges(start = livM_start, 
                                          end = livM_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABC_transmembrane_start <- 6373
ABC_transmembrane_end <- 7050
arrow5 <- GRanges("chr1", ranges = IRanges(start = ABC_transmembrane_start, 
                                          end = ABC_transmembrane_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABC_permease_start <- 7055
ABC_permease_end <- 8050
arrow6 <- GRanges("chr1", ranges = IRanges(start = ABC_permease_start, 
                                          end = ABC_permease_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABC_transporter_start <- 8537
ABC_transporter_end <- 9121
arrow7 <- GRanges("chr1", ranges = IRanges(start = ABC_transporter_start, 
                                          end = ABC_transporter_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABC_periplasmic3_start <- 9126
ABC_periplasmic3_end <- 10094
arrow8 <- GRanges("chr1", ranges = IRanges(start = ABC_periplasmic3_start, 
                                          end = ABC_periplasmic3_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

ABC_permease2_start <- 10099
ABC_permease2_end <- 11754
arrow9 <- GRanges("chr1", ranges = IRanges(start = ABC_permease2_start, 
                                          end = ABC_permease2_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

IM_protein_start <- 11759
IM_protein_end <- 11998
arrow1 <- GRanges("chr1", ranges = IRanges(start = IM_protein_start, 
                                          end = IM_protein_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

features <- c(arrow, arrow2, arrow3, arrow4, arrow5, arrow6, arrow7, arrow8, 
              arrow9, arrow1)
svg("ER_KAN_genes.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 12000)), 
          cex = 0.8)
dev.off()


## We plot sbmA individually
evol_rescue_kana <- read_xlsx("evolutionary_Rescue_all.xlsx",sheet = 8)
evol_rescue_kana <- evol_rescue_kana[evol_rescue_kana$Plot_gene == "sbmA", ]

evol_rescue_kana <- evol_rescue_kana %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


evol_rescue_kana$events <- ifelse(
  is.na(evol_rescue_kana$IS),
  evol_rescue_kana$Event,
  paste(evol_rescue_kana$Event, evol_rescue_kana$IS, sep = "_")
)


positions <- evol_rescue_kana %>%
  pull(Position)
events <- evol_rescue_kana %>%
  pull(events)
strain <- evol_rescue_kana %>%
  pull(STRAIN)
condition <- evol_rescue_kana %>%
  pull(Condition)
is <- evol_rescue_kana %>%
  pull(IS)


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP_I", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$label.parameter.rot <- 45


sbmA_start <- 0
sbmA_end <- 1220
arrow <- GRanges("chr1", ranges = IRanges(start = sbmA_start, end = sbmA_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

features <- c(arrow)
svg("ER_KAN_main.svg")
lolliplot(mut_gr, features, 
          ranges = GRanges("chr1", IRanges(0, 1260)), 
          cex = 0.8)
dev.off()


evol_rescue_kana <- read_xlsx("evolutionary_Rescue_all.xlsx",sheet = 8)
evol_rescue_kana <- evol_rescue_kana[!is.na(evol_rescue_kana$Event), ]

muts_counter <- evol_rescue_kana %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation) ~ {
      amino_change <- sub("\\(.*$", "", Annotation)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


muts_counter$events <- ifelse(
  is.na(muts_counter$IS),
  muts_counter$Event,
  paste(muts_counter$Event, muts_counter$IS, sep = "_")
)


muts_PF <- muts_counter[muts_counter$Condition=="PF", ]
muts_TC <- muts_counter[muts_counter$Condition=="TC", ]



counter_PF <- muts_PF %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_PF$events <- factor(counter_PF$events, levels = c("DEL", "INDEL", "SNP_I",
                                                          "SNP_N", "SNP", "NJ_N", 
                                                          "NJ_Y"))



PF <- ggplot(counter_PF, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") + 
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "SNP" = "#b7c4c8ff",
                               "DEL" = "#575757ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_KAN_PF.svg", plot = PF)



counter_TC <- muts_TC %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_TC$events <- factor(counter_TC$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", "NJ_Y",
                                                          "SNP"))



TC <- ggplot(counter_TC, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") + 
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "SNP" = "#b7c4c8ff",
                               "DEL" = "#575757ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_ER_KAN_TC.svg", plot = TC)
