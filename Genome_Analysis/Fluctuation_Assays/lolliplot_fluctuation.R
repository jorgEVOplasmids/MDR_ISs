#Definitive script for the visualization of the mutational events conferring resistance
#in the sequenced samples of the fluctuation assay (Fig.1) and Supplementary Figure 1

setwd("~/Fluctuation_analysis")
library(trackViewer)
library(readxl)
library(dplyr)
library(svglite)
library(ggplot2)


# Import our data resulted from analyzing the variant calling (breseq)
#of the sequenced samples of the fluctuation tests and we discard those sample
#without a mutation to explain the resistance phenotype.

muts_strain_TC_PF <- read_xlsx("fluctuation_test_all.xlsx", sheet = 3 )
muts_strain_TC_PF <- muts_strain_TC_PF[!is.na(muts_strain_TC_PF$Plot_gene), ]


#Detect whether the SNPs are non synonymous + non sense, synonymous or intergenic by 
#extracting the aminoacid change 

muts_strain_TC_PF <- muts_strain_TC_PF %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation_1) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation_1) ~ {
      amino_change <- sub("\\(.*$", "", Annotation_1)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))

# For better visualization, we only plot the genes with parallel evolution in Figure 1. Therefore, we
# exclude the "other" genes. Also, due to the small size of mgrB (144 bp), it will be plotted
# independently

mutsP <- muts_strain_TC_PF
mutsP <- mutsP[mutsP$Plot_gene != "other", ]
mutsP <- mutsP[mutsP$Plot_gene != "mgrB", ]

# Create a single category of the mutational events, adding the information of 
# the type of NJ (whether is mediated by an IS1 or by other element).
mutsP$events <- ifelse(
  is.na(mutsP$IS),
  mutsP$Event,                    
  paste(mutsP$Event, mutsP$IS, sep = "_")  
)

# Calculation of the position of the genes for representation according to 
# the reference position
mutsP <- mutsP %>%
  mutate(plot = case_when(
    Plot_gene == "arnB" ~ Position + 97,
    Plot_gene == "envZ" ~ Position + 1239,
    Plot_gene == "qseC" ~ Position + 2599,
    Plot_gene == "lpxM" ~ Position + 3953,
    Plot_gene == "dsbB" ~ Position + 4932,
    Plot_gene == "crrB" ~ Position + 5467,
    Plot_gene == "phoQ" ~ Position + 6533,
    Plot_gene == "phoP" ~ Position + 7999,
    Plot_gene == "pmrA" ~ Position + 8675,
    Plot_gene == "pmrB" ~ Position + 9346,
    Plot_gene == "dsbA" ~ Position + 10448
  ))

# Setting different heights or scores to differentiate between strains and conditions.
mutsP$score <- NA
mutsP <- mutsP %>%
  mutate(score = case_when(
    STRAIN == "KPN01" ~ 1,
    STRAIN == "KPN10" ~ 51,
    STRAIN == "KPN13" ~ 101,
    STRAIN == "KPN16" ~ 151,
    STRAIN == "KPN08" ~ 201, 
    STRAIN == "Delta_KPN01" ~ 251,
    STRAIN == "Delta_KPN10" ~ 301,
    STRAIN == "Delta_KPN13" ~ 351,
    STRAIN == "Delta_KPN16" ~ 401, 
    STRAIN == "Delta_KPN08" ~ 451
  ))

positions <- mutsP %>%
  pull(plot)
events <- mutsP %>%
  pull(events)
strain <- mutsP %>%
  pull(STRAIN)
condition <- mutsP %>%
  pull(Condition)
is <- mutsP %>%
  pull(IS)
score <- mutsP %>%
  pull(score)

## Set color and side for the events
color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)


side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

#Initialize the object and apply all the parameters
mut_gr <- GRanges("chr1", IRanges(positions, width = 1))

mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45


## Description of the start and end for each gen in order to represent them.
arnB_start<- 1
arnB_end<- 1236
arrow <- GRanges("chr1", ranges = IRanges(start = arnB_start, end = arnB_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")
envz_start<- 1239
envz_end<- 2594
arrow1 <- GRanges("chr1", ranges = IRanges(start = envz_start, end = envz_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")
qsec_start <- 2599
qsec_end <- 3948
arrow2 <- GRanges("chr1", ranges = IRanges(start = qsec_start, end = qsec_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")
lpxM_start <- 3953
lpxM_end<- 4927
arrow3 <- GRanges("chr1", ranges = IRanges(start = lpxM_start, end = lpxM_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")
dsbb_start <- 4932
dsbb_end <- 5462
arrow4 <- GRanges("chr1", ranges = IRanges(start = dsbb_start, end = dsbb_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

crrb_start <- 5467
crrb_end <- 6528
arrow5 <- GRanges("chr1", ranges = IRanges(start = crrb_start, end = crrb_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

phoQ_start <- 6533
phoQ_end <- 7999
arrow6 <- GRanges("chr1", ranges = IRanges(start = phoQ_start, end = phoQ_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

phoP_start <- 8000
phoP_end <- 8670
arrow7 <- GRanges("chr1", ranges = IRanges(start = phoP_start, end = phoP_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")
pmra_start <- 8675
pmra_End <- 9346
arrow8 <- GRanges("chr1", ranges = IRanges(start = pmra_start, end = pmra_End), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

pmrb_start <- 9347
pmrb_end <- 10443
arrow9 <- GRanges("chr1", ranges = IRanges(start = pmrb_start, end = pmrb_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")
dsba_start <- 10448
dsba_end <- 11071
arrow10 <- GRanges("chr1", ranges = IRanges(start = dsba_start, end = dsba_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

# Combine all the arrows as features and plotting
features <- c(arrow, arrow1, arrow2, arrow3, arrow4, arrow5, arrow6, arrow7, arrow8,
              arrow9, arrow10)

svg("fluctuationtest_withoutmgrb.svg")
lolliplot(mut_gr, features = features,
          ranges = GRanges("chr1", IRanges(0, 11072)), cex= 1)

dev.off()


#### We repeat the same process plotting only mgrB.
mgrB <- muts_strain_TC_PF[muts_strain_TC_PF$Plot_gene=="mgrB", ]
mgrB$Plot <- mgrB$Position+ 288

mgrB$event <- ifelse(
  is.na(mgrB$IS),
  mgrB$Event,
  paste(mgrB$Event, mgrB$IS, sep = "_")
)

mgrB$score <- NA
mgrB <- mgrB %>%
  mutate(score = case_when(
    STRAIN == "KPN01" ~ 1,
    STRAIN == "KPN10" ~ 51,
    STRAIN == "KPN13" ~ 101,
    STRAIN == "KPN16" ~ 151,
    STRAIN == "KPN08" ~ 201, 
    STRAIN == "Delta_KPN01" ~ 251,
    STRAIN == "Delta_KPN10" ~ 301,
    STRAIN == "Delta_KPN13" ~ 351,
    STRAIN == "Delta_KPN16" ~ 401, 
    STRAIN == "Delta_KPN08" ~ 451
    
  ))

positions <- mgrB %>%
  pull(Plot)
events <- mgrB %>%
  pull(event)
strain <- mgrB %>%
  pull(STRAIN)
condition <- mgrB %>%
  pull(Condition)
is <- mgrB %>%
  pull(IS)
score <- mgrB %>%
  pull(score)

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_Y", "#810f7c", color) 
color <- ifelse (events == "NJ_N", "#c6a8dfff", color) 
color <- ifelse(events == "SNP", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)  

side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_grs <- GRanges("chr1", IRanges(positions, width = 1))


side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")

mut_grs$color <- color
mut_grs$SNPsideID <-side
mut_grs$label.parameter.rot <- 45
mut_grs$score <- score

mgr_start <- 288
mgr_end <- 432


arrow_mgr <- GRanges("chr1", ranges = IRanges(start = mgr_start, end = mgr_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")


features <-  arrow_mgr
svg("lolliplot_mgrB_FLUCTUATION.svg")
lolliplot(mut_grs, features = features,
          ranges = GRanges("chr1", IRanges(0, 500)),
          cex = 0.8 )
dev.off()




##### We plot all the mutational evetns affecting all the strains,
#independently of parallel evolution. 
muts_counter<- read_xlsx("fluctuation_test_all.xlsx", sheet = 3 )
muts_counter <- muts_counter[!is.na(muts_counter$Event), ]

muts_counter <- muts_counter %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation_1) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation_1) ~ {
      amino_change <- sub("\\(.*$", "", Annotation_1)
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

muts_counter$Condition[grepl("^Delta_", muts_counter$STRAIN)] <- "Delta"


muts_PF <- muts_counter[muts_counter$Condition=="PF", ]
muts_TC <- muts_counter[muts_counter$Condition=="TC", ]
muts_Delta <- muts_counter[muts_counter$Condition == "Delta", ]


counter_PF <- muts_PF %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_PF$events <- factor(counter_PF$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", 
                                                          "NJ_Y"))


PF_strains <- ggplot(counter_PF, aes(x = "", y = Count, fill = events)) +
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

ggsave("events_fluctuation_PFstrains.svg", plot = PF_Strains)



counter_TC <- muts_TC %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_TC$events <- factor(counter_TC$events, levels = c("DEL", "INDEL", 
                                                          "SNP_I", "SNP_N", 
                                                          "NJ_N","NJ_Y"))


TC_strains <- ggplot(counter_TC, aes(x = "", y = Count, fill = events)) +
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

ggsave("events_fluctuation_TCstrains.svg", plot = TC_strains)



counter_delta <- muts_Delta %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_delta$events <- factor(counter_delta$events, levels = c("DEL", "INDEL", 
                                                          "SNP_I", "SNP_N", 
                                                          "NJ_N","NJ_Y"))

counter_delta$events <- factor(
  counter_delta$events,
  levels = c("NJ_Y", "NJ_N", "SNP_N", "SNP_I", "INDEL", "DEL")
)

Delta_strains <- ggplot(counter_delta, aes(x = "", y = Count, fill = events)) +
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

ggsave("events_fluctuation_Deltastrains.svg", plot = Delta_strains)





### For the supplementary Figure 1, we plot all the genes including those without parallel evolution. 
# For this, we are plotting the missing genes and adding them to the final svg. 
#We also add the information of the Type_IS column, differenciating the IS5.
#We also repeated the piecharts adding all the information of the type of ISs. 
muts_strain_TC_PF <- read_xlsx("fluctuation_test_all.xlsx", sheet = 3 )
muts_strain_TC_PF <- muts_strain_TC_PF[!is.na(muts_strain_TC_PF$Plot_gene), ]


muts_strain_TC_PF <- muts_strain_TC_PF %>%
  mutate(Event = case_when(
    Event != "SNP" ~ Event,
    Event == "SNP" & grepl("^intergenic", Annotation_1) ~ "SNP_I",
    Event == "SNP" & !grepl("^intergenic", Annotation_1) ~ {
      amino_change <- sub("\\(.*$", "", Annotation_1)
      first_aa <- substr(amino_change, 1, 1)
      last_aa <- substr(amino_change, nchar(amino_change), nchar(amino_change))
      ifelse(first_aa == last_aa, "SNP", "SNP_N")
    }
  ))


mutsP <- muts_strain_TC_PF
mutsP <- mutsP[mutsP$Plot_gene == "other", ]

mutsP$events <- ifelse(
  is.na(mutsP$Type_IS),
  mutsP$Event,                    
  paste(mutsP$Event, mutsP$Type_IS, sep = "_")  
)


mutsP <- mutsP %>%
  mutate(Position = as.numeric(Position)) %>%
  mutate(plot = case_when(
    Secondary_gene == "ompR" ~ Position + 190,
    Secondary_gene == "rpoD" ~ Position + 914,
    Secondary_gene == "hldE" ~ Position + 2760,
    Secondary_gene == "transporter" ~ Position + 4258,
    Secondary_gene == "lpxC" ~ Position + 4652
  ))

mutsP$score <- NA
mutsP <- mutsP %>%
  mutate(score = case_when(
    STRAIN == "KPN01" ~ 1,
    STRAIN == "KPN10" ~ 51,
    STRAIN == "KPN13" ~ 101,
    STRAIN == "KPN16" ~ 151,
    STRAIN == "KPN08" ~ 201, 
    STRAIN == "Delta_KPN01" ~ 251,
    STRAIN == "Delta_KPN10" ~ 301,
    STRAIN == "Delta_KPN13" ~ 351,
    STRAIN == "Delta_KPN16" ~ 401, 
    STRAIN == "Delta_KPN08" ~ 451
  ))

positions <- mutsP %>%
  pull(plot)
events <- mutsP %>%
  pull(events)
strain <- mutsP %>%
  pull(STRAIN)
condition <- mutsP %>%
  pull(Condition)
is <- mutsP %>%
  pull(IS)
score <- mutsP %>%
  pull(score)

color <- rep(NA, length(events))
color <- ifelse(events == "NJ_IS1", "#810f7c", color) 
color <- ifelse (events == "NJ_IS5", "#dd55ffff", color) 
color <- ifelse (events == "NJ_Other", "#c6a8dfff", color) 
color <- ifelse(events == "SNP", "#aad0e2ff", color)
color <- ifelse(events == "SNP_N", "#69acecff", color)
color <- ifelse(events == "INDEL",  "#edf8fbff", color)
color <- ifelse(events =="SNP", "#b7c4c8ff", color)


side <- replace (condition, condition == "TC", "top")
side <- replace (side, side == "PF", "bottom")


mut_gr <- GRanges("chr1", IRanges(positions, width = 1))
mut_gr$color <- color
mut_gr$SNPsideID <-side
mut_gr$score <- score
mut_gr$label.parameter.rot <- 45


ompR_start<- 190
ompR_end<- 719
arrow <- GRanges("chr1", ranges = IRanges(start = ompR_start, end = ompR_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")

rpoD_start <- 914
rpoD_end <- 2755
arrow1 <- GRanges("chr1", ranges = IRanges(start = rpoD_start, end = rpoD_end), 
                 height = 0.05, shape = "arrow", fill = "lightgrey")


hldE_start <- 2760
hldE_end <- 4193
arrow2 <- GRanges("chr1", ranges = IRanges(start = hldE_start, end = hldE_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

t_start <- 4258
t_end <- 4647
arrow3 <- GRanges("chr1", ranges = IRanges(start = t_start, end = t_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

lpxC_start <- 4652
lpxC_end <- 5569
arrow4 <- GRanges("chr1", ranges = IRanges(start = lpxC_start, end = lpxC_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")


features <- c(arrow, arrow1, arrow2, arrow3, arrow4)

svg("misinggenes_fluctuation.svg")
lolliplot(mut_gr, features = features,
          ranges = GRanges("chr1", IRanges(0, 5600)),
          cex = 0.8 )
dev.off()
