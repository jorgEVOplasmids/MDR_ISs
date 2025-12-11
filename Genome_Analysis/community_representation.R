############################## Community data ##################################
#We start representing KPN08

muts_strain_TC_PF <- read_xlsx("community.xlsx", sheet = 1 )
muts_strain_TC_PF <- muts_strain_TC_PF[!is.na(muts_strain_TC_PF$Plot_gene), ]
muts_strain_TC_PF <- muts_strain_TC_PF[muts_strain_TC_PF$STRAIN == "KPN08", ]

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

# For better visualization, we only plot the most repeated set of genes. Therefore, we
# exclude the "other" genes. Also, due to the small size of mgrB (144 bp), it will be plotted
# independently,

mutsP <- muts_strain_TC_PF
mutsP <- mutsP[mutsP$Plot_gene != "other", ]

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
    Plot_gene == "mgrB" ~ Position + 31 ,
    Plot_gene == "crrB" ~ Position + 180,
    Plot_gene == "phoQ" ~ Position + 1245,
    Plot_gene == "pmrA" ~ Position + 2720,
    Plot_gene == "pmrB" ~ Position + 3388
  ))

# Setting different highs or scores to differentiate between strains.
mutsP$score <- NA
mutsP <- mutsP %>%
  mutate(score = case_when(
    STRAIN == "KPN08" ~ 11
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
color <- ifelse(events =="DEL", "#b7c4c8ff", color)
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
mgrb_start <- 2841
mgrb_end<- 2985
arrow <- GRanges("chr1", ranges = IRanges(start = mgrb_start, end = mgrb_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")
crrb_start <- 3755
crrb_end <- 4816
arrow1 <- GRanges("chr1", ranges = IRanges(start = crrb_start, end = crrb_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

phoQ_start <- 5491
phoQ_end <- 6957
arrow2 <- GRanges("chr1", ranges = IRanges(start = phoQ_start, end = phoQ_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

pmra_start <- 6965
pmra_End <- 7636
arrow3 <- GRanges("chr1", ranges = IRanges(start = pmra_start, end = pmra_End), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

pmrb_start <- 7700
pmrb_end <- 8797
arrow4 <- GRanges("chr1", ranges = IRanges(start = pmrb_start, end = pmrb_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

# Combine all the arrows as features and plotting
features <- c(arrow, arrow1, arrow2, arrow3, arrow4)

svg("KPN08_Community.svg")
lolliplot(mut_gr, features = features,
          ranges = GRanges("chr1", IRanges(0, 9000)), cex= 1)

dev.off()



## Then we plot Citrobacter freundii (CF13)
muts_strain_TC_PF <- read_xlsx("community.xlsx", sheet = 1 )
muts_strain_TC_PF <- muts_strain_TC_PF[!is.na(muts_strain_TC_PF$Plot_gene), ]
muts_strain_TC_PF <- muts_strain_TC_PF[muts_strain_TC_PF$STRAIN == "CF13", ]

#Detect whether the SNPs are non synonymous + non sense, synonymous or intergenic by 
#extracting the aminoacid change 
mutsP <- muts_strain_TC_PF %>%
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
mutsP$events <- ifelse(
  is.na(mutsP$IS),
  mutsP$Event,                    
  paste(mutsP$Event, mutsP$IS, sep = "_")  
)


# Calculation of the position of the genes for representation according to 
# the reference position
mutsP <- mutsP %>%
  mutate(plot = case_when(
    Plot_gene == "glycosyl" ~ Position + 1058 ,
    Plot_gene == "pmrA" ~ Position + 2850,
    Plot_gene == "pmrB" ~ Position + 3519
  ))

# Setting different highs or scores to differentiate between strains.
mutsP$score <- NA
mutsP <- mutsP %>%
  mutate(score = case_when(
    STRAIN == "CF13" ~ 11
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
color <- ifelse(events =="DEL", "#b7c4c8ff", color)
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
rfaQ_start <- 1
rfaQ_end<- 1056
arrow1 <- GRanges("chr1", ranges = IRanges(start = rfaQ_start, end = rfaQ_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")
glyc1_start <- 1058
glyc1_end <- 1840
arrow2 <- GRanges("chr1", ranges = IRanges(start = glyc1_start, end = glyc1_end), 
                  height = 0.05, shape = "arrow", fill = "lightgrey")

IS1_start <- 1890
IS1_end <- 2595
arrow3 <- GRanges("chr1", ranges = IRanges(start = IS1_start, end = IS1_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

glyc2_start <- 2611
glyc2_End <- 2847
arrow4 <- GRanges("chr1", ranges = IRanges(start = glyc2_start, end = glyc2_End), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")
pmrA_start <- 2850
pmrA_end <- 3518
arrow5 <- GRanges("chr1", ranges = IRanges(start = pmrA_start, end = pmrA_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")

pmrb_start <- 3519
pmrb_end <- 4601
arrow6 <- GRanges("chr1", ranges = IRanges(start = pmrb_start, end = pmrb_end), 
                     height = 0.05, shape = "arrow", fill = "lightgrey")
# Combine all the arrows as features and plotting
features <- c(arrow1, arrow2, arrow3, arrow4, arrow5, arrow6)

svg("community_CF13_glycosyl.svg")
lolliplot(mut_gr, features = features,
          ranges = GRanges("chr1", IRanges(0, 4601)), cex= 1)

dev.off()



#We represent all the mutations in a piecharts
muts_strain_TC_PF <- read_xlsx("community.xlsx", sheet = 1 )
muts_strain_TC_PF <- muts_strain_TC_PF[!is.na(muts_strain_TC_PF$Plot_gene), ]
muts_counter <- muts_strain_TC_PF[muts_strain_TC_PF$STRAIN == "CF13", ]


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
muts_TC <- muts_TC[muts_TC$STRAIN != "Delta_IS", ]


counter_PF <- muts_PF %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_PF$events <- factor(counter_PF$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", 
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
                               "#a9a9a9ff"
  ))

ggsave("counts_CF13_community_PF.svg", plot = PF)



counter_TC <- muts_TC %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_TC$events <- factor(counter_TC$events, levels = c("DEL", "INDEL", 
                                                          "SNP_I", "SNP_N", 
                                                          "NJ_N","NJ_Y"))


TC <- ggplot(counter_TC, aes(x = "", y = Count, fill = events)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + 
  coord_polar("y", start = 0) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 5, color = "black", fontface = "bold") + 
  #labs(title = "Mutational Events per Strains in TC") +
  theme_void() + 
  scale_fill_manual(values = c("INDEL" = "#edf8fbff",  
                               "NJ_N" = "#c6a8dfff",    
                               "NJ_Y" = "#810f7c",     
                               "SNP_I" = "#aad0e2ff",  
                               "SNP_N" = "#69acecff", 
                               "DEL" = "#575757ff",
                               "#a9a9a9ff"
  ))

ggsave("counts_CF13_community_TC.svg", plot = TC)



# Lastly the piecharts for KPN08
muts_strain_TC_PF <- read_xlsx("community.xlsx", sheet = 1 )
muts_strain_TC_PF <- muts_strain_TC_PF[!is.na(muts_strain_TC_PF$Plot_gene), ]
muts_counter <- muts_strain_TC_PF[muts_strain_TC_PF$STRAIN == "KPN08", ]



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
muts_TC <- muts_TC[muts_TC$STRAIN != "Delta_IS", ]


counter_PF <- muts_PF %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_PF$events <- factor(counter_PF$events, levels = c("DEL", "INDEL", "SNP_I", 
                                                          "SNP_N", "NJ_N", 
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
                               "#a9a9a9ff"
  ))

ggsave("counts_KPN08_community_PF.svg", plot = PF)



counter_TC <- muts_TC %>%
  group_by(events) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = Count / sum(Count) * 100)
counter_TC$events <- factor(counter_TC$events, levels = c("DEL", "INDEL", 
                                                          "SNP_I", "SNP_N", 
                                                          "NJ_N","NJ_Y"))


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
                               "#a9a9a9ff"
  ))

ggsave("counts_KPN08_community_TC.svg", plot = TC)
