
### Analysis of MICs data from days 1, 2 and 3 at increasing concentrations
## This script takes as input the OD600 measurement at 22+/-2 hours
## of 96 well plates to then analyze the survivors of 2 K. pneumoniae
## strains both w/wo pOXA-48 plasmid under COL, CIP, CMP, KAN or RIF pressures


# Set wd to folder containing OD measurements from plate readers
setwd("~/experimental_IS1.2/CMIS")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(broom)

### Load data and add metadata (ab, ab concentration, strain and genotype)

# Loop through parent directory
parent_dir <- "~/experimental_IS1.2/CMIS"
subdirs <- list.dirs(parent_dir, recursive = FALSE)

MICs_data <- data.frame()

for (subdir in subdirs) {
  # Get the list of txt files in the current subdirectory
  txt_files <- list.files(subdir, pattern = "\\.txt$", full.names = TRUE)
  # Loop through each txt file in the current subdirectory
  for (filename in txt_files) {
    X <- read.table(filename, fill=T)[2:10,1:13] # Remove 1st empty row
    colnames(X) <- c(NA, seq(1,12))
    X <- X[-1,-1]
    # Reformat to longer format and Well column
    X <- pivot_longer(X, cols = 1:12, names_to = "Well", values_to = "OD")
    X <- as.data.frame(X)
    X$OD <- as.numeric(gsub(",","." , X$OD))
    X$Position <- c(rep("A", 12), rep("B", 12), rep("C", 12), rep("D", 12), rep("E", 12), rep("F", 12), rep("G", 12), rep("H", 12))
    X$Well <- gsub(" ", "", paste(X$Position, X$Well))
    X <- subset(X , select = -Position)
    # Filter and get Wells in which there's culture (checkerboard)
    header_info <- unlist(strsplit(gsub("\\.txt", "", basename(filename)), " ")[[1]], recursive = F) # Get Info from name
    X$Antibiotic <- header_info[1] # Antibiotic
    X$Concentration <- as.numeric(header_info[2]) # Ab concentration
    X$Strain <- header_info[3] # Strain name
    X$Day <- as.numeric(gsub("d", "",unlist(strsplit(dirname(filename), "/"), recursive = F)[7])) # Get day from folder (change to 7 if in Ubuntu)
    # Put Negative and Positive controls
    X$Positive <- as.numeric(X[72,2])
    if (header_info[1] == "Rifa") {
      X$Negative <- 0.045
    } else {
      X$Negative <- as.numeric(X[96,2])
    }
    # Correct OD removing negative control and remove neg values
    X$OD_corrected <- X$OD - X$Negative
    X[X<0] <- 0
    # Append the contents to the all_data data frame
    MICs_data <- rbind(MICs_data, X)
    rm(X) # Remove temporally stored data frame
  }
}
# Calculate positive control median of each day for each strain
mpos_rossett <- MICs_data %>%
  group_by(Strain, Day) %>%
  summarise(Median_Positive = median(Positive))
mpos_rossett <- as.data.frame(mpos_rossett)
mpos_rossett$S_D <- interaction(mpos_rossett$Strain, mpos_rossett$Day)
mpos_rossett <- mpos_rossett %>% select(Median_Positive, S_D)
MICs_data$S_D <- interaction(MICs_data$Strain, MICs_data$Day)
# Add to dataframe
MICs_data <- merge(MICs_data, mpos_rossett, by = "S_D")
# Remove temporal dataframe
rm(mpos_rossett)

# Now, calculate IC90 from positive control
MICs_data$IC90 <- 0.1*MICs_data$Median_Positive

# And get status column from IC90 (L <- alive cells; D <- dead cells)
MICs_data$Status <- with(MICs_data, ifelse(MICs_data$OD_corrected > MICs_data$IC90, "L", "D"))
# Keep out wells non cultured (checkerboard configuration starting in A1)
rosset_well <- c(paste0("A", seq(2,12, by =2)), paste0("B", seq(1,11,by = 2)),
                 paste0("C", seq(2,12, by =2)), paste0("D", seq(1,11,by = 2)),
                 paste0("E", seq(2,12, by =2)), paste0("F", seq(1,11,by = 2)),
                 paste0("G", seq(2,12, by =2)), paste0("H", seq(1,11,by = 2)))
rosset_well <- c(rosset_well, "F12", "G12")# Remove positive and negative controls
# Filter by well
MICs_data <- MICs_data %>% filter(!Well %in% rosset_well)

# Calculate median and mean OD600 values for each strain+antibiotic+concentration
# taking into account only L cells

rosset_median_L <- MICs_data %>%
  filter(Status == "L") %>%
  group_by(Strain, Antibiotic, Concentration) %>%
  summarise(Median_Live = median(OD_corrected), Mean_Live = mean(OD_corrected))
rosset_median_L <- as.data.frame(rosset_median_L)
rosset_median_L$S_D <- interaction(rosset_median_L$Strain, rosset_median_L$Antibiotic, rosset_median_L$Concentration)
rosset_median_L <- rosset_median_L %>% select(Median_Live, Mean_Live, S_D)
rosset_median_L$S_D <- as.vector(rosset_median_L$S_D)
# Include 0's for combinations with no L cells (i.e. K168 col 16, K253 cip 8, PF08 cip 128)
rosset_median_L <- rbind(rosset_median_L, c(0, 0, "K168.Colistina.16"))
rosset_median_L <- rbind(rosset_median_L, c(0, 0, "K253.Cipro.8"))
rosset_median_L <- rbind(rosset_median_L, c(0, 0, "PF08.Cipro.128"))

MICs_data$S_D <- interaction(MICs_data$Strain, MICs_data$Antibiotic, MICs_data$Concentration)
MICs_data <- merge(MICs_data, rosset_median_L, by = "S_D")
# And remove key interaction column
MICs_data <- subset(MICs_data, select = -S_D)
rm(rosset_median_L) # remove temporal data frame

MICs_data$Mean_Live <- as.numeric(MICs_data$Mean_Live)

########################## GLM models using OD data ########################## 

model_sum <- data.frame()

### Colistin model for K168/PF08

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Colistina"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Colistina"
logit_model$Strain <- "K168"
model_sum <- rbind(model_sum, logit_model)

# Generate prediction data for plotting K168_col model

new_data_model <- MICs_data %>% 
  filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Colistina")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K168_col <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "COL",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank()); logitmodelplot_K168_col


### Cipro model for K168/PF08

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Cipro"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Cipro"
logit_model$Strain <- "K168"
model_sum <- rbind(model_sum, logit_model)

# Generate prediction data for plotting K168_cipro model

new_data_model <- MICs_data %>% 
  filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Cipro")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K168_cipro <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "CIP",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank()); logitmodelplot_K168_cipro


### Cloranfenicol model for K168/PF08

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Cloranfenicol"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Cloranfenicol"
logit_model$Strain <- "K168"
model_sum <- rbind(model_sum, logit_model)

# Generate prediction data for plotting K168_cipro model

new_data_model <- MICs_data %>% 
  filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Cloranfenicol")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K168_cmp <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "CMP",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank()); logitmodelplot_K168_cmp


### Rifa model for K168/PF08

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Rifa"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Rifa"
logit_model$Strain <- "K168"
model_sum <- rbind(model_sum, logit_model)

# Generate prediction data for plotting K168_cipro model

new_data_model <- MICs_data %>% 
  filter(Strain == "K168" | Strain == "PF08", Antibiotic == "Rifa")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K168_rifa <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "RIF",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank()); logitmodelplot_K168_rifa

ggarrange(logitmodelplot_K168_col, logitmodelplot_K168_cipro,
          logitmodelplot_K168_cmp, logitmodelplot_K168_rifa, common.legend = TRUE, nrow = 1)

########### Models for K253 and PF16

### Colistin model for K253/PF16

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Colistina"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Colistina"
logit_model$Strain <- "K253"
model_sum <- rbind(model_sum, logit_model)

# Generate prediction data for plotting K253_col model

new_data_model <- MICs_data %>% 
  filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Colistina")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K253_col <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "COL",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank()); logitmodelplot_K253_col


### Cipro model for K253/PF16

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Cipro"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Cipro"
logit_model$Strain <- "K253"
model_sum <- rbind(model_sum, logit_model)

probit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                      filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Cipro"), family = binomial(link = "probit"))

summary(probit_model)
tidy(probit_model)

# Generate prediction data for plotting K253_cipro model

new_data_model <- MICs_data %>% 
  filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Cipro")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K253_cipro <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "CIP",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank()); logitmodelplot_K253_cipro


### Kanamicin model for K253/PF16

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Kana"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Kana"
logit_model$Strain <- "K253"
model_sum <- rbind(model_sum, logit_model)

probit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                      filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Kana"), family = binomial(link = "probit"))

summary(probit_model)
tidy(probit_model)

# Generate prediction data for plotting K253_cipro model

new_data_model <- MICs_data %>% 
  filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Kana")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

probit_pred <- predict(probit_model, new_data, type = "response", se.fit = TRUE)
new_data$probit_prob <- probit_pred$fit
new_data$probit_lower <- probit_pred$fit - 1.96 * probit_pred$se.fit
new_data$probit_upper <- probit_pred$fit + 1.96 * probit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K253_cmp <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "KAN",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank()); logitmodelplot_K253_cmp


### Rifa model for K253/PF16

logit_model <- glm(Status_binary ~ Concentration + Strain, data = MICs_data %>% 
                     filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Rifa"), family = binomial(link = "logit"))
summary(logit_model)
logit_model <- as.data.frame(tidy(logit_model))
logit_model$Antibiotic <- "Rifa"
logit_model$Strain <- "K253"
model_sum <- rbind(model_sum, logit_model)

# Generate prediction data for plotting K253_cipro model

new_data_model <- MICs_data %>% 
  filter(Strain == "K253" | Strain == "PF16", Antibiotic == "Rifa")

new_data <- expand.grid(
  Concentration = seq(min(new_data_model$Concentration), max(new_data_model$Concentration), length.out = 100),
  Strain = levels(factor(new_data_model$Strain))
)

logit_pred <- predict(logit_model, new_data, type = "response", se.fit = TRUE)
new_data$logit_prob <- logit_pred$fit
new_data$logit_lower <- logit_pred$fit - 1.96 * logit_pred$se.fit
new_data$logit_upper <- logit_pred$fit + 1.96 * logit_pred$se.fit

# Plotting observed data and model predictions
logitmodelplot_K253_rifa <- ggplot(new_data_model, aes(x = Concentration, y = logit_prob, color = Strain)) +
  #geom_point(alpha = 0.6, size = 2) +
  geom_line(data = new_data, aes(y = logit_prob), linetype = "solid", size = 1.5) +
  #geom_line(data = new_data, aes(y = probit_prob), linetype = "dashed", size = 1.5) +
  geom_ribbon(data = new_data, aes(ymin = logit_lower, ymax = logit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  #geom_ribbon(data = new_data, aes(ymin = probit_lower, ymax = probit_upper), linetype = 0, fill = c("lightsteelblue"), alpha = 0.2) +
  labs(
    title = "RIF",
    x = "Concentration",
    y = "Survival Probability",
    color = "Strain"
  ) +
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("#AA3377", "#BBBBBB"))+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank()); logitmodelplot_K253_rifa

ggarrange(logitmodelplot_K253_col, logitmodelplot_K253_cipro,
          logitmodelplot_K253_cmp, logitmodelplot_K253_rifa, common.legend = TRUE)

#write_xlsx(model_sum, "summary_GLM_stats.xlsx")

# Load data

MICs_data <- read.xlsx("MICs_data.xlsx", sheetIndex = 1)
MICs_data$Status_binary <- ifelse(MICs_data$Status == "L", 1, 0) # Add binary col 1 <- alive; 0 <- dead

result <- MICs_data %>%
  group_by(pf_name, Concentration, Antibiotic, Day) %>%
  summarise(L_count = sum(Status == "L", na.rm = TRUE), .groups = "drop")

result <- as.data.frame(result)

# Create the bar plot

# Add a dummy point for RIF
result <- result %>%
  bind_rows(data.frame(Concentration = 128,  # A new concentration value
                       pf_name = "KPN08",
                       Antibiotic = "RIF",
                       Day = 3,
                       L_count = 0  # Zero value for the dummy point
  ))

result <- result %>%
  bind_rows(data.frame(Concentration = 128,  # A new concentration value
                       pf_name = "KPN08p",
                       Antibiotic = "RIF",
                       Day = 3,
                       L_count = 0  # Zero value for the dummy point
  ))

result <- result %>%
  bind_rows(data.frame(Concentration = 128,  # A new concentration value
                       pf_name = "KPN16",
                       Antibiotic = "RIF",
                       Day = 3,
                       L_count = 0  # Zero value for the dummy point
  ))

result <- result %>%
  bind_rows(data.frame(Concentration = 128,  # A new concentration value
                       pf_name = "KPN16p",
                       Antibiotic = "RIF",
                       Day = 3,
                       L_count = 0  # Zero value for the dummy point
  ))

result

#write.xlsx(result, "results_MICs.xlsx")

# Plot the barplot shown in Fig. 2A

result <- read.xlsx("results_with_MICs.xlsx", sheetIndex = 1)

two_col_palette <-  c("#BBBBBB", "#AA3377")

KPN08_CMIs <- result %>%
  filter(pf_name == "KPN08" | pf_name == "KPN08p") %>%
  ggplot(aes(x = factor(MIC, levels = c("0.5xMIC", "MIC", "2xMIC", "4xMIC")), y = L_count, fill = factor(pf_name, levels = c("KPN08", "KPN08p")), col = factor(pf_name, levels = c("KPN08", "KPN08p")))) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.95), 
    alpha = 0.85, 
    size = 1.1) +
  theme_bw(base_size = 20) +
  ylim(0,50)+
  ylab("Surviving samples") +
  xlab("Concentration (μg/mL)")+
  facet_wrap(~ factor(Antibiotic, levels = c("COL", "CIP", "CMP", "RIF")), scales = "free_x", nrow = 1) +
  scale_fill_manual(values = two_col_palette) +
  scale_color_manual(values = two_col_palette) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        #axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank()); KPN08_CMIs

KPN16_CMIs <- result %>%
  filter(pf_name == "KPN16" | pf_name == "KPN16p") %>%
  ggplot(aes(x =factor(MIC, levels = c("0.5xMIC", "MIC", "2xMIC", "4xMIC")), y = L_count, fill = factor(pf_name, levels = c("KPN16", "KPN16p")), col = factor(pf_name, levels = c("KPN16", "KPN16p")))) +
  geom_bar(
    stat = "identity", 
    position = position_dodge(width = 0.95), 
    alpha = 0.85, 
    size = 1.1) +
  theme_bw(base_size = 20) +
  ylim(0,50)+
  ylab("Surviving samples") +
  xlab("Concentration (μg/mL)")+
  facet_wrap(~ factor(Antibiotic, levels = c("COL", "CIP", "KAN", "RIF")), scales = "free_x", nrow = 1) +
  scale_fill_manual(values = two_col_palette) +
  scale_color_manual(values = two_col_palette) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank()); KPN16_CMIs


