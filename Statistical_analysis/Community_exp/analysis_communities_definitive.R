
setwd("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/experimental_IS1.2/community/experiment_01052025")

library(xlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Import sheet with main statistics

communities_data <- read.xlsx("community experiment 01_05_25 segundo plaqueo_w_transconjugants.xlsx", sheetIndex = 9)

# Statistically compare the number of TCs in the whole population and the ColR population

table_stats <- read.xlsx("community experiment 01_05_25 segundo plaqueo_w_transconjugants.xlsx", sheetIndex = 10)

contingency_table_CF <- data.frame(TC = c(146, 103), PF = c(19,97))
rownames(contingency_table_CF) <- c("ColR", "Col_S")
contingency_table_EC <- data.frame(TC = c(6,0), PF = c(240,200))
rownames(contingency_table_EC) <- c("ColR", "Col_S")
contingency_table_KPN <- data.frame(TC = c(149,79), PF = c(77,121))
rownames(contingency_table_KPN) <- c("ColR", "Col_S")

# Statistic results shown in main Fig 4D

chisq.test(contingency_table_CF)
chisq.test(contingency_table_EC)
chisq.test(contingency_table_KPN)

# Now, plot the results
# First, the community composition as a stacked barplot for either pOXA+/pOXA-  (Main Fig. 4B without TC shading)

custom_palette <- c("#E4572E", "#F3A712", "#29335C")

communities_data %>%
  ggplot(aes(x = factor(Tube, levels = c("B2", "B3", "B4", "B5", "B6", "B7", "B10", "B11", "B12", "B13", "A2", "A3", "A4", "A5", "A6", "A7", "A10", "A11", "A12", "A13")),
             y = pop_size, fill = Strain, col = Strain)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Population (%)") +
  xlab("Tube") +
  scale_fill_manual(values = custom_palette)+
  scale_color_manual(values = custom_palette) +
  theme_bw(base_size = 12) +
  facet_wrap(~factor(pOXA, levels = c("YES", "NO")), nrow = 2, scales = "free_x") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        #aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = element_blank())

# Now plot the mutants/cell corrected by transconjugants frequency (which will be later overlayed to plot TCs in Fig. 4B)

communities_data %>%
  filter(pOXA == "YES") %>%
  ggplot(aes(x = factor(Tube, levels = c("B2", "B3", "B4", "B5", "B6", "B7", "B10", "B11", "B12", "B13", "A2", "A3", "A4", "A5", "A6", "A7", "A10", "A11", "A12", "A13")),
             y = freq_TC_ColR, fill = Strain, col = Strain)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Mutants/cell") +
  xlab("Tube") +
  scale_fill_manual(values = custom_palette)+
  scale_color_manual(values = custom_palette) +
  theme_bw(base_size = 12) +
  #scale_y_log10() +
  facet_wrap(~factor(pOXA, levels = c("YES", "NO")), nrow = 2, scales = "free_x") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        #aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = element_blank())

# Plot the mutants per treatment

two_col_palette <- c("#BBBBBB", "#AA3377")

# Perform statistics

communities_data_CF <- communities_data %>%
  filter(Strain == "CF13")

communities_data_EC <- communities_data %>%
  filter(Strain == "EC17")

communities_data_KPN <- communities_data %>%
  filter(Strain == "KPN08")

t.test(as.numeric(communities_data_CF$col_mut_prop_to_pop_size) ~ communities_data_CF$pOXA, alternative = "two.sided")
t.test(as.numeric(communities_data_EC$col_mut_prop_to_pop_size) ~ communities_data_EC$pOXA, alternative = "two.sided")
t.test(as.numeric(communities_data_KPN$col_mut_prop_to_pop_size) ~ communities_data_KPN$pOXA, alternative = "two.sided")

# First, summarize the data
summary_df <- communities_data %>%
  #filter(Strain != "EC17") %>%
  group_by(Strain, pOXA) %>%
  summarise(
    mean = mean(col_mut_prop_to_pop_size, na.rm = TRUE),
    se = sd(col_mut_prop_to_pop_size, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# And plot the COL resistant mutants per population cell (main fig. 4C)
communities_data %>%
  #filter(Strain != "EC17") %>%
  ggplot(aes(x = factor(Strain, levels = c("EC17", "CF13", "KPN08")), y = col_mut_prop_to_pop_size, fill = pOXA, color = pOXA)) +
  geom_bar(data = summary_df,
           aes(y = mean),
           stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.6,
           size = 1,
           alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.15),
              size = 3, shape = 21, stroke = 0.4, alpha = 0.9) +
  geom_errorbar(data = summary_df,
                aes(y = mean, ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.7),
                width = 0, col = "black", size = 1) +
  ylab("Colistin resistant mutants / cell") +
  xlab("Strain") +
  scale_y_continuous(
    limits = c(0, 2.5e-06),
    breaks = c(0, 5e-07, 1e-06, 1.5e-06, 2e-06, 2.5e-06),
    labels = scales::scientific
  )+
  theme_bw(base_size = 20) +
  scale_fill_manual(values = two_col_palette) +
  scale_color_manual(values = two_col_palette) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank())

# And compare frequency of COL mutants per treatment (stat in main Fig. 4C)
summary_df <- communities_data %>%
  group_by(pOXA) %>%
  summarise(
    mean = mean(col_mut_prop_to_pop_size, na.rm = TRUE),
    se = sd(col_mut_prop_to_pop_size, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

aov.community <- aov(col_mut ~ pOXA, data = communities_data)
summary(aov.community)


#### Include transconjugants data

setwd("/home/jorge/Documents/important_docs/paper_IS1.2/responses_reviewers/independent_conjugations")

conjugation_data <- read.xlsx("conjs_limpio.xlsx", sheetIndex = 3)

# Now plot independent conjugation stuff keeping colour palette

custom_palette <- c("#E4572E", "#F3A712", "#29335C")

### Plot

summary_df <- conjugation_data %>%
  group_by(strain, time) %>%
  summarise(
    mean = median(log(conj_freq), na.rm = TRUE),
    se = sd(log(conj_freq), na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Finally, plot independent conjugations results (frequency of transconjugants)

library(cfbplotR)
library(ggplot2)

conjugation_data %>%
  ggplot(aes(x = as.factor(time), 
             y = log(conj_freq), 
             col = strain, fill = strain)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.15),
              size = 3, shape = 21, stroke = 0.9, alpha = 0.4) +
  stat_summary(
    fun = median, 
    geom = "crossbar", 
    aes(group = interaction(strain)), 
    linewidth = 0.8,
    width = 0.4, alpha = 0.5
  ) +
  geom_hline(yintercept = log(0.000001), col = "darkgrey", linetype = "dashed")+
  scale_fill_manual(values = custom_palette) +
  scale_color_manual(values = custom_palette) +
  # --- y-axis as original-scale scientific notation ---
  scale_y_continuous(
    breaks = log(c(1e0, 1e-2, 1e-4, 1e-6, 1e-8)),
    labels = c("1e+00", "1e-02", "1e-04", "1e-06", "1e-08")
  ) +
  theme_bw(base_size = 16) +
  facet_wrap(~strain) +
  ylab("Conjugation frequency (TC/R)") +
  xlab("Conjugation time (hours)") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank()
  )
