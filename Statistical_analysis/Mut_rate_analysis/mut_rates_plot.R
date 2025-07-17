# Set path to folder with S. Table 2, by default, Downloads folder
setwd("~/Downloads")

# Load necessary libraries
library(xlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggpubr)
library(car)

# Set color palette
two_col_palette <- c("#BBBBBB", "#AA3377", "#286995")

# Load data, stored in the results sheet of the supplementary table 2
table_mut_rates <- read.xlsx("supplementary_table_2.xlsx", sheetName = "Results")

# Reorder table to plot as shown in the manuscript order
table_mut_rates$Strain <- factor(table_mut_rates$Strain, levels = c("KPN01", "KPN01p", "n1", "KPN10", "KPN10p", "n2", "KPN13", "KPN13p", "n3", "KPN16", "KPN16p", "n4", "KPN08", "KPN08p", "KPN08pΔΔIS1"))

# Plot Colistin results

COLPLOT <- table_mut_rates %>%
  filter(Antibiotic == "COL") %>%
  ggplot(aes(
    xmin = as.numeric(as.factor(Strain)) - 0.4,
    xmax = as.numeric(as.factor(Strain)) + 0.4,
    ymin = 10^-11,  # Log baseline for ymin
    ymax = Mutation.Rate,  # Bar height
    col = Genotype,
    fill = Genotype
  )) +
  geom_rect(
    aes(
      ymin = 10^-11,  # Correctly map ymin to the axis baseline
      ymax = Mutation.Rate  # Top of the bar
    ),
    position = position_dodge(0.8),
    alpha = 0.8,
    size = 1
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  geom_errorbar(
    aes(
      x = as.numeric(as.factor(Strain)),
      ymin = X84..IC.Lower.Limit,
      ymax = X84..IC.Upper.Limit
    ),
    width = 0.0,
    col = "black",
    position = position_dodge(0.7)
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(table_mut_rates_RIF$Strain)),  # Match strain levels
    labels = ifelse(
      grepl("^n\\d+$", levels(table_mut_rates_RIF$Strain)), 
      "",  # Skip dummy strain labels
      levels(table_mut_rates_RIF$Strain)
    ),
    expand = expansion(mult = c(0.05, 0.05))  # Add spacing
  ) +
  xlab("Strain") +
  ylab("Phenotypic COL resistant mutation rate") +
  scale_fill_manual(values = two_col_palette) +
  scale_color_manual(values = two_col_palette) +
  coord_cartesian(ylim = c(10^-11, 10^-5)) +
  theme_bw(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot Rifampicin results
RIFPLOT <- table_mut_rates %>%
  filter(Antibiotic == "RIF") %>%
  ggplot(aes(
    xmin = as.numeric(as.factor(Strain)) - 0.4,
    xmax = as.numeric(as.factor(Strain)) + 0.4,
    ymin = 10^-11,  # Log baseline for ymin
    ymax = Mutation.Rate,  # Bar height
    col = Genotype,
    fill = Genotype
  )) +
  geom_rect(
    aes(
      ymin = 10^-11,  # Correctly map ymin to the axis baseline
      ymax = Mutation.Rate  # Top of the bar
    ),
    position = position_dodge(0.8),
    alpha = 0.8,
    size = 1
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  geom_errorbar(
    aes(
      x = as.numeric(as.factor(Strain)),
      ymin = X84..IC.Lower.Limit,
      ymax = X84..IC.Upper.Limit
    ),
    width = 0.0,
    col = "black",
    position = position_dodge(0.7)
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(table_mut_rates_RIF$Strain)),  # Match strain levels
    labels = ifelse(
      grepl("^n\\d+$", levels(table_mut_rates_RIF$Strain)), 
      "",  # Skip dummy strain labels
      levels(table_mut_rates_RIF$Strain)
    ),
    expand = expansion(mult = c(0.05, 0.05))  # Add spacing
  ) +
  xlab("Strain") +
  ylab("Phenotypic RIF resistant mutation rate") +
  scale_fill_manual(values = two_col_palette) +
  scale_color_manual(values = two_col_palette) +
  coord_cartesian(ylim = c(10^-11, 10^-5)) +
  theme_bw(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggarrange(plotlist = list(COLPLOT, RIFPLOT), nrow = 2, common.legend =  TRUE, legend = "bottom")
