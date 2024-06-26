library(decoupleR)
library(tidyverse)
library(readxl)
library(devtools)
library(circlize)
library(ggpubr)
library(ggplot2)
#library(limma)
#library(edgeR)
library(OmnipathR)


# Read data from Excel and CSV files
KD <- as.data.frame(read_xlsx("RP5 KO.xlsx"))
KD <- KD[, c(1,4)]
limma <- read.csv2("Scramble VS Gapmer 2 and 5 RP5 KD.csv", as.is = TRUE, check.names = FALSE, row.names = 1)
limma <- merge(KD, limma, by="Identifier")
colnames(limma)

# Extract relevant columns from limma data
l <- as_tibble(limma) %>%
  select(Name, Scramble_vs_Gapmer_2_5_LogFC, t, P.Value) %>%
  filter(!is.na(t)) %>%
  column_to_rownames(var = "Name") %>%
  as.matrix()

# Get a protein-protein interaction network
net <- get_collectri(organism = 'human', split_complexes = FALSE)

# Run weighted mean of contrasts analysis
contrast_acts <- run_wmean(mat = l[, 't', drop = FALSE], net = net, .source = 'source', .target = 'target',
                           .mor = 'mor', times = 100, minsize = 5)


#The aim of this code is to perform a weighted mean of contrasts analysis on a specified column ('t') of a data frame l, considering a network structure defined in the net variable. 
#This type of analysis is commonly used in systems biology to infer the activity of biological entities (e.g., transcription factors) based on their relationships with other entities in a network. 
#The weighted mean takes into account the strength of connections between nodes in the network.

# Filter and rank top transcription factors
f_contrast_acts <- contrast_acts %>%
  filter(statistic == 'norm_wmean') %>%
  filter(p_value < 0.05)%>%
  mutate(rnk = NA)

# Filter top TFs in both signs
#In summary, the code is creating a ranking column ('rnk') in the data frame f_contrast_acts. 
#The ranking is based on the 'score' column, and it ensures that higher absolute scores receive higher ranks, 
#regardless of their sign. The ranks for positive scores are assigned as the negative of their actual values to ensure the correct ordering. 
#This type of ranking is often used in analyses where both positive and negative values need to be considered in a unified manner.

msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

n_tfs <- 25
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)

f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

# Plot the top 25 TF activity inference
plot_top25 <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) +
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "#2E5B88", high = "indianred", mid = "whitesmoke", midpoint = 0) +
  theme_pubr() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("Transcription Factors (TF's)") + labs(title = "Top 25 Significant TF's Activity Inference", subtitle = "HMLE-mes Gapmer 2 and 5 RP5 KD vs Control")

# Save the plot to PDF
pdf("TF Filter analyses for Gapmer 2&5 Rp5 KD.pdf", width = 8, height = 6)
print(plot_top25)
dev.off()

# Select specific TFs for EMT analysis
tf <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")
tf_f <- contrast_acts %>%
  filter(statistic == 'norm_wmean') %>%
  filter(source %in% tf) %>%
  arrange(score)
tf_f$Significant <- tf_f$p_value < 0.05

# Plot EMT TF activity inference
plot_emt <- ggplot(tf_f, aes(x = reorder(source, score), y = score, fill = Significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#2E5B88", "indianred")) +
  theme_pubr() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("EMT Transcription Factors") + labs(title = "EMT TF Activity Inference", subtitle = "HMLE-mes Gapmer 2 and 5 RP5 KD vs Control")

# Save the plot to PDF
pdf("EMT TF analyses for Gapmer 2 and 5 Rp5 KD.pdf", width = 7, height = 5)
print(plot_emt)
dev.off()
