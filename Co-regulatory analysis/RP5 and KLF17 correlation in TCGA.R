library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(readxl)
library(ggpubr)


#expression #DATA PROCESSING
expr <- read.csv2("TCGA_log.csv", sep = ";", as.is = T, check.names = F)
clinical <- read.csv2("TCGA_BRCA_Clinical.csv", as.is = T, check.names = F, row.names = 1)

x <- expr %>% filter(gene %in%  "KLF17")
rownames(x) <- x$gene
x <- x[,-c(1,2)]
x<- 2^(x)-1
x$gene <- rownames(x)
z <- expr %>% filter(id %in%  "ENSG00000230615")
rownames(z) <- z$id
z<- z[,-c(1,2)]
z<- 2^(z)-1
z$gene <- rownames(z)
KLF17 <- pivot_longer(x, cols= 1:1041, names_to = "Patient_ID", values_to = "KLF17")

RP5 <- pivot_longer(z, cols= 1:1041, names_to = "Patient_ID", values_to = "ENSG00000230615")

y <-merge(RP5, KLF17, by="Patient_ID")
y <- y[, c(1,2,14)]

y <- y[, c(1,3,5)]

###create a scatter correlation graph

ggscatter(y, x = "KLF17", y = "ENSG00000230615",
          title = "RP5 and KLF17 correlation in TCGA BRCA (n = 1,001)",
          color = "black", shape = 21, size = 2,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = "KLF17", ylab = "ENSG00000230615")+ 
 theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=20, ), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 18,face = "bold"), axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 20))

#open pdf
pdf("RP5 and klf17 Correlation in TCGA.pdf",  width=10,height=7)

dev.off()
#########################################################################################################
#####SLC6A9
x <- expr %>% filter(gene %in%  "SLC6A9")
SLC6A9 <- pivot_longer(x, cols= 3:1043, names_to = "Patient_ID", values_to = "SLC6A9")

y <-merge(clinical, SLC6A9, by="Patient_ID")
y <- y[, c(1,2,14)]

###create a scatter correlation graph

ggscatter(y, x = "SLC6A9", y = "ENSG00000230615",
          title = "RP5 and SLC6A9 correlation in TCGA BRCA (n = 1,001)",
          color = "black", shape = 21, size = 2,fill="#D73027",alpha=0.4, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#D73027", fill = "#D73027"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE,
          cor.coef.size = 5,# Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = "SLC6A9", ylab = "ENSG00000230615")+ 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face ="italic"), #Italic if it is a gene. 
        axis.text.x = element_text(size=20, ), axis.ticks.x=element_blank(), 
        axis.title.x = element_text(size = 18,face = "bold"), axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 20))

#open pdf
pdf("RP5 and SLC6A9 Correlation in TCGA.pdf",  width=9,height=7)

dev.off()

#####################################################################################################################
#####################################################################################################################
