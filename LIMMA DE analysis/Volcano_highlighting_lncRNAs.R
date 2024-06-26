library(EnhancedVolcano)
library(gplots)
library(circlize)
library(dplyr)
library(ggplot2)
library(readxl) 
library(tidyverse)

#Importing Data
expr1 <- as.data.frame(read_xlsx("HMLE EP-MES v104.xlsx"))
sigf <- as.data.frame(read.csv2("16 common LNCRNA HMLE, MCF10A, D492.csv", as.is = T, check.names = F ))
rownames(sigf) <- sigf$Identifier


expr1 <- expr1[,c(4,7,9)]
expr1 <- as.data.frame(expr1)
rownames(expr1) <- expr1$Identifier

#Removing NA values

sum(is.na(expr1))
expr1 <- na.omit(expr1)
#CHANGE COLUMN NAMES.

colnames(expr1)[colnames(expr1)== "Mes vs. Ep - FDR p-value"] <- "FDR p-value"
colnames(expr1)[colnames(expr1)== "Mes vs. Ep - Log fold change"] <- "Log2FC"


#We will export the heatmap in a pdf for better view
pdf("HMLE-EP_MES V104 volcano_plot.pdf",  width=11, height=7)

# Define your lncRNAs
lncs <- rownames(sigf)

#If you want to show lncRNAs using both different shapes and colors
# Assign custom values for lncRNAs
custom_col <- rep("grey", length(rownames(expr1)))
custom_shape <- rep(16, length(rownames(expr1)))  # Default shape (e.g., circle)

# Assign custom values for lncRNAs
custom_col[rownames(expr1) %in% lncs] <- "red"  # Solid red for lncRNAs
names(custom_col)[custom_col == 'red'] <- 'lncRNAs'
names(custom_col)[custom_col != 'red'] <- 'Genes'

#If you only want to use different size for lncRNAs
custom_col2 <- ifelse(expr1$Log2FC > 1.5 & expr1$`FDR p-value` < 0.05, "red", ifelse(expr1$Log2FC < -1.5 & expr1$`FDR p-value` < 0.05, "royalblue", "black"))
names(custom_col2)[custom_col2 == 'red'] <- 'Upregulated'
names(custom_col2)[custom_col2 == 'royalblue'] <- 'Downregulated'
names(custom_col2)[custom_col2 == 'black'] <- 'Not Significant'

# EnhancedVolcano plot
ht <- EnhancedVolcano(expr1,
                      lab = rownames(expr1),
                      x = 'Log2FC',
                      y = 'FDR p-value',
                      title = 'Differential Expression of HMLE Mes vs Ep, Illumina',
                      subtitle = expression('Log'[2]*'FC cutoff: 1.5, FDR p-value cutoff: 0.05'),
                      subtitleLabSize = 16,
                      selectLab = "ENSG00000230615",
                      xlab = expression('Log'[2]*'Fold Change'),
                      ylab=expression('-Log'[10]*'P'[adj]),
                      ylim = c(0,40),
                      caption = paste("Total number of significant genes(FDR p-value): ", sum(expr1$`FDR p-value` < 0.05)),
                      captionLabSize = 20,
                      pointSize = c(ifelse(rownames(expr1) %in% lncs, 5, 3)), #For changing size of the points
                      pCutoff = 0.05,
                      FCcutoff = 1.5,
                      #pointSize = 3.0,
                      labSize = 12.0,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = c(ifelse(rownames(expr1) %in% lncs, 1.5, 0.1)), #for changing transparency of the points
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      border = 'full',
                      typeConnectors = "closed",
                      endsConnectors = "first",
                      lengthConnectors = unit(0.03, "npc"),
                      borderWidth = 1.0,
                      borderColour = 'black',
                      colCustom = custom_col,
                      axisLabSize = 22,
                      #shapeCustom = custom_shape,  # Different shape for lncRNAs
                      legendPosition = 'right',
                      legendLabSize = 25,
                      legendIconSize = 10)  # Legend position


# Print the plot
print(ht)

dev.off()

housekeeping_genes <- c("ENSG00000111640","ENSG00000075624","ENSG00000166710","ENSG00000162881","ENSG00000165704","ENSG00000112592","ENSG00000102144","ENSG00000196262","ENSG00000169919","ENSG00000164924")
  SDHA = "ENSG00000073578",
  TFRC = "ENSG00000072274",
  RPL13A = "ENSG00000147403",
  HMBS = "ENSG00000244734",
  TUBB = "ENSG00000196230",
  UBC = "ENSG00000150991",
  ATP5F1 = "ENSG00000261731",
  RPS18 = "ENSG00000288647",
  EIF4A2 = "ENSG00000161960",
  PSMC4 = "ENSG00000105974"
)

