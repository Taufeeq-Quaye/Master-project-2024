library(EnhancedVolcano)
library(gplots)
library(circlize)
library(dplyr)
library(ggplot2)
library(readxl) 
library(tidyverse)

#Importing Data
expr1<- as.data.frame(read_xlsx("HMLE, D492, and MCF10a GRCh38.104 October 2021.xlsx"))
table(expr1$`Gene type`)

expr1 <- select(expr1,contains(c("Identifier","Gene type", "Mes vs. Ep, HMLE","Mes vs. Ep, D492", "Mes vs. Ep, MCF10A")))
expr1 <- as.data.frame(expr1)
rownames(expr1) <- expr1$Identifier

sigf <- as.data.frame(read.csv2("16 common LNCRNA HMLE, MCF10A, D492.csv", as.is = T, check.names = F ))
rownames(sigf) <- sigf$Identifier

#Removing NA values
sum(is.na(expr1))
expr1 <- na.omit(expr1)

HMLE<- expr1[,c(1,2,3,4)]
table(HMLE$`Gene type`)
D492 <- expr1[,c(1,2,5,6)]

MCF10A <- expr1[,c(1,2,7,8)]

#LOG 2 TRANSFORMATION OF FOLD CHANGE VALUES
HMLE$log_fold_change <- log2(abs(HMLE$`Mes vs. Ep, HMLE - Fold change`)) * sign(HMLE$`Mes vs. Ep, HMLE - Fold change`)
D492$log_fold_change <- log2(abs(D492$`Mes vs. Ep, D492 - Fold change`)) * sign(D492$`Mes vs. Ep, D492 - Fold change`)
MCF10A$log_fold_change <- log2(abs(MCF10A$`Mes vs. Ep, MCF10A - Fold change`)) * sign(MCF10A$`Mes vs. Ep, MCF10A - Fold change`)

#CHANGE COLUMN NAMES.
colnames(HMLE)[colnames(HMLE)== "Mes vs. Ep, HMLE - FDR p-value"] <- "FDR p-value"
colnames(HMLE)[colnames(HMLE)== "log_fold_change"] <- "Log2FC"

colnames(D492)[colnames(D492)== "Mes vs. Ep, D492 - FDR p-value"] <- "FDR p-value"
colnames(D492)[colnames(D492)== "log_fold_change"] <- "Log2FC"


colnames(MCF10A)[colnames(MCF10A)== "Mes vs. Ep, MCF10A - FDR p-value"] <- "FDR p-value"
colnames(MCF10A)[colnames(MCF10A)== "log_fold_change"] <- "Log2FC"

##

hmle_1 <- HMLE %>% filter(`FDR p-value` < 0.05)
x <- as.data.frame(table(hmle_1$`Gene type`))

mcf10a_1 <- MCF10A %>% filter(`FDR p-value` < 0.05)
x <- as.data.frame(table(mcf10a_1$`Gene type`))

d492_1 <- D492 %>% filter(`FDR p-value` < 0.05)
x <- as.data.frame(table(d492_1$`Gene type`))
##Volcano plot HMLE#############################################################################

#We will export the heatmap in a pdf for better view
pdf("HMLE-EP_MES_SOLid V104 volcano_plot.pdf",  width=11, height=7)

# Define your lncRNAs
lncs <- rownames(sigf)
custom_col <- rep("grey", length(rownames(HMLE)))
custom_shape <- rep(16, length(rownames(HMLE)))  # Default shape, e.g., circle


# Assign custom values for lncRNAs
custom_col[rownames(HMLE) %in% lncs] <- "red"  # Solid red for lncRNAs
#custom_shape[rownames(HMLE) %in% lncs] <- 17  # Custom shape (e.g., triangle) for lncRNAs
names(custom_col)[custom_col == 'red'] <- 'lncRNAs'
names(custom_col)[custom_col != 'red'] <- 'Others'
#names(custom_shape)[custom_shape == 17] <- 'lncRNAs'
#names(custom_shape)[custom_shape == 16] <- 'Others'

#If you only want to use different size for lncRNAs
custom_col2 <- ifelse(expr1$Log2FC > 1.5 & expr1$`FDR p-value` < 0.05, "red", ifelse(expr1$Log2FC < -1.5 & expr1$`FDR p-value` < 0.05, "royalblue", "black"))
names(custom_col2)[custom_col2 == 'red'] <- 'Upregulated'
names(custom_col2)[custom_col2 == 'royalblue'] <- 'Downregulated'
names(custom_col2)[custom_col2 == 'black'] <- 'Not Significant'

#EnhancedVolcano plot
ht <- EnhancedVolcano(HMLE,
                      lab = rownames(HMLE),
                      x = "Log2FC",
                      y = "FDR p-value",
                      title = 'Differential Expression of HMLE Mes vs Ep, Solid',
                      subtitle = expression('Log'[2]*'FC cutoff: 1.5, FDR p-value cutoff: 0.05'),
                      subtitleLabSize = 16,
                      selectLab = "ENSG00000230615",
                      xlab = expression('Log'[2]*'Fold Change'),
                      ylab=expression('-Log'[10]*'P'[adj]),
                      ylim = c(0,40),
                      caption = paste("Total number of significant genes(FDR p-value): ", sum(HMLE$`FDR p-value`< 0.05)),
                      captionLabSize = 20,
                      pointSize = c(ifelse(rownames(HMLE) %in% lncs, 5, 3)), #For changing size of the points
                      pCutoff = 0.05,
                      FCcutoff = 1.5,
                      #pointSize = 3.0,
                      labSize = 12.0,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = c(ifelse(rownames(HMLE) %in% lncs, 1.5, 0.1)), #for changing transparency of the points
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      typeConnectors = "closed",
                      endsConnectors = "first",
                      lengthConnectors = unit(0.03, "npc"),
                      border = 'full',
                      borderWidth = 1.0,
                      borderColour = 'black',
                      colCustom = custom_col,  # Highlight lncRNAs
                      axisLabSize = 22,
                      #shapeCustom = custom_shape,  # Different shape for lncRNAs
                      legendPosition = 'right',
                      legendLabSize = 25,
                      legendIconSize = 10)  # Legend position

# Print the plot
print(ht)

dev.off()
##volcano plot for MCF10A###############################################################################################################

#We will export the heatmap in a pdf for better view
pdf("MCF10A-EP_MES_SOLid V104 volcano_plot.pdf",  width=11, height=7)

# Define your lncRNAs
lncs <- rownames(sigf)
custom_col <- rep("grey", length(rownames(MCF10A)))
custom_shape <- rep(16, length(rownames(MCF10A)))  # Default shape, e.g., circle


# Assign custom values for lncRNAs
custom_col[rownames(MCF10A) %in% lncs] <- "red"  # Solid red for lncRNAs
#custom_shape[rownames(MCF10A) %in% lncs] <- 17  # Custom shape (e.g., triangle) for lncRNAs
names(custom_col)[custom_col == 'red'] <- 'lncRNAs'
names(custom_col)[custom_col != 'red'] <- 'Others'
#names(custom_shape)[custom_shape == 17] <- 'lncRNAs'
#names(custom_shape)[custom_shape == 16] <- 'Others'

#If you only want to use different size for lncRNAs
custom_col2 <- ifelse(expr1$Log2FC > 1.5 & expr1$`FDR p-value` < 0.05, "red", ifelse(expr1$Log2FC < -1.5 & expr1$`FDR p-value` < 0.05, "royalblue", "black"))
names(custom_col2)[custom_col2 == 'red'] <- 'Upregulated'
names(custom_col2)[custom_col2 == 'royalblue'] <- 'Downregulated'
names(custom_col2)[custom_col2 == 'black'] <- 'Not Significant'

# EnhancedVolcano plot
ht <- EnhancedVolcano(MCF10A,
                      lab = rownames(MCF10A),
                      x = "Log2FC",
                      y = "FDR p-value",
                      title = 'Differential Expression of MCF10A Mes vs Ep, Solid',
                      subtitle = expression('Log'[2]*'FC cutoff: 1.5, FDR p-value cutoff: 0.05'),
                      subtitleLabSize = 16,
                      selectLab = "ENSG00000230615",
                      xlab = expression('Log'[2]*'Fold Change'),
                      ylab=expression('-Log'[10]*'P'[adj]),
                      ylim = c(0,40),
                      caption = paste("Total number of significant genes(FDR p-value): ", sum(MCF10A$`FDR p-value`< 0.05)),
                      captionLabSize = 20,
                      pointSize = c(ifelse(rownames(MCF10A) %in% lncs, 5, 3)), #For changing size of the points
                      pCutoff = 0.05,
                      FCcutoff = 1.5,
                      #pointSize = 3.0,
                      labSize = 12.0,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = c(ifelse(rownames(MCF10A) %in% lncs, 1.5, 0.1)), #for changing transparency of the points
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      typeConnectors = "closed",
                      endsConnectors = "first",
                      lengthConnectors = unit(0.03, "npc"),
                      border = 'full',
                      borderWidth = 1.0,
                      borderColour = 'black',
                      colCustom = custom_col,  # Highlight lncRNAs
                      axisLabSize = 22,
                      #shapeCustom = custom_shape,  # Different shape for lncRNAs
                      legendPosition = 'right',
                      legendLabSize = 25,
                      legendIconSize = 10)  # Legend position

# Print the plot
print(ht)

dev.off()

##volcano plot for d492#############################################################################################################################################
##############################################################################################################################################################################

#We will export the heatmap in a pdf for better view
pdf("D492-EP_MES_SOLid V104 volcano_plot.pdf",  width=11, height=7)

# Define your lncRNAs
lncs <- rownames(sigf)
custom_col <- rep("grey", length(rownames(D492)))
custom_shape <- rep(16, length(rownames(D492)))  # Default shape, e.g., circle


# Assign custom values for lncRNAs
custom_col[rownames(D492) %in% lncs] <- "red"  # Solid red for lncRNAs
#custom_shape[rownames(D492) %in% lncs] <- 17  # Custom shape (e.g., triangle) for lncRNAs
names(custom_col)[custom_col == 'red'] <- 'lncRNAs'
names(custom_col)[custom_col != 'red'] <- 'Others'
#names(custom_shape)[custom_shape == 17] <- 'lncRNAs'
#names(custom_shape)[custom_shape == 16] <- 'Others'

#If you only want to use different size for lncRNAs
custom_col2 <- ifelse(expr1$Log2FC > 1.5 & expr1$`FDR p-value` < 0.05, "red", ifelse(expr1$Log2FC < -1.5 & expr1$`FDR p-value` < 0.05, "royalblue", "black"))
names(custom_col2)[custom_col2 == 'red'] <- 'Upregulated'
names(custom_col2)[custom_col2 == 'royalblue'] <- 'Downregulated'
names(custom_col2)[custom_col2 == 'black'] <- 'Not Significant'

# EnhancedVolcano plot

ht <- EnhancedVolcano(D492,
                      lab = rownames(D492),
                      x = "Log2FC",
                      y = "FDR p-value",
                      title = 'Differential Expression of D492 Mes vs Ep, Solid',
                      subtitle = expression('Log'[2]*'FC cutoff: 1.5, FDR p-value cutoff: 0.05'),
                      subtitleLabSize = 16,
                      selectLab = "ENSG00000230615",
                      xlab = expression('Log'[2]*'Fold Change'),
                      ylab=expression('-Log'[10]*'P'[adj]),
                      ylim = c(0,40),
                      caption = paste("Total number of significant genes(FDR p-value): ", sum(D492$`FDR p-value`< 0.05)),
                      captionLabSize = 20,
                      pointSize = c(ifelse(rownames(D492) %in% lncs, 5, 3)), #For changing size of the points
                      pCutoff = 0.05,
                      FCcutoff = 1.5,
                      #pointSize = 3.0,
                      labSize = 12.0,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      colAlpha = c(ifelse(rownames(D492) %in% lncs, 1.5, 0.1)), #for changing transparency of the points
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      typeConnectors = "closed",
                      endsConnectors = "first",
                      lengthConnectors = unit(0.03, "npc"),
                      border = 'full',
                      borderWidth = 1.0,
                      borderColour = 'black',
                      colCustom = custom_col,  # Highlight lncRNAs
                      axisLabSize = 22,
                      #shapeCustom = custom_shape,  # Different shape for lncRNAs
                      legendPosition = 'right',
                      legendLabSize = 25,
                      legendIconSize = 10)  # Legend position

# Print the plot
print(ht)

dev.off()

#c("ENSG00000230615"," ENSG00000130600", " ENSG00000225746", "ENSG00000163597")
