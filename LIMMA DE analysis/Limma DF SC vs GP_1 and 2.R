library(dplyr)
library(ggpubr)
library(limma)
library(magrittr)
library(edgeR)
library(readxl)
#library(DESeq2)

#Importing dataset
KD <- as.data.frame(read_xlsx("RP5 KO.xlsx"))

rownames(KD) <- KD$Identifier
df <- dplyr::select(KD, contains("Total"))

#Get the metadata
# Create a sample info dataframe.
info <- data.frame(
  Sample = colnames(df),
  State = c("Treated", "Treated", "Treated","Treated", "Treated", "Treated", "Untreated", "Untreated","Untreated"),
  Replicate = c("Rep 1", "Rep 2", "Rep 3","Rep 4", "Rep 5", "Rep 6", "Rep 1", "Rep 2","Rep 3"))
 
rownames(info) <- info$Sample
info <- subset(info, (rownames(info) %in% colnames(df)))
info <- as.data.frame(info)

#now we can rename the samples for better visualization
rownames(info) <- gsub(" - Total counts", "", rownames(info))
colnames(df) <- gsub(" - Total counts", "", colnames(df))


#Check if the column of counts corresponds to the metadata rows
all(colnames(df) == rownames(info))
group <- as.factor(info$State)

# Create a Differential Gene Expression List (DGEList) object in the limma pipeline
# - This is a crucial step in the limma pipeline for RNA-seq data analysis.
# - The DGEList object contains the count data, information about experimental groups,
#   and provides a platform for downstream analyses, such as normalization, differential expression analysis, etc.


d0 <- DGEList(df, group = group) #Converting to DGE list for edgeR

d0$genes <- data.frame(Symbol=rownames(d0))

cpm <- cpm(d0)
lcpm <- cpm(d0, log=TRUE)

keep.exprs <- filterByExpr(d0, group=group) #By default it filters out genes that has total counts below 10
d0 <- d0[keep.exprs,, keep.lib.sizes=FALSE]


dim(d0)


#Continue analysis after Filtering
d0 <- calcNormFactors(d0, method = "TMM")
d0$samples$norm.factors


#Visualize using MDS plot
#From the sample report, we saw that there are some factors that may affect the data


replicates <- as.factor(c("Rep_1", "Rep_2", "Rep_3","Rep_4", "Rep_5", "Rep_6", "Rep_1", "Rep_2","Rep_3"))
d0$samples$replicates <- replicates


d0$samples$State <- as.factor(info$State)
State <- as.factor(c("Treated", "Treated", "Treated", "Treated", "Treated", "Treated", "Untreated", "Untreated","Untreated"))


#Visualize using MDS plot

lcpm <- cpm(d0, log=TRUE)
cpm <- cpm(d0)
s
pdf("MDS PLOT Scramble vs gapmer_2_5.pdf", height = 6, width = 7) #creating a pdf file

plotMDS(lcpm, labels=d0$samples$replicates, cex=1.5, dim = c(1,2), col = as.numeric(d0$samples$State)) #Please check with different labels and col. You should also look for other dimensions(i,e dim=c(2,3))
dev.off()

design <- model.matrix(~0+group)
colnames(design) <- gsub("group","", colnames(design))
design
contr.matrix <- makeContrasts(TreatedvsUntreated = Treated-Untreated, levels = colnames(design)) #
contr.matrix


v <- voom(d0, design, plot=TRUE) #DATA DISTRIBUTION VISUALIZATION
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
summary(decideTests(efit))

save(d0, vfit,v, file = "limma df with batch.RDATA")
#diff expr matrix table
# efit' is the result obtained from the limma pipeline
# 'TreatedvsUntreated' is the coefficient of interest


top.table <- topTable(efit, n = Inf, coef = "TreatedvsUntreated")
colnames( top.table)

# Customize column names
new_column_names <- c( "Identifier","Scramble_vs_Gapmer_2_5_LogFC", "AveExpr", "t", "P.Value", "Scramble_vs_Gapmer_2_5_adj_P_Val", "B")
colnames(top.table) <- new_column_names

# Print the modified differential expression matrix table with new column names
print(top.table)

write.csv2(top.table,"Scramble VS Gapmer 2 and 5 RP5 KD.csv")






