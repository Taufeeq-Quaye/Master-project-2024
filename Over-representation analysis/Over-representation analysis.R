library(readxl)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(tidyverse)
library(GseaVis)
library(readxl)

#Importing data frame
gapmer5<- as.data.frame(read.csv2("Scramble VS Gapmer 5 RP5 KD.csv", as.is = T, check.names = F))
gapmer5 <- gapmer5[,-1]

gapmer2<- as.data.frame(read.csv2("Scramble VS Gapmer 2 RP5 KD.csv", as.is = T, check.names = F))
gapmer2 <- gapmer2[,-1]

#filtering out significant Genes

gapmer2 <- gapmer2 %>% filter(Scramble_vs_Gapmer_2_adj_P_Val < 0.05)
gapmer5 <- gapmer5 %>% filter(Scramble_vs_Gapmer_5_adj_P_Val < 0.05)

#stratifying data into positive and negatively regulated based on LogFC

gp2_positive <- gapmer2 %>% filter(Scramble_vs_Gapmer_2_LogFC > 0)
gp5_positive <- gapmer5 %>% filter(Scramble_vs_Gapmer_5_LogFC > 0)

gp2_negative <- gapmer2 %>% filter(Scramble_vs_Gapmer_2_LogFC < 0)
gp5_negative <- gapmer5 %>% filter(Scramble_vs_Gapmer_5_LogFC < 0)

#Merge positively regulated genes together in order to find the common genes.
gp2_gp5_positive <- merge.data.frame(gp2_positive, gp5_positive, by="Identifier")
write.csv2(gp2_gp5_positive, "GP1 AND 2 Positve lOgFC common genes ")

gp2_gp5_negative <- merge.data.frame(gp2_negative,gp5_negative, by="Identifier")
write.csv2(gp2_gp5_negative, "GP1 AND 2 Negative lOgFC common genes")

#perform the enrichment analyses using the enrich r pipeline.
 #import data for background genes. in this case RP_KD GENES

#FIRST LETS START WITH THE POSITVE LOGFC GENES
KD <- as.data.frame(read_xlsx("RP5 KO.xlsx"))
rownames(KD) <- KD$Identifier

gp2_gp5_positive<- distinct(gp2_gp5_positive,Identifier,.keep_all = T) #remove duplicates
rownames(gp2_gp5_positive) <- gp2_gp5_positive$Identifier


significant_genes <- rownames(gp2_gp5_positive)
background_genes <- rownames(KD)

m_t2g <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%   #For Hallmark
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name, ensembl_gene)

##We will export the dotplot in a pdf for better view 
pdf("EnrichR Positive GP1,GP2 LOGFC Genes.pdf",  width=10,height=10)

enricher_hallmark <- enricher(significant_genes, TERM2GENE=m_t2g, 
                              universe = background_genes,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# Generate the dot plot
dot_plot <- dotplot(enricher_hallmark)

# Add a title
dot_plot <- dot_plot + ggtitle("Enrichment Analysis for Gapmer 1 & 2 Positive LogFC ") + theme(plot.title = element_text(face = "bold", size = 16))

# Display the plot
print(dot_plot)

dev.off()

#SECOND LETS MAKE WITH THE NEGATIVE LOGFC GENES

gp2_gp5_negative<- distinct(gp2_gp5_negative,Identifier,.keep_all = T) #remove duplicates
rownames(gp2_gp5_negative) <- gp2_gp5_negative$Identifier


significant_genes <- rownames(gp2_gp5_negative)
background_genes <- rownames(KD)

m_t2g <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%   #For Hallmark
  mutate(gs_name=gsub("HALLMARK_","",gs_name)) %>%
  dplyr::select(gs_name, ensembl_gene)

##We will export the dotplot in a pdf for better view 
pdf("EnrichR Negative GP1,GP2 LOGFC Genes.pdf",  width=10,height=10)

enricher_hallmark <- enricher(significant_genes, TERM2GENE=m_t2g, 
                              universe = background_genes,pvalueCutoff = 0.05,qvalueCutoff = 0.05)

# Generate the dot plot
dot_plot <- dotplot(enricher_hallmark)

# Add a title
dot_plot <- dot_plot + ggtitle("Enrichment Analysis for Gapmer 1 & 2 Negative LogFC ") + theme(plot.title = element_text(face = "bold", size = 16))

# Display the plot
print(dot_plot)

dev.off()
