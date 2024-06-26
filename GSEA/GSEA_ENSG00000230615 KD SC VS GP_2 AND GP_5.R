library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(dplyr)
library(stats)
library(readxl)
library(ggplot2)
library(viridis)
library(GseaVis)
library(enrichplot)
library(stringr)

res <- read.csv2("Scramble VS Gapmer 2 and 5 RP5 KD BATCH.csv",as.is = T, check.names = F, row.names = 1)
rownames(res) <- res$Identifier

significant_genes<- res %>%
  as.data.frame() %>%
  filter(Scramble_vs_Gapmer_2_5_adj_P_Val<= 0.05) %>% 
  rownames()


significant_genes_map<- clusterProfiler::bitr(geneID = significant_genes,
                                              fromType="ENSEMBL", toType="ENTREZID",
                                              OrgDb="org.Hs.eg.db")


## background genes are genes that are detected in the RNAseq experiment 
KD <- as.data.frame(read_xlsx("RP5 KO.xlsx"))
rownames(KD) <- KD$Identifier

background_genes<- KD %>% 
  as.data.frame() %>% 
  pull(Identifier)


background_genes_map<- bitr(geneID = background_genes, 
                            fromType="ENSEMBL", 
                            toType="ENTREZID",
                            OrgDb="org.Hs.eg.db")


### Gene set enrichment analysis#########################################################################################################
### HALLMARK DATABASE###################################################################################################################

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name,ensembl_gene)


res<- res %>% 
  mutate(signed_rank_stats = sign(Scramble_vs_Gapmer_2_5_LogFC) * -log10(Scramble_vs_Gapmer_2_5_adj_P_Val)) %>%
  arrange(desc(signed_rank_stats))

gene_list<- res$signed_rank_stats
names(gene_list)<- res$Identifier
set.seed(123)
em2 <- GSEA(gene_list, TERM2GENE=m_t2g, seed = 123)

dotplotGsea(em2)

nes <- em2[,c(1,5,7)]
nes <- nes %>% mutate(ID = gsub("HALLMARK_", "", ID),      # Remove "HALLMARK_" prefix
                      ID = gsub("_", " ", ID))
plot <- ggplot(nes, aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill= `p.adjust`)) +           # Fill by actual p.adjust values
  scale_fill_gradient(low = "#347eba", high="#df6664",
                      guide = "colorbar") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="ENSG00000230615 KD Scrambe VS Gapmer 2 and 5", fill="Adjusted p-value") + 
  theme_minimal()

plot

pdf("GSEA Hallmark ENSG00000230615 KD Scrambe VS Gapmer2 and 5.pdf",  width=8,height=5)
print(plot)
dev.off()

##REACTOME DATABASE.##########################################################################################################
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name,ensembl_gene)

res<- res %>% 
  mutate(signed_rank_stats = sign(Scramble_vs_Gapmer_2_5_LogFC) * -log10(Scramble_vs_Gapmer_2_5_adj_P_Val)) %>%
  arrange(desc(signed_rank_stats))

gene_list<- res$signed_rank_stats
names(gene_list)<- res$Identifier
set.seed(123)
em2 <- GSEA(gene_list, TERM2GENE=m_t2g, seed = 123)

dotplotGsea(em2)

nes <- em2[,c(1,5,7)]
nes <- nes %>% mutate(ID = gsub("REACTOME_", "", ID),      # Remove "REACTOME_" prefix
                      ID = gsub("_", " ", ID))

top_bottom_rows <- nes %>%  arrange(desc(NES)) #%>%  bind_rows(slice_head(n = 10), slice_tail(n = 10))
head <- slice_head(top_bottom_rows, n = 10)
tail <- slice_tail(top_bottom_rows, n= 10)
nes1 <- rbind(head,tail)

# Assuming nes1$ID contains the pathway names
nes1$Wrapped_ID <- str_wrap(nes1$ID, width = 40)  # Adjust width as needed

# Create the plot with wrapped pathway names
plot <- ggplot(nes1, aes(reorder(Wrapped_ID, NES), NES)) +
  geom_col(aes(fill = `p.adjust`)) +
  scale_fill_gradient(low = "#347eba", high = "#df6664", guide = "colorbar") +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "REACTOME DB ENSG00000230615 KD Scrambe VS Gapmer 2 and 5",
       fill = "Adjusted p-value") +
  theme_minimal()
#+
  #theme(axis.text.y = element_text(size = 12))

plot



pdf("GSEA REACTOME ENSG00000230615 KD Scrambe VS Gapmer 2 and 5.pdf",  width=10,height=7)
print(plot)
dev.off()

##BIOCARTA DATABASE ##########################################################################################################
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") %>% 
  dplyr::select(gs_name,ensembl_gene)


em2 <- GSEA(gene_list, TERM2GENE=m_t2g, seed = 123)

dotplotGsea(em2)

nes <- em2[,c(1,5,7)]
nes <- nes %>% mutate(ID = gsub("BIOCARTA_", "", ID),      # Remove "REACTOME_" prefix
                      ID = gsub("_", " ", ID))


plot <- ggplot(nes, aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill= `p.adjust`)) +           # Fill by actual p.adjust values
  scale_fill_gradient(low = "#347eba", high="#df6664",
                      guide = "colorbar") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="BIOCARTA DB ENSG00000230615 KD Scrambe VS Gapmer 2 and 5", fill="Adjusted p-value") + 
  theme_minimal()

plot

pdf("GSEA BIOCARTA ENSG00000230615 KD Scrambe VS Gapmer 2 and 5.pdf",  width=8,height=5)
print(plot)
dev.off()


##KEGG DATABASE#####################################################################################################################
#Perform GSEA using KEGG database
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name,ensembl_gene)


em2 <- GSEA(gene_list, TERM2GENE=m_t2g, seed = 123)

dotplotGsea(em2)

nes <- em2[,c(1,5,7)]
nes <- nes %>% mutate(ID = gsub("KEGG_", "", ID),      # Remove "REACTOME_" prefix
                      ID = gsub("_", " ", ID))

top_bottom_rows <- nes %>%  arrange(desc(NES)) #%>%  bind_rows(slice_head(n = 10), slice_tail(n = 10))
head <- slice_head(top_bottom_rows, n = 10)
tail <- slice_tail(top_bottom_rows, n= 10)
nes1 <- rbind(head,tail)

plot <- ggplot(nes1, aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill= `p.adjust`)) +           # Fill by actual p.adjust values
  scale_fill_gradient(low = "#347eba", high="#df6664",
                      guide = "colorbar") +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG DB ENSG00000230615 KD Scrambe VS Gapmer 2 and 5", fill="Adjusted p-value") + 
  theme_minimal()

plot

pdf("GSEA KEGG ENSG00000230615 KD Scrambe VS Gapmer 2 and 5.pdf",  width=10,height=7)
print(plot)
dev.off()

######################################################################################################################################################




em <- enricher(significant_genes, TERM2GENE=m_t2g,
               universe = background_genes)
dotplot(em)

#Visualization using cnet
ek <- enrichGO(gene = significant_genes,
               universe = background_genes,
               keyType = "ENSEMBL",
               ont = "BP",
               OrgDb = "org.Hs.eg.db")
dotplot(ek)
dot_plot <- dotplot(ek, title="PLD2 OE",showCategory=20, font.size=10) + theme(legend.position = "right",
                                                                               axis.text.x = element_text(size = 10),
                                                                               axis.text.y = element_text(size = 10),
                                                                               axis.title.x = element_text(size = 11),
                                                                               axis.title.y = element_text(size = 11),
                                                                               plot.title = element_text(hjust = 0.3, size = 12),
                                                                               panel.border = element_rect(colour = "black", fill = NA, 
                                                                                                           linewidth = 1))
dot_plot
pdf("PLD2 OE EnrichGO.pdf",  width=8,height=9)
print(dot_plot)
dev.off()
ggplot(nes %>% mutate(ID = gsub("HALLMARK_", "", ID)), aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill= `p.adjust` < 0.05)) +
  coord_flip() +
  scale_fill_manual(values= c("#347eba", "#df6664")) + 
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="D492 PLD2 OE vs Control") + 
  theme_minimal()
