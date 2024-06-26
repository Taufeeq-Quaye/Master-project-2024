library(Seurat)
library(tidyverse)
library(limma)
library(clustree)
library(cowplot)
library(sctransform)
library(scCustomize)
library(SCpubr)
library(patchwork)

sessionInfo()
#install.packages("igraph")
set.seed(123)
#set directory to save files
plot_dir <- "/home/erikk/data1/Taufeeq/Patient data/"
save_dir <- "/home/erikk/data1/Taufeeq/Patient data/14 BC/"

dose.data <- Read10X("/home/erikk/data1/Taufeeq/Patient data/14 BC")
#load metadata and insert into seurat object
metadata <- read.csv2("metadata_E-MTAB-8107.csv", sep = ",", as.is = T, check.names = F, row.names = 1)
seurat <- CreateSeuratObject(dose.data, project = "14_BC_PT", assay = "RNA", min.cells = 5, min.features = 200, meta.data = metadata)

#ADD QC METRICEES
seurat <- Add_Cell_QC_Metrics(seurat, add_top_pct = F, species = "human")

VlnPlot(seurat, features = c("CDH1","RP5-1198O20.4"))
Idents(seurat) <- seurat
VlnPlot(seurat, features = "percent_mito")

#The purpose of this code snippet is to compute a threshold for mitochondrial gene percentage in scRNA-seq data. 
#Cells with mitochondrial gene percentages exceeding this threshold may be flagged as potentially low-quality or damaged cells and could be filtered out from downstream analysis. 
#Identifying and removing such cells helps improve the quality of the scRNA-seq dataset and enhances the reliability of downstream analyses, such as cell clustering, trajectory inference, and differential expression analysis.
median_percent_MT <- median(seurat$percent_mito)
mad_percent_MT <- mad(seurat$percent_mito)
high_threshold_MT <- median_percent_MT + 5*mad_percent_MT

Idents(seurat) <- seurat@meta.data$TumorType

# QC plots and justification for cut off point
qc_plot <- QC_Plots_Combined_Vln(seurat_object = seurat, mito_cutoffs = high_threshold_MT, pt.size = 0.1) 
qc_plot_2 <- QC_Plot_UMIvsGene(seurat_object = seurat, meta_gradient_name = "percent_mito", meta_gradient_low_cutoff = high_threshold_MT)

histo_qc <- QC_Histogram(seurat_object = seurat, features = "percent_mito", low_cutoff = high_threshold_MT)
  
  # Save QC plots
  pdf(paste0(plot_dir, "14-BC_QC.pdf"), height = 7, width = 13)
print(qc_plot)
print(qc_plot_2)
print(histo_qc)
dev.off()

plot(qc_plot)
plot(qc_plot_2)
plot(histo_qc)

#WE SKIPPED THIS STEP BECAUSE THE DATA HAS ALREADY BEEN FILTERED
# Subset and process the data
#seurat <- subset(seurat, subset = percent_mito <= high_threshold_MT)


#The scDblFinder package gathers various methods for the detection and handling of doublets/multiplets in single-cell sequencing data 
#(i.e. multiple cells captured within the same droplet or reaction volume). 
library(scDblFinder)

#creating a single cell experiment object
nsce <- as.SingleCellExperiment(DietSeurat(seurat))

#annotation of singlet and doublets  in our Seurat object
sce <- scDblFinder(nsce)
seurat$Multiplet <- sce$scDblFinder.class
#Removal of of doublets
seurat <- subset(seurat, subset = Multiplet == "singlet")

# Save the processed Seurat object
saveRDS(seurat, file = paste0(save_dir, "Filtered_14-BC_PT.RDS"))
rm(seurat)
gc()
