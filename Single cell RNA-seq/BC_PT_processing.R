library(Seurat)
library(tidyverse)
library(limma)
library(clustree)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(scCustomize)
library(SCpubr)
library(patchwork)
library(scCustomize)
library(reticulate)
library(harmony)
library(glmGamPoi)
#BiocManager::install('glmGamPoi')

ngc()
set.seed(123)

plot_dir <- "/home/erikk/data1/Taufeeq/Patient data/"
save_dir <- "/home/erikk/data1/Taufeeq/Patient data/14 BC/"

#Import QC filtered seurat object
seurat <- readRDS(paste0(save_dir, "Filtered_14-BC_PT.RDS"))


seurat <- SCTransform(seurat, vars.to.regress = c("percent_mito", "S.Score", "G2M.Score", "nCount_RNA"), verbose = F)
table(seurat@meta.data$subtype)
#checking for batch effects with seurat object
seurat <- RunPCA(seurat, verbose = F)

#Use harmony to remove batch effects.
#seurat <- RunHarmony(seurat, c("Sample"))

seurat <- RunUMAP(seurat, dims = 1:30, 
                  seed.use = 123,
                  #reduction = "harmony", 
                  verbose = F,
                  reduction = "pca")

#PCA plot
#open pdf file
pdf(paste0(plot_dir, "2nd14_BC_PCA.pdf"), height = 4, width = 6)

DimPlot_scCustom(seurat, group.by = "PatientNumber", reduction = "pca", label = F, figure_plot = T)
DimPlot_scCustom(seurat, group.by = "CellType", reduction = "pca",figure_plot = T)
DimPlot_scCustom(seurat, group.by = "subtype", reduction = "pca", figure_plot = T)
dev.off()
#UMAP 
#open pdf file
#pdf(paste0(plot_dir, "14_BC_Umap.pdf"), height = 6, width = 10)

patient_u_plot <- DimPlot_scCustom(seurat, group.by = "PatientNumber", reduction = "umap", figure_plot = T)
subtype_u_plot <- DimPlot_scCustom(seurat, group.by = "CellType", reduction = "umap", figure_plot = T)
celltype_u_plot <- DimPlot_scCustom(seurat, group.by = "subtype", reduction = "umap", figure_plot = T)

feature_plot <- FeaturePlot_scCustom(seurat,features = "RP5-1198O20.4", figure_plot = T)


pdf(paste0(plot_dir, "2nd14-bc_UMAP.pdf"), height = 4, width = 6)

#visualizing expression of gene of interest 
#open pdf file
pdf(paste0(plot_dir, "14_BC_RP5.pdf"), height = 4, width = 8)

x <- VlnPlot_scCustom(seurat, features = "RP5-1198O20.4", group.by = "CellType")
y <- VlnPlot_scCustom(seurat, features = "RP5-1198O20.4", group.by = "subtype")

w <- DotPlot_scCustom(seurat, features = "RP5-1198O20.4", group.by = "subtype")
v <- DotPlot_scCustom(seurat, features = "RP5-1198O20.4", group.by = "CellType")


pdf(paste0(plot_dir, "2nd14_BC_RP5.pdf"), height = 4, width = 6)
print(x)
print(y)
print(w )
print(v)
dev.off()


