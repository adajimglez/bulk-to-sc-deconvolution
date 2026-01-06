#######################################################################
# Code adapted from:
# Marquez-Galera et al., 2022
# Title: A protocol to extract cell-type-specific signatures from differentially expressed genes in bulk-tissue RNA-seq
# DOI: https://doi.org/10.1016/j.xpro.2022.101121
#
# Modifications and updates have been implemented to:
# - Adapt the workflow to current Seurat versions
# - Generate correlation heatmaps and export gene counts per cell
# - Customize plotting for specific dataset
#######################################################################


library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(cowplot)
library(data.table)
library(ggplot2)
library(gplots)
library(patchwork)
library(loomR)

setwd("D:/Data/RNAseq/")


############### Import 10X data after conversion with python instead of converting with SeuratDisk #############

raw_data<- Read10X(data.dir = "matrix_files_5dpf")
metadata <- read.csv("metadata.csv")
head(metadata)

rownames(metadata)<-metadata$cell_id

sobj <- CreateSeuratObject(counts = raw_data, meta.data = metadata)

sobj[[]]


#read DEG list
DEG<-read.csv("Up_64uM_ZnCl2_fullset.csv")
View(DEG)

#rank based on logfc
DEG<- DEG[order(-DEG$log2FoldChange),]

#extract gene symbols
DEG<-DEG$SYMBOL

#consult the number of cells per category
table(sobj$zebrafish_anatomy_ontology_class)

Idents(sobj) <- "zebrafish_anatomy_ontology_class"
Idents(sobj)<-factor(Idents(sobj), levels = sort(levels(sobj)))

#normalize gene counts
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify genes that are outliers on a "mean variability plot" (look fot this)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)

#Scale and center variable genes
sobj <- ScaleData(sobj)

#perform linear dimensionaltiy reduction by PCA over the variable genes
sobj <- RunPCA(sobj)

#calculate the number of principal components that are informative using the elbow plot (here I concluded around 25)
ElbowPlot(sobj, ndims = 50)

#perform non linear dimensionality redution over the PCs determined in the previous step (25) and UMAP
sobj<- RunTSNE(sobj, dims = 1:25)

sobj<- RunUMAP(sobj, dims = 1:25)

# #at this point the scell dataset can be saved...
# saveRDS(sobj, "Single-cell_custom_subset.rds")
# 
# #... and reloaded when needed
# sobj <- readRDS("Single-cell_custom_subset.rds")

#then plots can be created to show the single cell data and the dif populations

#without labels
pca_plot<-DimPlot(sobj, reduction = "pca", pt.size = 0.1, label = F, cols = topo.colors(38))
ggsave(filename= "pca.png", plot = pca_plot, width = 15, height = 10, dpi = 300)

tsne_plot<-DimPlot(sobj, reduction = "tsne", pt.size = 0.1, label = F, cols = rainbow(38))
ggsave(filename= "tsne.png", plot = tsne_plot, width = 15, height = 10, dpi = 300)

umap_plot<-DimPlot(sobj, reduction = "umap", pt.size = 0.1, label = F, cols = rainbow(38))
ggsave(filename= "umap.png", plot = umap_plot, width = 15, height = 10, dpi = 300)


#with labels
pca_plot2<-DimPlot(sobj, reduction = "pca", pt.size = 0.1, label = TRUE, cols = topo.colors(38))
ggsave(filename= "pca2.png", plot = pca_plot2, width = 15, height = 10, dpi = 300)

tsne_plot2<-DimPlot(sobj, reduction = "tsne", pt.size = 0.1, label = TRUE, cols = rainbow(38))
ggsave(filename= "tsne2.png", plot = tsne_plot2, width = 15, height = 10, dpi = 300)

umap_plot2<-DimPlot(sobj, reduction = "umap", pt.size = 0.1, label = TRUE, cols = rainbow(38))
ggsave(filename= "umap2.png", plot = umap_plot2, width = 15, height = 10, dpi = 300)

#the variable "cols" can be defined beforehand assigning colours to cell types. I set it up to the predefined topo.colors or rainbow

#### extract gene names from single cell data
gene_names <- rownames(sobj)

#### intersect bulk DEGs with single cell data
DEG <- DEG[DEG %in% gene_names]
DEG <- head(DEG,153)
topDEGs_list<- DEG
sobj <- ScaleData(sobj, features = DEG)
sobj <- RunPCA(sobj, features = DEG, approx=FALSE)
topDEGs_list<- rownames(sobj@reductions[["pca"]]@feature.loadings)
topdeg <-sobj@reductions[["pca"]]@feature.loadings

pca_plot3<-DimPlot(sobj, reduction = "pca", pt.size = 0.1, label = F, cols = rainbow(38))
pca_plot4<-DimPlot(sobj, reduction = "pca", pt.size = 0.1, label = T, cols = rainbow(38))
legend2<- get_legend(pca_plot3)
pca_plot3<- pca_plot3 & NoLegend()
ggsave(filename = "pca_up_64uMZnCl2.png", plot = pca_plot4, width = 15, height = 10, dpi = 300)
ggsave(filename = "pca_up_64uMZnCl2_legend.png", plot = legend2, width = 2, height = 3.75, dpi = 300)

sobj<- RunTSNE(sobj, features = DEG, approx=FALSE, check_duplicates = FALSE)

sobj<- RunUMAP(sobj, features = DEG, approx=FALSE, check_duplicates = FALSE)

#TSNE after intersection with labels
tsne_plot3<-DimPlot(sobj, reduction = "tsne", pt.size = 0.1, label = TRUE, cols = rainbow(38))
ggsave(filename= "tsne_up_64uMZnCl2.png", plot = tsne_plot3, width = 20, height = 15, dpi = 300)

#TSNE after intersection without labels
tsne_plot4<-DimPlot(sobj, reduction = "tsne", pt.size = 0.1, label = F, cols = rainbow(38))
ggsave(filename= "tsne_up_64uMZnCl2_nolab.png", plot = tsne_plot4, width = 20, height = 15, dpi = 300)

#UMAP after intersection with labels
umap_plot3<-DimPlot(sobj, reduction = "umap", pt.size = 0.1, label = TRUE, cols = rainbow(38))
ggsave(filename= "umap_up_64uMZnCl2.png", plot = umap_plot3, width = 20, height = 15, dpi = 300)

#UMAP after intersection without labels
umap_plot4<-DimPlot(sobj, reduction = "umap", pt.size = 0.1, label = F, cols = rainbow(38))
ggsave(filename= "umap_up_64uMZnCl2_nolab.png", plot = umap_plot4, width = 20, height = 15, dpi = 300)


### Correlation plots ####

random.matrix <- matrix(runif(500, min = -1, max = 1), nrow = 50)
quantile.range <- quantile(random.matrix, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["35%"], quantile.range["83%"], 0.06)
color.palette <- colorRampPalette(c("#0571b0","#f7f7f7","#ca0020"))(length(palette.breaks)-1)
clustFunction<-function(x)
hclust(as.dist(1-cor(t(as.matrix(x)), method = "pearson")), method = "average")


heatmapPearson <- function(correlations) {
  heatmap.2(
    x = correlations,
    col = color.palette,
    breaks = palette.breaks,
    trace = "none",
    symm = TRUE,
    hclustfun = clustFunction
  )
}

expr <- GetAssayData(sobj, assay = "RNA", layer = "data")


correlations_DEGs_log <- cor(
  log2(t(as.matrix(expr[topDEGs_list, ])) + 1),
  method = "pearson"
)

pdf(file = "Corr_plot.pdf", width = 30, height = 30)
heatmapPearson(correlations_DEGs_log)
dev.off()


#Getting metadata associated to each cell identifier to replace it by tissue name
expr_DEGs <- GetAssayData(sobj, assay = "RNA", slot = "counts")[topDEGs_list, ]
expr_DEGs <- as.data.frame(expr_DEGs)


expr_DEGs <- cbind(tissue = sobj@meta.data[colnames(expr_DEGs), "zebrafish_anatomy_ontology_class"], t(expr_DEGs))

write.table(expr_DEGs, file = 'Gene_Count_per_Cell_with_Tissue.tsv',
            quote = FALSE, sep = '\t', col.names = NA)
