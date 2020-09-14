

# ************************** scRNA-seq Analysis **********************
# ************************** 14 September 2020 ************************

# Tutorial: https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
# Seurat user manual: https://cran.r-project.org/web/packages/Seurat/Seurat.pdf
# 10x datasets: 

setwd("/home/sarwar/scRNASeq/") # set your own directory
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "hg19/") # set your own data directory
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size

dense.size/sparse.size
# ************************* Preprocessing ************
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# -------- Normalizing the data ----------
pbmc <- NormalizeData(pbmc)

# ------Identification of highly variable features (feature selection)---
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# ------- Scaling the data
 all.genes <- rownames(pbmc)
 pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- ScaleData(pbmc) # sacling for highly varibale features only
pbmc[["RNA"]]@scale.data[1:2,1:3]

# --- Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) ## npcs = 50
DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 12:17, cells = 500, balanced = TRUE)

#pbmc <- JackStraw(pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

# --------- Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# ----- Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap", label = T, label.size = 7) + NoLegend()

pbmc <- RunTSNE(object = pbmc,reduction = "pca", dims.use = 1:10, do.fast = TRUE)
#head(x = Idents(object = pbmc), 5)
DimPlot(pbmc, reduction = "tsne", label = T, label.size = 7, pt.size = 1)

# Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
nrow(cluster1.markers)
head(cluster1.markers, n = 5)
cluster1.markers.sort.fc<- cluster1.markers[order(-cluster1.markers$avg_logFC),]
head(cluster1.markers.sort.fc, n = 5)

VlnPlot(pbmc, features = c("CCR7", "LDHB"))

nrow(filter(cluster1.markers.sort, p_val_adj<=0.05))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(pbmc, features = top5$gene) + NoLegend()

VlnPlot(pbmc, features = c("MS4A1", "CD79A", "PPBP"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

VlnPlot(pbmc, features = c("CD14","CD8A","CD79A"))
FeaturePlot(pbmc, features = c("CD14","CD8A","CD79A"))

FeatureScatter(pbmc, feature1 = "CD14", feature2 = "CD8A")
