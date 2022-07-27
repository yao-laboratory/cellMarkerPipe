library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
#data.dir = args[2]

print(work.dir)
#print(data.dir)

if(length(args) < 3){
	n.PCA = 10} else{
	n.PCA = strtoi(args[3]) }
#"/work/yaolab/yinglu/project/selection/seurat"

# here we use data treated with 10X, barcodes.tsv  genes.tsv  matrix.mtx
# here the pbmc.data is a sparse matrix of saving gene*cell

matrix.dir <- file.path(work.dir, "/data/10x")
cluster.dir <- file.path(work.dir, "/data/cluster_labels.csv")

print("the 10x data is here:")
print(matrix.dir)
# read 10x data directly
pbmc.raw <- Read10X(data.dir = matrix.dir)
rownames(pbmc.raw) <- toupper(rownames(pbmc.raw))
pbmc <- CreateSeuratObject(counts = pbmc.raw)

df.markers <- read.table(file.path(work.dir, "marker/marker_genes.txt"), sep="\t", col.names='markers')
markers <- df.markers$markers
# because Comet automatically change all genes to upper cases, so I need change all genes in the original matrix into upper case.
markers <- toupper(markers)

print(markers)
# Initialize the Seurat project object with the raw (non-normalized data).
# here the count matrix is save as pbmc@assays$RNA@counts
# the min # of cell and min # of features are also defined here
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# A-> ln(A/total_counts*scale.factor+1)
# the normalized data is saved in slot data (pbmc@assays$RNA@data) as a spase matrix
# or Normalized values are stored in pbmc[["RNA"]]@data
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# cluster with KNN
# make a directory to save re-cluster files
	cluster.dir <- file.path(work.dir, "re-cluster/")
	dir.create(cluster.dir)
#high cell-to-cell variation
#	pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# dimensional reduction like PCA
# if this step takes too long,use ScaleData(pbmc), only scale on high variable features
	pbmc <- ScaleData(pbmc, features = markers)
# PCA
	pbmc <- RunPCA(pbmc, features = markers)
# Visualize top genes associated with PCAs
	plot_dim_reduction <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
	ggsave(file.path(cluster.dir, "top_genes_of_PCA.png"), plot_dim_reduction, width = 15, height = 10)
# Determine the dimension of PCA used for cluster 
	#pbmc <- JackStraw(pbmc, num.replicate = 100)
	#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
	#print("still go")
#a ranking of principle components based on the percentage of variance explained by each one
# choosing the features before elbows
	#plot_elbow <- ElbowPlot(pbmc)
	#ggsave(file.path(cluster.dir, "elbow_rank_PCA.png"), plot_elbow, width = 15, height = 10)
# number of PCA should be determined based on the elbow plots, defalting value is 10
# KNN graph

	pbmc <- FindNeighbors(pbmc, dims = 1:n.PCA)
	pbmc <- FindClusters(pbmc, resolution = 0.5)
	print("still go")
# visualize cluster using UMAP
	pbmc <- RunUMAP(pbmc, dims = 1:n.PCA)
	print("I am alive")
	plot_UMAP <- DimPlot(pbmc, reduction = "umap")
	ggsave(file.path(cluster.dir, "UMAP.png"), plot_UMAP, width = 15, height = 10)

# export cluster result and UMAP axis
# save cluster name for each cell
	write.csv(pbmc@ meta.data$ seurat_clusters,file = file.path(cluster.dir, "re_cluster_labels.csv")) # keeps the rownames
	write.csv(pbmc@ active.ident,file = file.path(cluster.dir, "re_cluster_names.csv")) # keeps the rownames
	# save UMAP embedding for plotting
	vis_data <- Embeddings(pbmc, reduction = "umap")
	write.table(vis_data, file.path(cluster.dir, "vis.txt"), sep="\t", col.names=FALSE)
