library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DropletUtils)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
data.dir = args[2]
n.variable = args[3]
n.variable = as.integer(n.variable)
if.known.marker = args[4]
if.known.marker = as.logical(if.known.marker)
if.cluster = args[5]
if.cluster = as.logical(if.cluster)

print("The work directory is:")
print(work.dir)
print("The data should be found at:")
print(data.dir)
print("The number of high variable features to choose:")
print(n.variable)
print("Whether do cluster?")
print(if.cluster)
print("Whether test if the known high variable genes are all included into the high variable genes")
print(if.known.marker)

#print(args[5])
if(length(args) < 5){
	n.PCA = 10 }else{
	n.PCA = strtoi(args[5])
	print("If choose to cluster, the number of PCA is going to be used:")
	print(n.PCA)
}
#"/work/yaolab/yinglu/project/selection/seurat"

# here we use data treated with 10X, barcodes.tsv  genes.tsv  matrix.mtx
# here the pbmc.data is a sparse matrix of saving gene*cell
pbmc.data <- Read10X(data.dir = data.dir)
# Initialize the Seurat project object with the raw (non-normalized data).
# here the count matrix is save as pbmc@assays$RNA@counts
# the min # of cell and min # of features are also defined here
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

if(!if.cluster){
	cluster.dir <- file.path(data.dir, "/groups.csv")
	cluster <- read.csv(
  	file = cluster.dir,
  	as.is = TRUE
	)

	Idents(pbmc) <- cluster[2]	
}

pbmc
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

Wdata.dir <- file.path(work.dir, "data")
dir.create(Wdata.dir)
print("Clearing data...")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

plot_distribution <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(Wdata.dir, "distribution_of_features_counts.png"), plot_distribution, width = 15, height = 10)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


dim(Idents(pbmc))
write.table(Idents(pbmc), file = file.path(Wdata.dir, "cluster_labels.csv"), sep="\t", col.names=FALSE)

# A-> ln(A/total_counts*scale.factor+1)
# the normalized data is saved in slot data (pbmc@assays$RNA@data) as a spase matrix
# or Normalized values are stored in pbmc[["RNA"]]@data

mat <- GetAssayData(object = pbmc, slot = 'counts')
print("Normalize data...")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## save the normalized matrix and output to file
#saveRDS(pbmc[["RNA"]]@data, file = file.path(work.dir, "data/norm_data.rds"))
#mat <- GetAssayData(object = pbmc, slot = 'counts')
#Scale the data
print("Scaling the data...")
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


#mat <- pbmc[["RNA"]]@data
#write.csv(mat,file = file.path(work.dir, "data/norm_data.csv")) # keeps the rownames
mat_scale <- GetAssayData(object = pbmc, slot = 'scale.data')
#write.table(mat,file = file.path(Wdata.dir, "norm_data.csv"),sep="\t") # keeps the rownames

#high variabel gene from cell to cell
print("Find High Variable Features...")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = n.variable)
mat_keep_rows <- head(VariableFeatures(object = pbmc),n.variable)
print("The number of high variable genes used:")
print(length(mat_keep_rows))

# Here to test whether the known genes are included into the high_variable genes


if(if.known.marker){
        include <- TRUE
        known.marker.dir <- file.path(data.dir, "Known_marker.csv")
        known.marker <- read.csv(
        file = known.marker.dir,
        as.is = TRUE,
        header = FALSE
        )

	print(known.marker)
        for (i in 1:nrow(known.marker)){
            genes <-strsplit(known.marker[i,2], ", ")[[1]]
	    print(genes)
            all_test <- genes %in% mat_keep_rows
            print(all_test)
	    print(genes %in% all.genes)
            if (FALSE %in% all_test){
                include <- FALSE
                break()
            }
        }

        # If not, need to ajusted the number of high_variable genes
        if(!include){
                stop("The known markers are not totally included into the high_variable genes, need to adjust the number of high_variable genes")
        } else {
                print("The known markers are well included into the high_variable genes")
        }
    }


mat_subset <- mat[rownames(mat) %in% mat_keep_rows, ]
mat_scale_subset <- mat_scale[rownames(mat_scale) %in% mat_keep_rows, ]
print("The dimension of the selected matrix is:")
print(dim(mat_subset))
print("Save the normalized data of high variables")
start_time <- Sys.time()
write.table(mat_subset,file = file.path(Wdata.dir, "counts_high_variable.csv"),sep="\t") # keeps the rownames
end_time <- Sys.time()
print("time used:")
print(end_time - start_time)
print("Save the scaled data of high variables")
start_time <- Sys.time()
write.table(mat_scale_subset,file = file.path(Wdata.dir, "scale_data_high_variable.csv"),sep="\t") # keeps the rownames
end_time <- Sys.time()
print("time used:")
print(end_time - start_time)
print("Save the 10x format of normalized data of high variables")
start_time <- Sys.time()
write10xCounts(path = file.path(Wdata.dir, "10x"), mat_subset)
end_time <- Sys.time()
print("time used:")
print(end_time - start_time)
print("Obtained High Variable Features!")

# cluster with KNN
if(if.cluster){
# make a directory to save cluster files
	print("Find cluster....")
	cluster.dir <- file.path(work.dir, "cluster/")
	dir.create(cluster.dir)
# dimensional reduction like PCA
# run PCA
	pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Visualize top genes associated with PCAs
	plot_dim_reduction <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
	ggsave(file.path(cluster.dir, "top_genes_of_PCA.png"), plot_dim_reduction, width = 15, height = 10)
# Determine the dimension of PCA used for cluster 
	pbmc <- JackStraw(pbmc, num.replicate = 100)
	pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#a ranking of principle components based on the percentage of variance explained by each one
# choosing the features before elbows
	plot_elbow <- ElbowPlot(pbmc)
	ggsave(file.path(cluster.dir, "elbow_rank_PCA.png"), plot_elbow, width = 15, height = 10)

# number of PCA should be determined based on the elbow plots, defalting value is 10
# KNN graph

	pbmc <- FindNeighbors(pbmc, dims = 1:n.PCA)
	pbmc <- FindClusters(pbmc, resolution = 0.5)

	print("Finish Cluster!!!")
	print("Plotting Figures & save cluster result...")
# visualize cluster using UMAP
	pbmc <- RunUMAP(pbmc, dims = 1:n.PCA)
	plot_UMAP <- DimPlot(pbmc, reduction = "umap")
	ggsave(file.path(cluster.dir, "UMAP.png"), plot_UMAP, width = 15, height = 10)

# export cluster result and UMAP axis
# save cluster name for each cell
	#write.csv(pbmc@ meta.data$ seurat_clusters,file = file.path(cluster.dir, "pbmc3k_cluster_labels.csv"), sep="\t", col.names=FALSE) # keeps the rownames
	write.table(pbmc[['seurat_clusters']],file = file.path(cluster.dir, "tabcluster.csv"), sep="\t", col.names=FALSE) # keeps the rownames
#	write.csv(pbmc@ active.ident,file = file.path(cluster.dir, "pbmc3k_cluster_names.csv"), sep="\t", col.names=FALSE) # keeps the rownames
	write.table(pbmc@ active.ident,file = file.path(cluster.dir, "tabnames.csv"), sep="\t", col.names=FALSE) # keeps the rownames
	# save UMAP embedding for plotting
	vis_data <- Embeddings(pbmc, reduction = "umap")
	write.table(vis_data, file.path(cluster.dir, "tabvis.txt"), sep="\t", col.names=FALSE)
	# save Seurat project
#	saveRDS(pbmc, file =file.path(work.dir, "data/pbmc3k_cluster.rds"))
}
