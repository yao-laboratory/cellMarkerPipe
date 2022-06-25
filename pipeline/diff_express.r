library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
n.marker = args[2]

print("The work directory is:")
print(work.dir)
print("The number of marker that will be find out for each group:")
n.marker <- as.integer(n.marker)
print(n.marker)
#pbmc <- readRDS(file.path(work.dir, "data/pbmc3k_cluster.rds"))

# reading 1G matrix data is still a time consuming problem
#matrix.dir <- file.path(work.dir, "/data/norm_data_high_variable.csv")
#pbmc.raw <- read.table(
#  file = martix.dir,
#  as.is = TRUE
#)

# read 10x data directly
matrix.dir <- file.path(work.dir, "/data/10x")
pbmc.raw <- Read10X(data.dir = matrix.dir)

# here need to prove, the seperation is specific for my example
# should be standarlized using sep as "\t"
cluster.dir <- file.path(work.dir, "/data/cluster_labels.csv")
cluster <-  read.table(
  file = cluster.dir,
  sep = '\t',
  as.is = TRUE
)

# using the matrix and cluster data to create a Seurat Object
pbmc <- CreateSeuratObject(counts = pbmc.raw)
Idents(pbmc) <- cluster[2]

print("Creating Seurat Object:")
pbmc

# scaling data; have to scale data to process FindMarkers

print("Scaling data...")
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# KNN graph

#pbmc <- FindNeighbors(pbmc, dims = 1:10)
#pbmc <- FindClusters(pbmc, resolution = 0.5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
print("Try to find markers...")
start_time <- Sys.time()
pbmc.markers <- FindAllMarkers(pbmc,logfc.threshold = 0.25)
end_time <- Sys.time()
sprintf("Excuting time:%s second", end_time - start_time)

#print(pbmc.markers)
# top 10 for each group
print("Finish searching!!!")
topones <- pbmc.markers %>% group_by(cluster) %>% top_n(n = n.marker, wt = avg_log2FC)

print("Save figures...")
plot_hotmap <- DoHeatmap(pbmc, features = topones$gene) + NoLegend()

marker.dir <- file.path(work.dir, "marker")
dir.create(marker.dir)
ggsave(file.path(marker.dir, "marker_hot_map.png"), plot_hotmap, width = 15, height = 10)

print("The founded markers are:")
markers <- unique(topones$gene)
print(markers)

# save markers
print("Save the markers...")
write.csv(topones, file.path(marker.dir, "top_ones_per_group.cvs"))
write.table(markers, file.path(marker.dir, "marker_genes.txt"), sep="\t", col.names=FALSE)
df.per.group <- data.frame(Cluster=character(0), Marker=character(0))
for (group in unique(Idents(pbmc))){
	print(group)
	marker_per_group <- topones[which(topones$cluster == group), ]$gene
	print(marker_per_group)
	df.per.group[nrow(df.per.group)+1, ] <- c(group, paste(marker_per_group, collapse = ", "))
}
write.csv(df.per.group, file.path(marker.dir, "marker_gene_per_group.csv"),row.names = FALSE)
