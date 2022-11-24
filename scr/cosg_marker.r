# need to load COSG first, then Seurat. Ortherwise it is unable to work. Weired though...
library(COSG)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
n.marker = args[2]

print("The work directory is:")
print(work.dir)
print("The number of marker that will be find out for each group:")
n.marker <- as.integer(n.marker)
print(n.marker)

matrix.dir <- file.path(work.dir, "data/scale_data_high_variable.csv")
cluster.dir <- file.path(work.dir, "data/cluster_labels.csv")

print("Loading data")
pbmc_raw <- read.table(
		       file = matrix.dir,
		       as.is = TRUE
)

cluster <- read.table(
  file = cluster.dir,
    sep="\t",
  as.is = TRUE
)

print(cluster)
print("data loaded")

# create a seurat object. COSG is depended on the data structure of Seurat object.
print("Create a Seurat object...")
pbmc <- CreateSeuratObject(counts = pbmc_raw)
Idents(pbmc) <- cluster[2]
str(pbmc)
group_info <- Seurat::Idents(object = pbmc)
print(group_info)

# selection of marker genes
print("Selecting marker genes...")
#n.marker <- 10
start_time <- Sys.time()
marker_cosg<-cosg(
  pbmc,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=n.marker)

end_time <- Sys.time()
sprintf("Excuting time:%s seconds", end_time - start_time)
# get the top genes
print("Saving the top genes...")
marker.dir <- file.path(work.dir, "marker")
print(marker.dir)
if (!file.exists(marker.dir)) {
dir.create(marker.dir)
}
top_list<-c()
df.per.group <- data.frame(Cluster=character(0), Marker=character(0))
for (group in colnames(marker_cosg$names)){
    top_i<-marker_cosg$names[group][1:n.marker,1]

    df.per.group[nrow(df.per.group)+1, ] <- c(group, paste(top_i, collapse = ", "))
    top_list<-c(top_list,top_i)
}
print("The top marker genes are:")
markers <- unique(top_list)
print(markers)
# save marker
write.csv(marker_cosg, file.path(marker.dir, "top_ones_per_group.csv"))
write.table(markers, file.path(marker.dir, "marker_genes.txt"), sep="\t", col.names=FALSE)
write.csv(df.per.group, file.path(marker.dir, "marker_gene_per_group.csv"),row.names = FALSE)
