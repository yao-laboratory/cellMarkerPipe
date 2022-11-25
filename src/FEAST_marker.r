library(Seurat)
library(FEAST)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
n.marker = args[2]

martix.dir <- file.path(work.dir, "/data/norm_data_high_vriable.csv")
cluster.dir <- file.path(work.dir, "cluster/tabcluster.txt")

pbmc_raw <- read.table(
  file = martix.dir,
  as.is = TRUE
)

cluster <- read.table(
  file = cluster.dir,
    sep="\t",
  as.is = TRUE
)

F_res = cal_F2(pbmc_raw, cluster)
ixs = order(F_res$F_scores, decreasing = T) # order the features
top_list = rownames(pbmc_raw)[ixs][1:n.marker]

# save markers
#write.csv(marker_cosg, file.path(marker.dir, "top_ones_per_group.cvs"))
write.table(top_list, file.path(marker.dir, "marker_genes.txt"), sep="\t", col.names=FALSE)

