library(Seurat)
library(SCMarker)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
data.dir = args[2]
n.marker = args[3]

print("The work directory is:")
print(work.dir)

#print("The original 10x count data is at:")
#print(data.dir)

print("The number of marker chosen for each group:")
n.marker <- as.integer(n.marker)
print(n.marker)
data.dir <- file.path(work.dir, "/data/10x")
# this code will filter gene and cell by itself. Should use the 10x count directly.
pbmc.data <- Read10X(data.dir = data.dir)

# get highest expressive gene for each cluster
print("Create a Seurat object...")
cluster.dir <- file.path(work.dir, "data/cluster_labels.csv")
cluster <- read.table(
  file = cluster.dir,
    sep="\t",
  as.is = TRUE
)

marker.dir <- file.path(work.dir, "marker")
if (!file.exists(marker.dir)) {
dir.create(marker.dir)
}

pbmc <- CreateSeuratObject(counts = pbmc.data)
Idents(pbmc) <- cluster[2]
group_name <- unique(cluster[,2])

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

avg_express <- AverageExpression(object = pbmc, slot = 'counts')
df <- data.frame(avg_express)
print("The average expression of selected marker genes")
#print(df)

print("For each column (cluster), Order the genes by decreasing the expression")
res_exp <- lapply(1:dim(df)[2], function(col.number) df[order(df[, col.number], decreasing=TRUE)[1:n.marker], col.number, drop = FALSE])

print("The ordered genes for each group")
#print(res_exp)

all.marker <- c()
df.per.group <- data.frame(Cluster=character(0), Marker=character(0))
for (group in 1:dim(df)[2]) {
    res_g <- data.frame(res_exp[group])

    top_i <- rownames(res_g)
    print(group_name[group])
    print(top_i)
    df.per.group[nrow(df.per.group)+1, ] <- c(group_name[group], paste(top_i, collapse = ", "))
    all.marker <- c(all.marker, top_i)
}

all.marker <- unique(all.marker)
# this algorithm does not have a ranking for the top ones.
# this algorithm also does not have a parameter to decide
# smaller number of k gives less markers
write.table(all.marker, file.path(marker.dir, "marker_genes.txt"), sep="\t", col.names=FALSE)
write.csv(res_exp, file.path(marker.dir, "top_ones_per_group.cvs"))
write.csv(df.per.group, file.path(marker.dir, "marker_gene_per_group.csv"),row.names = FALSE)
