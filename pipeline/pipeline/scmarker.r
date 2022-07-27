library(Seurat)
library(SCMarker)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
#data.dir = args[2]
kval = args[2]
nval = args[3]
n.marker = args[4]

print("The work directory is:")
print(work.dir)

#print("The original 10x count data is at:")
#print(data.dir)

print("The number of marker chosen for each group:")
n.marker <- as.integer(n.marker)
print(n.marker)
data.dir <- file.path(work.dir, "/data/10x")
print("For nearest neighbour algorithm, the used k & n values are:")
print(paste(kval, nval))
# this code will filter gene and cell by itself. Should use the 10x count directly.
pbmc.data <- Read10X(data.dir = data.dir)
#pbmc.data
data <- as.matrix(pbmc.data)
#head(data)
start_time <- Sys.time()

# reference:  https://doi.org/10.1371/journal.pcbi.1007445
# The first step is to fit the data into a multi-modal, which is find out by the code itself, expressed as the form of Guassian Kernal function.
# Formular (1)
# width is a smoothing parameter called bandwidth: h
# firstly, the genes and cells are filtered
# for each gene, at least in cellK number of cell this gene expresses
# for each cell, at least geneK genes express
# only genens with multinormial distribution are considered as markers.eg. more than 1 peakg
print("Filter the genes that are not multinormial...")
res=ModalFilter(data=data,geneK=10,cellK=10,width=2)
#str(res)
# The second step is just to filter some genes has expression higher than max-expression criteria.
print("Filter the genes that have high expression")
res=GeneFilter(obj=res)
#str(res)

# the third step will plot KNCEN or the KNMEN graph
# based on the 2 graphs, select the ones that are strongly correlated or mutually exclusive.
# k is the k number of nearest neighbour
# n is the criteria of S. S<=n in formular (3) 
print("Run KNN to decide genes...")
res=getMarker(obj=res,k=kval,n=nval)

print("Finish searching!!")
end_time <- Sys.time()
sprintf("Excuting time:%s second", end_time - start_time)
#str(res)
head(res$marker)

# select the result of high ranking genes
# save markers
marker.dir <- file.path(work.dir, "marker")
if (!file.exists(marker.dir)) {
dir.create(marker.dir)
}
# get highest expressive gene for each cluster
print("Create a Seurat object...")
cluster.dir <- file.path(work.dir, "data/cluster_labels.csv")
cluster <- read.table(
  file = cluster.dir,
    sep="\t",
  as.is = TRUE
)

pbmc <- CreateSeuratObject(counts = pbmc.data)
Idents(pbmc) <- cluster[2]
group_name <- unique(cluster[,2])

avg_express <- AverageExpression(object = pbmc, features=res$marker)
df <- data.frame(avg_express)
print(df)
res_exp <- lapply(1:dim(df)[2], function(col.number) df[order(df[, col.number], decreasing=TRUE)[1:n.marker], col.number, drop = FALSE])

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
