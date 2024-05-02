library(Seurat)
library(SCMarker)

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
n.marker = args[2]
kval = args[3]
nval = args[4]

print("For nearest neighbour algorithm, the used k & n values are:")
print(paste(kval, nval))

print("The work directory is:")
print(work.dir)

#print("The original 10x count data is at:")
#print(data.dir)

print("The number of marker chosen for each group:")
n.marker <- as.integer(n.marker)
print(n.marker)
######################## using 10x data
#data.dir <- file.path(work.dir, "/data/10x")
# this code will filter gene and cell by itself. Should use the 10x count directly.
#pbmc.data <- Read10X(data.dir = data.dir)
#pbmc.data
#data <- as.matrix(pbmc.data)
#head(data)
#data = log(data+1)

################# using scaled data of seurat
matrix.dir <- file.path(work.dir, "data/scale_data_high_variable.csv")
cluster.dir <- file.path(work.dir, "data/cluster_labels.csv")

print("Loading data")
pbmc.data <- read.table(
                       file = matrix.dir,
                       as.is = TRUE
)
data <- as.matrix(pbmc.data)

####################

print("the optimized bw is:")
thumb_one_bw = bw.nrd0(data)
print(thumb_one_bw)
#bw = 0.05/thumb_one_bw
bw = 1
print(bw)
start_time <- Sys.time()

# reference:  https://doi.org/10.1371/journal.pcbi.1007445
# The first step is to fit the data into a multi-modal, which is find out by the code itself, expressed as the form of Guassian Kernal function.
# Formular (1)
# width is a smoothing parameter called bandwidth: h
# firstly, the genes and cells are filtered
# for each gene, at least in cellK number of cell this gene expresses
# for each cell, at least geneK genes express
# only genens with multinormial distribution are considered as markers.eg. more than 1 peakg
# larger width is going to filter more gene and left with fewer genes
print("Filter the genes that are not multinormial...")
#head(data)
res=ModalFilter(data=data,geneK=10,cellK=10,width=bw)
# output is a object with atribute of geneSum, which contains the selected genes and number of peaks
str(res)
# The second step is just to filter some genes has expression higher than max-expression criteria.
print("Filter the genes that have high expression")
res=GeneFilter(obj=res)
#str(res)
#print(res$geneSumm)
#print(nrow(res$geneSumm))
# the third step will plot KNCEN or the KNMEN graph
# based on the 2 graphs, select the ones that are strongly correlated or mutually exclusive.
# k is the k number of nearest neighbour, smaller k give more marker genes
# then calculate Sij, which means for each pair of gene_i and gene_j, how many cells they co-occur.
# n is the criteria of S. S>=n in formular (3), which means the gene co-occurs in at least n cells
# among all genes occured (connected nodes), select the genes with largest Sij value for each gene gi).
# n is a special parameter, when n is large enough, smaller n value will provide stronger limits, and result
# in smaller pool. However, it seems the number of selected genes does not change too much with n.
# but if n is small enough, it seems the limits become stronger, and less genes are selected. And smaller
# value of n does not change the pool any more
print("Run KNN to decide genes...")
res=getMarker(obj=res,k=kval,n=nval)

print("Finish searching!!")
end_time <- Sys.time()
sprintf("Excuting time:%s second", end_time - start_time)
str(res)

# select the result of high ranking genes
# save markers
marker.dir <- file.path(work.dir, "marker")
if (!file.exists(marker.dir)) {
dir.create(marker.dir)
}

#print(res)
#head(res$marker)
#print(res$marker)
write.table(res$marker, file.path(marker.dir, "all_genes.txt"), sep="\t", col.names="markers")

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

pbmc <- ScaleData(pbmc, features = res$marker)

avg_express <- AverageExpression(object = pbmc, slot = 'counts', features=res$marker)
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
