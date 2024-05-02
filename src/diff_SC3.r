# SC3 is installed in yingluR4.1
library(Seurat)
library(SC3)

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

print("Loading data...")
pbmc_raw <- read.table(
  file = matrix.dir,
  as.is = TRUE
)

cluster <- read.table(
  file = cluster.dir,
    sep="\t",
  as.is = TRUE
)

# selecting marker genes by SC3

print("Selecting markers...")

start_time <- Sys.time()
group_name <- unique(cluster[,2])
mapped_values <- as.numeric(factor(cluster[,2], levels = group_name))

d <- get_marker_genes(as.matrix(pbmc_raw),  mapped_values)
head(d)
rownames(d) <- rownames(pbmc_raw)
########################################################################################
# this is a function of selecting markers from the original code os SC
# to select the genes based on cretiria of p_val and auroc
# and to rank the genes based on auroc
#' @param object an object of \code{SingleCellExperiment} class
#' @param k number of cluster
#' @param p_val p-value threshold
#' @param auroc area under the ROC curve threshold
#'
organise_marker_genes <- function(dat, p_val=0.05, auroc=0.4) {
   # dat <- rowData(object)[, c(paste0("sc3_", k, "_markers_clusts"), paste0("sc3_", k,
   #     "_markers_auroc"), paste0("sc3_", k, "_markers_padj"), "feature_symbol")]
    dat <- dat[dat[, 'pvalue'] < p_val & !is.na(dat[, 'pvalue']), ]
    dat <- dat[dat[, "auroc"] > auroc, ]

    d <- NULL

    for (i in sort(unique(dat[, "clusts"]))) {
        tmp <- dat[dat[, "clusts"] == i, ]
        tmp <- tmp[order(tmp[, "auroc"], decreasing = TRUE), ]
        d <- rbind(d, tmp)
    }

    if(nrow(dat) > 0) {
        return(d)
    } else {
        return(NULL)
    }
}

###############################################################################################
# get marker genes
print("Select marker genes...")
marker_sc3 <- organise_marker_genes(d)
head(marker_sc3)

end_time <- Sys.time()
sprintf("Excuting time:%s seconds", end_time - start_time)

# get top markers

marker.dir <- file.path(work.dir, "marker")
if (!file.exists(marker.dir)) {
dir.create(marker.dir)
}
top_list<-c()
df.per.group <- data.frame(Cluster=character(0), Marker=character(0))
for (group in unique(marker_sc3$clusts)){
    print(group)
    tmp <- marker_sc3[marker_sc3[, "clusts"] == group, ]
    tmp <- rownames(tmp)

    top_i<-tmp[1:n.marker]
    print(group_name[group])
    print(top_i)
    df.per.group[nrow(df.per.group)+1, ] <- c(group_name[group], paste(top_i, collapse = ", "))

    top_list<-c(top_list,top_i)
}
print("The top marker genes are:")
top_list <- top_list[!is.na(top_list)]
markers <- unique(top_list)
print(markers)

# save markers
write.csv(marker_sc3, file.path(marker.dir, "top_ones_per_group.csv"))
write.table(markers, file.path(marker.dir, "marker_genes.txt"), sep="\t", col.names=FALSE)
write.csv(df.per.group, file.path(marker.dir, "marker_gene_per_group.csv"),row.names = FALSE)
