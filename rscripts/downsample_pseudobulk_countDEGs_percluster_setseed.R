#-----------------------------------------------------------------
# Script performs DGEA on DOWNSAMPLED sample x cluster pseudobulk data, per cluster
# No interaction term
# This version runs anaysis on all clusters with > 1000 cells and removes low count features first
#
# Usage: supply seed as input variable when running script/submitting job to cluster
# use different seeds to get different downsampling runs
#
#
# Note: Returns a tibble with clusters as rows and reps as columns, each cell containing the number of genes differentially expressed between *ecotypes* (G genes)
#-----------------------------------------------------------------

#We're using a bash variable (that was given to the sbatch call as an argument. complex stuff.)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

print(args)
seed <- as.double(args[1])
print(seed)


library(Seurat)
library(Signac)
library(tidyverse)
library(DESeq2)

SeuratObjectPath <- "/path/to/data.rds"

#-----------------------------------------------------------------
#### Define functions ####

# function returns number of DEGs for ecotype contrast for a specified cluster 
# Note: no interaction term 
# Note: biological replicate is included as covariate
# Note: sex/sex ratio is *not* included

count_DGEs <- function(cluster=NA, pseudobulk,sample_info,filter=T){
    if(is.na(cluster)){message("No cluster specified")}
    message("Counting DGEs for cluster ", cluster)

    col_names <- colnames(pseudobulk)

    #### Subset pseudobulk count table
    columns_select <- col_names[str_detect(col_names, paste0("_",cluster,"$"))]
    pseudobulk_sub <- as.data.frame(pseudobulk[, columns_select])
    sample_info_sub <- sample_info[which(sample_info$cluster == cluster),]

    pseudobulk_sub <- pseudobulk_sub[which(rowSums(pseudobulk_sub) > 15),]

    ### Create dds object
    dds <- DESeqDataSetFromMatrix(countData = pseudobulk_sub, colData = sample_info_sub, design = ~ ecotype + treatment + bio_rep)
    mod_mat <- model.matrix(design(dds), colData(dds))

    dds$ecotype <- relevel(dds$ecotype, ref="marine")
    dds$treatment <- relevel(dds$treatment, ref="HS")

    #### Run DESeq
    dds <- DESeq(dds, test="Wald")
    res_ecotype <- lfcShrink(dds, contrast=c("ecotype","freshwater","marine"), type="ashr")
    summary(res_ecotype, alpha=0.05) 
    
    if(filter==T){
	nz <- pseudobulk_sub > 0
        nonzero_per_group <- sapply(
              split(colnames(pseudobulk_sub), sample_info_sub[["group"]]),
              function(cols) rowSums(nz[, cols, drop = FALSE])
            ) %>% as.data.frame()

    	filter_ecotype <- case_when(nonzero_per_group$marine_HS >= 3 | nonzero_per_group$marine_LS >= 3 | nonzero_per_group$freshwater_HS >=3 | nonzero_per_group$freshwater_LS >=3 ~ 1, TRUE ~ NA)
    	message(sum(is.na(filter_ecotype)), " genes filtered due to too few non-zero samples in ecotype contrast")
    	res_ecotype$padj <- res_ecotype$padj * filter_ecotype
    }

    x <- sum(res_ecotype$padj < 0.05, na.rm=T)
    return(x)
}

# function downsamples Seurat Object, aggregates data and runs count_DGEs() across all specified clusters
# returns a vector
run_percluster <- function(seed=1, clusters=NA, data, sample_info){ # clusters needs to a vector
    # subset to clusters of interest (e.g. those with > 1000 cells)
    data <- subset(data, idents=clusters)
    data_downsampled <- subset(data, downsample=1000, seed=seed)
    # aggregate data into pseudobulk by sample and cluster
    pseudobulk <- AggregateExpression(data_downsampled, group.by=c("orig.ident","seurat_clusters"), assays=c("RNA"))
    sample_info_percluster <- data.frame(sample_cluster=colnames(pseudobulk$RNA))
    sample_info_percluster[c("sample","cluster")] <- str_split_fixed(sample_info_percluster$sample_cluster, '_', 2)
    sample_info_percluster <- left_join(sample_info_percluster, sample_info, by="sample")

    result <- purrr::map_dbl(clusters, count_DGEs, pseudobulk=pseudobulk$RNA, sample_info=sample_info_percluster)
    names(result) <- clusters
    return(result)
}

#-----------------------------------------------------------------
#### Load and process data ####

# load integrated Seurat object
data.integrated <- readRDS(SeuratObjectPath) 
# load sample info
sample_info <- read.table(file="sample_info.tsv", header=T)
sample_info$group <- paste0(sample_info$ecotype,"_",sample_info$treatment)
sample_info$sample <- as.factor(sample_info$sample)

# set which clusters to include
table(data.integrated$seurat_clusters) # note that even if > 1000 cells, the lower the number the more often the same cells are resampled!
clusters <- 2:37

df <- run_percluster(seed=seed, clusters=clusters,data=data.integrated, sample_info=sample_info) 

save(df, file=paste0("tmp",seed,".RObj"))

