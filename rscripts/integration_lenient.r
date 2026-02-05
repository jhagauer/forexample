#------------------------------------
# Updated integration script
# Takes in more leniently filtered data.list but with doublets removed
# Binarises ATAC-seq data 
#------------------------------------
# Load packages
library(Seurat)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
options(future.globals.maxSize = 100 * 1024^3)

# set path
projectpath="/path/to/project/"
# load functions
source(paste0(projectpath,"scripts/fullintegration_function.R"))
source(paste0(projectpath,"scripts/addlabels_function.R"))
# load filtered data list
data.list <- readRDS(paste0(projectpath,"output/data_filtered_doubletsremoved_lenientmt.rds"))

# merge technical replicates
LSup_a <- merge(data.list[["sample02"]],data.list[["sample08"]])
LSup_a <- JoinLayers(LSup_a)
LSdown_a <- merge(data.list[["sample04"]],data.list[["sample07"]])
LSdown_a <- JoinLayers(LSdown_a)
HSup_a <- merge(data.list[["sample03"]],data.list[["sample05"]])
HSup_a <-JoinLayers(HSup_a)
HSdown_a <- merge(data.list[["sample01"]],data.list[["sample06"]])
HSdown_a <- JoinLayers(HSdown_a)

LSup_b <- merge(data.list[["sample12"]],data.list[["sample15"]])
LSup_b <- JoinLayers(LSup_b)
LSdown_b <- merge(data.list[["sample10"]],data.list[["sample14"]])
LSdown_b <- JoinLayers(LSdown_b)
HSup_b <- merge(data.list[["sample09"]],data.list[["sample11"]])
HSup_b <- JoinLayers(HSup_b)
HSdown_b <- merge(data.list[["sample13"]],data.list[["sample16"]])
HSdown_b <- JoinLayers(HSdown_b)

LSup_c <- merge(data.list[["sample19"]],data.list[["sample22"]])
LSup_c <- JoinLayers(LSup_c)
LSdown_c <- merge(data.list[["sample17"]],data.list[["sample24"]])
LSdown_c <- JoinLayers(LSdown_c)
HSup_c <- merge(data.list[["sample21"]],data.list[["sample23"]])
HSup_c <- JoinLayers(HSup_c)
HSdown_c <- merge(data.list[["sample18"]],data.list[["sample20"]])
HSdown_c <- JoinLayers(HSdown_c)

# make new data list with the merged objects
data.list <- list("LSup_a" = LSup_a, "LSdown_a" = LSdown_a, "HSup_a" = HSup_a, "HSdown_a"=HSdown_a, "LSup_b" = LSup_b, "LSdown_b"=LSdown_b, "HSup_b"=HSup_b, "HSdown_b"= HSdown_b, "LSup_c"=LSup_c, "LSdown_c"=LSdown_c, "HSup_c"=HSup_c, "HSdown_c"=HSdown_c)

# define sample.tree that specifies the order in which samples will be integrated
sample.tree <- matrix(c(
    -1, -5,
    1, -9,
    -2,-6,
    3,-10,
    -3,-7,
    5,-11,
    -4,-8,
    7,-12,
    2,4,
    6,8,
    9,10
), nrow=2)
sample.tree <- t(sample.tree)

# run function
data.integrated <- RunFullIntegration(data.list=data.list, sample.tree=sample.tree, dims1=50, dims2=50, dims3=50, seed=404)
data.integrated <- add_labels(data.integrated) # adds metadata column for bio_rep, ecotype, treatment, etc.
# save as RDS object
saveRDS(data.integrated, paste0(projectpath,"output/data_integrated_updated_atacbinarized_lenientmt.rds"))