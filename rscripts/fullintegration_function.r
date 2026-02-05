# Function to run the entire integration
# takes merged.data.list and sample.tree as input

# dims1 = gex separate samples, dims2 = gex integrated, dims3 = atac lsi
RunFullIntegration <- function (samples = all, data.list, seed=404, sample.tree = NULL, dims1 = 50, dims2=50, dims3=50){
    # Part 1: Processing GEX data
    message("Start processing for GEX data")
    # Normalise data with sctransform, returning all features
    for(sample in names(data.list)){DefaultAssay(data.list[[sample]]) <- "RNA"}
    data.list <- lapply(data.list, SCTransform, variable.features.n = 5000, return.only.var.genes=F) #Seurat 5
    message("GEX data normalised.")

    # Run PCA
    for(sample in names(data.list)){DefaultAssay(data.list[[sample]]) <- "SCT"}
    data.list <- lapply(data.list, RunPCA, npcs=dims1, assay="SCT")
    message("GEX PCA finished.")

    # Integration across samples for GEX
    for(sample in names(data.list)){DefaultAssay(data.list[[sample]]) <- "SCT"}

    gex.features <- SelectIntegrationFeatures(object.list = data.list, verbose=T, nfeatures = 3000, method="sct")
    # the code below is technically only necessary if only variable genes were returned when running SCTransform, and therefore the sct slot may not contain the same genes in every sample
    shared.features <- Reduce(intersect, lapply(data.list, function(x) rownames(GetAssayData(x, slot = "scale.data", assay = "SCT"))))
    anchor.features <- intersect(gex.features, shared.features)

    data.list <- PrepSCTIntegration(
        object.list = data.list,
        anchor.features = anchor.features,
        verbose = TRUE
    )

    gex.anchors <- FindIntegrationAnchors(
        object.list = data.list,
        normalization.method = "SCT",
        anchor.features = anchor.features,
        reduction = "rpca",
        dims = 1:dims1, # previously did not specify so defaulted to 30, gave ok results so maybe lower # dimensions is acceptable
        reference=NULL
    )

    if(is.null(sample.tree)){
        message("Warning: No sample tree used")
        data.integrated <- IntegrateData(
                anchorset = gex.anchors,
                dims = 1:dims1,
                #features.to.integrate = anchor.features,    #Reduce(intersect,lapply(data.list,rownames)),
                normalization.method = "SCT"
        )
        }
    else{
        data.integrated <- IntegrateData(
                anchorset = gex.anchors,
                dims = 1:dims1,
                #features.to.integrate = Reduce(union, lapply(data.list, function(x) VariableFeatures(x, assay = "SCT"))),
                #features.to.integrate = anchor.features,    #Reduce(intersect,lapply(data.list,rownames)),
                normalization.method = "SCT",
                sample.tree = sample.tree
        )
    }
    message("GEX integration finished")

    # Run PCA on integrated assay to use in next step
    DefaultAssay(data.integrated) <- "integrated"
    data.integrated <- RunPCA(data.integrated, npcs=dims2)

    # also run dimensional reduction and clustering in case I want to check it later
    data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:dims2)
    data.integrated <- FindClusters(data.integrated, resolution = 2, algorithm = 4, reduction = "pca")
    data.integrated <- RunUMAP(data.integrated, reduction = "pca", seed.use = seed, dims = 1:dims2, reduction.key = "umapgex", reduction.name="UMAPgex")
    data.integrated$gex_clusters <- data.integrated$seurat_clusters

    message("GEX processing finished")    

    # Part 2: Processing ATAC data
    message("Start processing for ATAC data")

    # Binarize counts (now even simpler because there is a function in Signac package thankssss)
    for(sample in names(data.list)){DefaultAssay(data.list[[sample]]) <- "ATAC"}
    data.list <- lapply(data.list, BinarizeCounts, assay="ATAC")
    message("ATAC counts binarized")
    # I will skip re-calling of my peaks since I already called my peaks on pseudobulk data using macs2
    data.list <- lapply(data.list, FindTopFeatures, assay = "ATAC")
    data.list <- lapply(data.list, RunTFIDF, assay = "ATAC")
    data.list <- lapply(data.list, RunSVD, n = dims3)

    # Integrate across samples in ATAC 
    message("Dropping GEX assays from data.list!")
    for(sample in names(data.list)){data.list[['RNA']] <- NULL}
    for(sample in names(data.list)){data.list[['SCT']] <- NULL}

    message("Starting ATAC integration")
    # merge samples for ATAC data
    atac.merged <- merge(data.list[[1]],y=data.list[2:length(data.list)])

    DefaultAssay(atac.merged) <- "ATAC"
    atac.merged <- FindTopFeatures(atac.merged)
    atac.merged <- RunTFIDF(atac.merged)
    atac.merged <- RunSVD(atac.merged, n = dims3)

    atac.anchors <- FindIntegrationAnchors(
        object.list = data.list,
        anchor.features=Reduce(intersect,lapply(data.list,rownames)),
        reduction = "rlsi",
        dims = 2:dims3,
        reference = NULL
    )
    message("ATAC integration anchors set")

    if(is.null(sample.tree)){
        message("Warning: sample.tree not set")
        atac.integrated <- IntegrateEmbeddings(
            anchorset = atac.anchors,
            reductions = atac.merged[["lsi"]],
            new.reduction.name = "integrated_lsi",
            dims.to.integrate = 1:dims3
        )
    }
    else{
        atac.integrated <- IntegrateEmbeddings(
            anchorset = atac.anchors,
            reductions = atac.merged[["lsi"]],
            new.reduction.name = "integrated_lsi",
            dims.to.integrate = 1:dims3,
            sample.tree = sample.tree
        )
    }

    data.integrated@reductions$integrated_lsi <- atac.integrated@reductions$integrated_lsi
    data.integrated@reductions$lsi <- atac.merged@reductions$lsi
    message("ATAC integration finished")

    DefaultAssay(data.integrated) <- "ATAC"
    data.integrated <- RunUMAP(data.integrated, reduction = "integrated_lsi", seed.use=seed, dims = 2:dims3, reduction.key = "umapatac", reduction.name="UMAPatac")
    data.integrated <- FindNeighbors(data.integrated, reduction = "integrated_lsi", dims = 1:dims3)
    data.integrated <- FindClusters(data.integrated, resolution = 1, algorithm=1, reduction='integrated_lsi')
    data.integrated$atac_clusters <- data.integrated$seurat_clusters
    message("ATAC processing finished")

    # cleaning up env
    rm(data.list, atac.integrated, atac.merged)

    # Part 5: Integration across modalities
    message("Start multi-modal integration")
    data.integrated <- FindMultiModalNeighbors(
        object = data.integrated,
        reduction.list = list("pca", "integrated_lsi"),
        dims.list = list(1:dims2, 2:dims3),
        #modality.weight.name = "RNA.weight",
        verbose = TRUE
    )

    data.integrated <- RunUMAP(
        object = data.integrated,
        nn.name = "weighted.nn",
        assay = "integrated",
        seed.use = seed,
        verbose = TRUE,
        reduction.key="umapwnn",
        reduction.name="UMAPwnn"
    )

    data.integrated <- FindClusters(
        object=data.integrated,
        resolution = 1.0,
        algorithm=4,
        random.seed=404,
        graph.name = "wsnn"
    )

    data.integrated$wnn_clusters <- data.integrated$seurat_clusters
    message("Multi-modal integration finished.")
    message("Full integration complete")

    return(data.integrated)
}