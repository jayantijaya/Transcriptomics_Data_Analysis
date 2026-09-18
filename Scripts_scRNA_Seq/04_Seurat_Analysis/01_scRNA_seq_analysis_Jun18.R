##################################################################################################################
# This workflow performs downstream preprocessing and clustering analysis of CellBender-corrected 10x Genomics 
# Fixed RNA Profiling single-cell RNA sequencing data.
#   1. CellBender-corrected count matrix
#   2. Quality control filtering
#   3. Doublet detection using DoubletFinder
#   4. SCTransform normalization
#   5. PCA, UMAP, and tSNE dimensionality reduction
#   6. Cluster resolution optimization using clustree and silhouette score (for automated selection of best resolution)
#   7. Graph-based Leiden clustering using Seurat
#   8. Marker gene identification & Heatmap generation
# Date: 18/06/2026
##################################################################################################################

library(parallel)
library(Seurat)
library(scCustomize)
library(dplyr)
library(DoubletFinder)
library(Matrix)
library(gridExtra)
library(clustree)
library(ggplot2)
library(cluster)
set.seed(8)
options(mc.cores = detectCores() - 1)
args <- commandArgs(trailingOnly = TRUE)

input_h5  <- args[1]
sample_id <- args[2]
out_dir   <- args[3]

dir.create(out_dir,
           recursive = TRUE,
           showWarnings = FALSE)
################################ Load CellBender filtered h5 file ###############################################

cell_bender_mat <- Read_CellBender_h5_Mat(
    file_name = input_h5
)

################################ Create Seurat object ###########################################################

hf_S2 <- CreateSeuratObject(
    counts = cell_bender_mat,
    project = sample_id,
    min.features = 200,
    min.cells = 3
)

hf_S2[["percent.mt"]] <- PercentageFeatureSet(
    hf_S2,
    pattern = "^MT-"
)

hf_S2_filtered <- subset(
    hf_S2,
    subset = nFeature_RNA >= 200 &
             nFeature_RNA < 2500 &
             percent.mt < 15
)

################################ Violin plots ###################################################################

png("VlnPlot_before_filter.png",
    width = 2000,
    height = 800,
    res = 300)

VlnPlot(
    hf_S2,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
)

dev.off()

png("VlnPlot_after_filter.png",
    width = 2000,
    height = 800,
    res = 300)

VlnPlot(
    hf_S2_filtered,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
)

dev.off()

################################ Scatter plots ##################################################################

plot1_1 <- FeatureScatter(
    hf_S2_filtered,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    raster = FALSE
)

plot2_1 <- FeatureScatter(
    hf_S2_filtered,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    raster = FALSE
)

png("Scatter_after_filter.png",
    width = 2000,
    height = 1000,
    res = 300)

grid.arrange(plot1_1, plot2_1, ncol = 2)

dev.off()

plot1_2 <- FeatureScatter(
    hf_S2,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    raster = FALSE
)

plot2_2 <- FeatureScatter(
    hf_S2,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    raster = FALSE
)

png("Scatter_before_filter.png",
    width = 2000,
    height = 1000,
    res = 300)

grid.arrange(plot1_2, plot2_2, ncol = 2)

dev.off()

################################ DoubletFinder preprocessing ####################################################

hf_S2_filtered <- NormalizeData(
    hf_S2_filtered,
    verbose = TRUE
)

hf_S2_filtered <- ScaleData(hf_S2_filtered)

hf_S2_filtered <- FindVariableFeatures(
    hf_S2_filtered,
    selection.method = "vst",
    nfeatures = 2000
)

hf_S2_filtered <- RunPCA(
    hf_S2_filtered,
    verbose = TRUE
)

################################ Significant PC selection #######################################################

stdv <- hf_S2_filtered[["pca"]]@stdev
sum.stdv <- sum(stdv)

percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)

co1 <- which(cumulative > 90 & percent.stdv < 5)[1]

co2 <- sort(
    which(
        (percent.stdv[1:length(percent.stdv)-1] -
         percent.stdv[2:length(percent.stdv)]) > 0.1
    ),
    decreasing = TRUE
)[1] + 1

min.pc <- min(co1, co2)

################################ Pre-clustering #################################################################

hf_S2_filtered <- RunUMAP(
    hf_S2_filtered,
    dims = 1:min.pc
)

hf_S2_filtered <- FindNeighbors(
    hf_S2_filtered,
    dims = 1:min.pc
)

hf_S2_filtered <- FindClusters(
    hf_S2_filtered,
    resolution = 0.1
)

################################ pK identification ##############################################################

sweep.list <- paramSweep(
    hf_S2_filtered,
    PCs = 1:min.pc,
    num.cores = detectCores() - 1
)

sweep.stats <- summarizeSweep(sweep.list)

bcmvn <- find.pK(sweep.stats)

bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]

optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

################################ Homotypic doublet proportion ###################################################

annotations <- hf_S2_filtered@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)

doublet_rate <- 0.055

nExp.poi <- round(
    doublet_rate * nrow(hf_S2_filtered@meta.data)
)

nExp.poi.adj <- round(
    nExp.poi * (1 - homotypic.prop)
)

################################ DoubletFinder ##################################################################

hf_S2_filtered <- DoubletFinder::doubletFinder(
    seu = hf_S2_filtered,
    PCs = 1:min.pc,
    pK = optimal.pk,
    nExp = nExp.poi.adj,
    pN = 0.25
)

metadata <- hf_S2_filtered@meta.data
hf_S2_filtered@meta.data <- metadata

df_col <- grep(
    "DF.classifications",
    colnames(hf_S2_filtered@meta.data),
    value = TRUE
)

hf_S2_filtered$doublet_finder <- hf_S2_filtered@meta.data[[df_col]]

################################ Filter singlets ################################################################

hf_S2_filtered.singlets <- subset(
    hf_S2_filtered,
    doublet_finder == "Singlet"
)

################################ SCTransform ####################################################################

DefaultAssay(hf_S2_filtered.singlets) <- "RNA"

hf_S2_filtered.singlets <- SCTransform(
    hf_S2_filtered.singlets,
    vst.flavor = "v2",
    vars.to.regress = "percent.mt"
)

################################ PCA ###########################################################################

hf_S2_filtered.singlets <- RunPCA(
    hf_S2_filtered.singlets,
    npcs = 30,
    verbose = TRUE,
    features = VariableFeatures(hf_S2_filtered.singlets)
)

png("ElbowPlot.png",
    width = 1200,
    height = 1000,
    res = 300)

ElbowPlot(hf_S2_filtered.singlets)

dev.off()

################################ Variable feature plot ##########################################################


top10 <- head(
    VariableFeatures(hf_S2_filtered.singlets),
    10
)

plot1 <- VariableFeaturePlot(
    hf_S2_filtered.singlets
)

plot2 <- LabelPoints(
    plot = plot1,
    points = top10,
    repel = TRUE,
    xnudge = 0,
    ynudge = 0
)

png(
    "VariableFeatures.png",
    width = 1600,
    height = 1200,
    res = 300
)

print(plot2)

dev.off()

################################ Find neighbors #################################################################

hf_S2_filtered.singlets <- FindNeighbors(
    hf_S2_filtered.singlets,
    dims = 1:10
)

################################ Resolution optimization ########################################################

################################ Resolution optimization ########################################################

resolutions <- seq(0.1,1.0,0.01)

cluster_counts <- c()
sil_scores <- c()

pca_mat <- Embeddings(
    hf_S2_filtered.singlets,
    "pca"
)[,1:10]

for(res in resolutions){

    message("Testing resolution ", res)

    hf_S2_filtered.singlets <- FindClusters(
        hf_S2_filtered.singlets,
        resolution = res,
        algorithm = 4,
        verbose = FALSE
    )

    cluster_col <- paste0(
        "SCT_snn_res.",
        res
    )

    clusters <- hf_S2_filtered.singlets@meta.data[[cluster_col]]

    cluster_counts <- c(
        cluster_counts,
        length(unique(clusters))
    )

    sil <- silhouette(
        as.numeric(as.factor(clusters)),
        dist(pca_mat)
    )

    sil_scores <- c(
        sil_scores,
        mean(sil[,3])
    )
}

resolution_metrics <- data.frame(
    Resolution = resolutions,
    Clusters = cluster_counts,
    Silhouette = sil_scores
)

write.csv(
    resolution_metrics,
    file.path(
        out_dir,
        paste0(sample_id,
               "_resolution_metrics.csv")
    ),
    row.names = FALSE
)
p_sil <- ggplot(
    resolution_metrics,
    aes(
        Resolution,
        Silhouette
    )
) +
    geom_line() +
    geom_point() +
    theme_bw() +
    ggtitle(
        paste(
            sample_id,
            "Silhouette Scores"
        )
    )

ggsave(
    file.path(
        out_dir,
        paste0(sample_id,
               "_silhouette.png")
    ),
    p_sil,
    width = 6,
    height = 4
)

png(
    file.path(
        out_dir,
        paste0(
            sample_id,
            "_clustree.png"
        )
    ),
    width = 2000,
    height = 1600,
    res = 300
)
clust_tree <- clustree(
    hf_S2_filtered.singlets@meta.data,
    prefix = "SCT_snn_res."
)

print(clust_tree)

dev.off()

################################ Final clustering ###############################################################

best_res <- resolutions[
    which.max(sil_scores)
]

message(
    "Selected resolution = ",
    best_res
)

hf_S2_filtered.singlets$final_cluster <-hf_S2_filtered.singlets@meta.data[[paste0("SCT_snn_res.",best_res)]]

Idents(
    hf_S2_filtered.singlets
) <- "final_cluster"

################################ UMAP ###########################################################################

hf_S2_filtered.singlets <- RunUMAP(
    hf_S2_filtered.singlets,
    dims = 1:10
)

png("UMAP.png",
    width = 1600,
    height = 1400,
    res = 300)

DimPlot(
    hf_S2_filtered.singlets,
    reduction = "umap",
    group.by = "final_cluster",
    label = TRUE,
    pt.size = 0.5
)

dev.off()

################################ tSNE ###########################################################################

hf_S2_filtered.singlets <- RunTSNE(
    hf_S2_filtered.singlets,
    dims = 1:10
)

png("TSNE.png",
    width = 1600,
    height = 1400,
    res = 300)

DimPlot(
    hf_S2_filtered.singlets,
    reduction = "tsne",
    group.by = "final_cluster",
    label = TRUE,
    pt.size = 0.5
)

dev.off()

################################ Marker identification ##########################################################

hf_S2.markers <- FindAllMarkers(
    hf_S2_filtered.singlets,
    only.pos = TRUE
)

write.csv(
    hf_S2.markers,
    file.path(
        out_dir,
        paste0(
            sample_id,
            "_markers.csv"
        )
    ),
    row.names = FALSE
)

top_markers <- hf_S2.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01) %>%
    slice_head(n = 15) %>%
    ungroup()

################################ Heatmap ########################################################################

png("Heatmap.png",
    width = 3000,
    height = 2500,
    res = 300)

DoHeatmap(
    hf_S2_filtered.singlets,
    features = top_markers$gene
)

dev.off()
#################################Saved the selected resolution ##################################################
write.table(
    data.frame(
        Sample = sample_id,
        Resolution = best_res,
        Silhouette = max(sil_scores)
    ),
    file.path(
        out_dir,
        paste0(sample_id,
               "_selected_resolution.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
################################ Save Seurat object #############################################################

saveRDS(
    hf_S2_filtered.singlets,
    file.path(
        out_dir,
        paste0(
            sample_id,
            "_seurat.rds"
        )
    )
)

#######################Scrattch-Hicat Clustering ##################################################################
library(Seurat)
library(edgeR)
library(Matrix)
library(scrattch.hicat)

# Load Seurat object##
hf_S2_scrattch <- hf_S2_filtered.singlets

# Get SCT-derived highly variable genes####
hv.genes <- VariableFeatures(hf_S2_scrattch)

# Remove mitochondrial and ribosomal genes#

exclude.genes <- grep("^MT-|^RPL|^RPS", hv.genes, value = TRUE)
hv.genes <- hv.genes[!hv.genes %in% exclude.genes]

# Extract counts matrix####################

counts_mat <- GetAssayData(hf_S2_scrattch,layer = "counts")

# CPM normalization########################

norm.dat <- cpm(counts_mat)
norm.dat <- as.matrix(norm.dat)
norm.dat <- Matrix(norm.dat, sparse = TRUE)
# log transform
norm.dat@x <- log2(norm.dat@x + 1)

# Keep only selected features##############
#Used the marker genes from seurat markers#

top_markers.genes<-top_markers$gene
features.use <- intersect(unique(c(hv.genes, top_markers.genes)),rownames(norm.dat))
norm.dat <- norm.dat[features.use, ]

# Remove low-expression genes##############

norm.dat <- norm.dat[rowSums(norm.dat > 0) > 10,]

# Select cells############################

select.cells <- colnames(norm.dat)

# DE parameters###########################

strict.param <- de_param(de.score.th = 60)

# One-step clustering##################### 

onestep.result <- onestep_clust(norm.dat,select.cells = select.cells,dim.method = "pca",max.dim = 30, k.nn = 15,max.cl.size = 300,
                                de.param = strict.param,rm.eigen = NULL,verbose = TRUE)

# Add clusters to Seurat metadata ########

hf_S2_scrattch <- AddMetaData(hf_S2_scrattch,metadata = onestep.result$cl,col.name = "scrattch_cluster")

# Add clusters to Seurat metadata
#hf_S2_scrattch <- AddMetaData(hf_S2_scrattch,metadata = onestep.result$cl,col.name = "scrattch_cluster")

# Visualization ##########################

#DimPlot(hf_S2_scrattch, group.by = "scrattch_cluster",label = TRUE)

################################ Scrattch-Hicat Visualization ################################


# Set identities
Idents(hf_S2_scrattch) <- "scrattch_cluster"

################################ Save UMAP ###################################################

p_scrattch <- DimPlot(
    hf_S2_scrattch,
    reduction = "umap",
    group.by = "scrattch_cluster",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.5
) +
    ggtitle("Scrattch-Hicat Clusters")

ggsave(
    filename = file.path(
        out_dir,
        paste0(sample_id, "_Scrattch_HiCAT_UMAP.png")
    ),
    plot = p_scrattch,
    width = 6,
    height = 5,
    dpi = 300
)

################################ Optional TSNE ###############################################

if ("tsne" %in% Reductions(hf_S2_scrattch)) {

    p_tsne <- DimPlot(
        hf_S2_scrattch,
        reduction = "tsne",
        group.by = "scrattch_cluster",
        label = TRUE,
        repel = TRUE,
        pt.size = 0.5
    ) +
        ggtitle("Scrattch-HiCAT Clusters")

    ggsave(
        filename = file.path(
            out_dir,
            paste0(sample_id, "_Scrattch_HiCAT_TSNE.png")
        ),
        plot = p_tsne,
        width = 6,
        height = 5,
        dpi = 300
    )
}

################################ Save Seurat Object ##########################################

saveRDS(
    hf_S2_scrattch,
    file = file.path(
        out_dir,
        paste0(sample_id, "_scrattch_hicat_seurat.rds")
    )
)

################################ Save Cluster Assignments ####################################

write.csv(
    data.frame(
        Cell = colnames(hf_S2_scrattch),
        Scrattch_HiCAT_Cluster = hf_S2_scrattch$scrattch_cluster
    ),
    file = file.path(
        out_dir,
        paste0(sample_id, "_scrattch_clusters.csv")
    ),
    row.names = FALSE
)

