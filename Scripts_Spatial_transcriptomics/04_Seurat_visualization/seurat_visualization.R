#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
###########################################################
# Arguments
###########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript seurat_visualization.R input.rds output_prefix")
}

input_rds <- args[1]
output_prefix <- args[2]

cat("Input :", input_rds, "\n")
cat("Output:", output_prefix, "\n")

if (!file.exists(input_rds)) {
    stop(paste("Input file does not exist:", input_rds))
}

###########################################################
# Load Seurat object
###########################################################

seu <- readRDS(input_rds)

cat("Metadata columns:\n")
print(colnames(seu@meta.data))

###########################################################
# Detect cluster column
###########################################################

cluster_col <- NULL

if ("spatial_leiden" %in% colnames(seu@meta.data)) {
    cluster_col <- "spatial_leiden"
} else if ("leiden" %in% colnames(seu@meta.data)) {
    cluster_col <- "leiden"
}

###########################################################
# UMAP
###########################################################

if (!"umap" %in% names(seu@reductions)) {
    stop("UMAP reduction not found in Seurat object.")
}

if (!is.null(cluster_col)) {

    p1 <- DimPlot(
        seu,
        reduction = "umap",
        group.by = cluster_col,
        label = TRUE,
        repel = TRUE
    )
} else {


    p1 <- DimPlot(seu, reduction = "umap")
}

ggsave(
    filename = paste0(output_prefix, "_UMAP.png"),
    plot = p1,
    width = 7,
    height = 6,
    dpi = 300
)

###########################################################
# Spatial
###########################################################

if ("spatial" %in% names(seu@reductions)) {

    if (!is.null(cluster_col)) {

        p2 <- DimPlot(
            seu,
            reduction = "spatial",
            group.by = cluster_col,
            label = TRUE,
            pt.size = 1
        )
    } else {

        p2 <- DimPlot(seu, reduction = "spatial")
    }

    ggsave(
        filename = paste0(output_prefix, "_Spatial.png"),
        plot = p2,
        width = 8,
        height = 7,
        dpi = 300
    )

} else {
    warning("Spatial reduction not found. Spatial plot was skipped.")
}
###########################################################
# Marker gene identification
###########################################################

if (!is.null(cluster_col)) {

    cat("Using cluster column:", cluster_col, "\n")

    # Set identities
    Idents(seu) <- cluster_col

    # Find markers
    markers <- FindAllMarkers(
        seu,
        only.pos = TRUE,
        min.pct = 0.10,
        logfc.threshold = 0.25
    )

    # Save all markers
    write.csv(
        markers,
        paste0(output_prefix, "_all_markers.csv"),
        row.names = FALSE
    )

    # Top 5 markers per cluster
    top5 <- markers %>%
        group_by(cluster) %>%
        slice_max(avg_log2FC, n = 5, with_ties = FALSE)

    write.csv(
        top5,
        paste0(output_prefix, "_top5_markers.csv"),
        row.names = FALSE
    )

    #######################################################
    # FeaturePlots
    #######################################################

    plot_dir <- paste0(output_prefix, "_FeaturePlots")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

    for (i in seq_len(nrow(top5))) {

        gene <- top5$gene[i]
        clust <- top5$cluster[i]

        p <- FeaturePlot(
            seu,
            features = gene,
            reduction = "spatial"
        ) +
        ggtitle(paste("Cluster", clust, "-", gene))

        ggsave(
            filename = file.path(
                plot_dir,
                paste0("Cluster_", clust, "_", gene, ".png")
            ),
            plot = p,
            width = 6,
            height = 5,
            dpi = 300
        )
    }

    cat("Marker analysis completed.\n")

} else {

    warning("No Leiden clusters found. Marker analysis skipped.")

}
cat("Seurat visualization completed successfully.\n")
