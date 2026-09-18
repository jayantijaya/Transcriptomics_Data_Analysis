################################################################################
#Reference Based Cell Type Annotation
################################################################################

################################################################################
# Load libraries
################################################################################
library(ggplot2)
library(Seurat)
library(anndata)
library(dplyr)
library(Matrix)
library(future)
library(parallel)
library(AnnotationDbi)
library(org.Hs.eg.db)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
################################################################################
# Options
################################################################################

options(mc.cores = detectCores() - 1)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

################################################################################
# Command-line arguments
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage:\n",
    "Rscript celltype_annotation.R ",
    "<query.rds> <reference.h5ad> <sample_id> <output_dir>"
  )
}

query_rds     <- args[1]
reference_h5ad <- args[2]
sample_id     <- args[3]
out_dir       <- args[4]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# Output files
################################################################################

annotated_rds <- file.path(out_dir, paste0(sample_id, "_annotated.rds"))
markers_csv   <- file.path(out_dir, paste0(sample_id, "_markers.csv"))
prediction_csv <- file.path(out_dir, paste0(sample_id, "_prediction.csv"))
umap_png      <- file.path(out_dir, paste0(sample_id, "_celltype_umap.png"))




################################################################################
# 2. Define helper functions
################################################################################

## 2.1 Downsample reference by cell type
downsample_reference_by_celltype <- function(seurat_obj, target_n = 50000, 
                                             min_cells_keep = 500, 
                                             group_vars = c("cell_type"),
                                             seed = 123) {
  set.seed(seed)
  meta <- seurat_obj@meta.data
  meta$cell <- rownames(meta)
  
  # Compute how many to sample per group
  meta_grouped <- meta %>% 
    group_by(across(all_of(group_vars))) %>% 
    summarise(n = n(), .groups = "drop")
  total_cells <- sum(meta_grouped$n)
  meta_grouped <- meta_grouped %>%
    mutate(prop = n / total_cells,
           n_to_sample = ifelse(n < min_cells_keep, n, round(prop * target_n)))
  
  # Join back to meta to get n_to_sample per cell row
  meta <- left_join(meta, meta_grouped, by = group_vars)
  
  # Sample per group safely
  sampled_cells <- meta %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(~ {
      n_keep <- unique(.x$n_to_sample) # single value per group
      .x[sample(1:nrow(.x), min(n_keep, nrow(.x))), ]
    }) %>% 
    pull(cell)
  
  subset(seurat_obj, cells = sampled_cells)
}
## 2.2 Convert Ensembl IDs to Gene Symbols
map_ensembl_to_symbol_seurat <- function(
    srt,
    assay = "RNA",
    organism = c("human","mouse"),
    keep_ensembl_in_name = FALSE,     # TRUE: "ENSG00000141510_TP53", FALSE: "TP53" (aggregates duplicates)
    drop_unmapped = FALSE,           # Drop unmapped features if TRUE
    duplicate_agg = c("sum","mean"), # How to aggregate duplicates if keep_ensembl_in_name=FALSE
    strip_version = TRUE,
    verbose = TRUE
) {
  stopifnot(inherits(srt, "Seurat"))
  organism <- match.arg(organism)
  duplicate_agg <- match.arg(duplicate_agg)
  
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    if (organism == "human") { library(org.Hs.eg.db); orgdb <- org.Hs.eg.db }
    else { library(org.Mm.eg.db); orgdb <- org.Mm.eg.db }
    library(Matrix)
  })
  
  # 1) Get raw counts
  counts <- Seurat::GetAssayData(srt, assay = assay, layer = "counts")
  ens <- rownames(counts)
  if (strip_version) ens <- sub("\\..*$", "", ens)
  
  # 2) Map Ensembl -> SYMBOL
  syms <- AnnotationDbi::mapIds(orgdb, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  
  # 3) Handle unmapped
  new_names <- syms
  if (!drop_unmapped) {
    unm <- is.na(new_names) | new_names == ""
    new_names[unm] <- ens[unm]  # keep Ensembl if no mapping
  } else {
    keep_idx <- !(is.na(new_names) | new_names == "")
    counts <- counts[keep_idx, , drop = FALSE]
    ens <- ens[keep_idx]
    new_names <- new_names[keep_idx]
    syms <- syms[keep_idx]
  }
  
  # 4) Build final rownames
  if (keep_ensembl_in_name) {
    final_names <- paste(ens, new_names, sep = "_")
    rownames(counts) <- final_names
  } else {
    rownames(counts) <- new_names
    if (anyDuplicated(new_names)) {
      if (verbose) message("Aggregating duplicate symbols by ", duplicate_agg, ".")
      sym_fac <- factor(rownames(counts), levels = unique(rownames(counts)))
      G <- Matrix::sparse.model.matrix(~ sym_fac - 1)
      counts_agg <- Matrix::t(G) %*% counts
      rownames(counts_agg) <- levels(sym_fac)
      if (duplicate_agg == "mean") {
        sizes <- as.numeric(table(sym_fac))
        counts_agg <- counts_agg / sizes
      }
      counts <- counts_agg
    }
  }
  
  # 5) Replace counts in Seurat object
  new_assay <- Seurat::CreateAssayObject(counts = counts)
  srt[[assay]] <- new_assay     
  # 6) Build mapping table aligned with final rownames
  final_genes <- rownames(counts)
  if (keep_ensembl_in_name) {
    ensembl_final <- sub("_.*$", "", final_genes)
    symbol_final <- sub("^[^_]+_", "", final_genes)
  } else {
    ensembl_final <- NA_character_
    symbol_final <- final_genes
  }
  
  map_tbl <- data.frame(
    ensembl = ensembl_final,
    symbol = symbol_final,
    final = final_genes,
    stringsAsFactors = FALSE
  )
  
  attr(srt, "ensembl_symbol_map") <- map_tbl
  
  # 7) Report
  n_total <- length(syms)
  n_na <- sum(is.na(syms) | syms == "")
  n_dup <- sum(duplicated(new_names))
  if (verbose) {
    message(sprintf("Mapped: %.1f%% (%d/%d). Unmapped: %d.", 
                    100*(n_total - n_na)/n_total, n_total - n_na, n_total, n_na))
    if (!keep_ensembl_in_name) message("Duplicate symbols (pre-agg): ", n_dup)
  }
  
  return(srt)
}

################################################################################
# Load query
################################################################################

query <- readRDS(query_rds)

################################################################################
# Load reference (H5AD)
################################################################################

reference_ad <- read_h5ad(reference_h5ad)
counts <- t(reference_ad$X)
rownames(counts) <- rownames(reference_ad$var)
colnames(counts) <- rownames(reference_ad$obs)
reference <- CreateSeuratObject(
       counts = counts,
       meta.data = reference_ad$obs,min.features = 200,
  min.cells = 3
)
#reference <- CreateSeuratObject(
#    counts = t(reference_ad$X),
#    meta.data = reference_ad$obs,
#    min.cells = 3,
#    min.features = 200
#)
################################################################################
# 6. Downsample References
################################################################################

ref_down1 <- downsample_reference_by_celltype(
  reference,
  target_n = 10000,
  seed = 101
)


################################################################################
# 7. Convert Ensembl IDs to Gene Symbols
################################################################################

reference1 <- map_ensembl_to_symbol_seurat(
  ref_down1,
  organism = "human"
)

################################################################################
# 8. Select Reference
################################################################################

reference <- reference1

################################################################################
# 9. Normalize Reference
################################################################################

reference <- SCTransform(reference, verbose = FALSE)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference, npcs = 50)

reference <- RunUMAP(
  reference,
  reduction = "pca",
  dims = 1:30,
  return.model = TRUE
)

################################################################################
# 10. Normalize Query
################################################################################

query <- SCTransform(query, verbose = FALSE)

################################################################################
# 11. Find Transfer Anchors
################################################################################

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
)

################################################################################
# 12. Transfer Cell Type Labels
################################################################################

query <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = query,
  refdata = list(
    celltype = "cell_type"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

################################################################################
# 13. Visualize Predictions
################################################################################

p2<-DimPlot(
  query,
  reduction = "umap",
  group.by = "predicted.celltype"
)
filtered_umap_png <- file.path(
    out_dir,
    paste0(sample_id, "_celltype_umap.png")
)
ggsave(
  filename = filtered_umap_png,
  plot = p2,
  width = 8,
  height = 6,
  dpi = 300
)
################################################################################
# 14. Filter High-confidence Predictions
################################################################################

query_high <- subset(
  query,
  subset = predicted.celltype.score > 0.4
)

Idents(query_high) <- "predicted.celltype"

################################################################################
# 15. Remove Small Cell Types
################################################################################

cluster_counts <- table(Idents(query_high))

keep_clusters <- names(
  cluster_counts[cluster_counts >= 10]
)

query_filtered <- subset(
  query_high,
  idents = keep_clusters
)

################################################################################
# 16. Visualize Filtered Cells
################################################################################

p1<-DimPlot(
  query_filtered,
  reduction = "umap",
  group.by = "predicted.celltype"
)
filtered_umap_png <- file.path(
    out_dir,
    paste0(sample_id, "_celltype_filtered_umap.png")
)
ggsave(
  filename = filtered_umap_png,
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)
################################################################################
# 17. Marker Gene Detection
################################################################################

Idents(query_filtered) <- "predicted.celltype"

markers <- FindAllMarkers(query_filtered)

top_markers <-
  markers %>%
  group_by(cluster) %>%
  filter(
    avg_log2FC > 1,
    p_val_adj < 0.01
  ) %>%
  slice_head(n = 15) %>%
  ungroup()

################################################################################
# 18. Save Results
################################################################################
annotated_rds <- file.path(
    out_dir,
    paste0(sample_id, "_annotated.rds")
)

markers_csv <- file.path(
    out_dir,
    paste0(sample_id, "_markers.csv")
)

write.csv(
  top_markers,
  markers_csv,
  row.names = FALSE
)

saveRDS(
  query_filtered,
  annotated_rds
)
