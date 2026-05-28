## Tidyomics Beginner's Guide — tidySeurat -- May 2026
## Single-cell RNA-seq: tidy interface over Seurat objects
## Data: 416B mouse cell line (Lun et al. 2017) via scRNAseq::LunSpikeInData()

# Install if needed:
# BiocManager::install(c("scRNAseq", "tidySeurat"))
# install.packages("Seurat")

# ============================================================
# 0. Build Seurat object from LunSpikeInData
# ============================================================

library(scRNAseq)
library(Seurat)
library(tidySeurat)
library(ggplot2)

# 416B mouse cell line; two conditions (wild type vs. oncostatin-M-induced)
sce <- LunSpikeInData(which = "416b")

# Features are Ensembl gene IDs; ERCC spike-in rows are separate
counts_mat <- counts(sce)
seurat_obj <- CreateSeuratObject(
  counts   = counts_mat[!grepl("^ERCC-", rownames(counts_mat)), ],
  meta.data = as.data.frame(colData(sce))
)

# ============================================================
# 1. Tidy inspection
# ============================================================

# tidySeurat: Seurat objects print and behave like tibbles
seurat_obj

# Select metadata columns of interest
seurat_obj |>
  select(.cell, phenotype, block, nCount_RNA, nFeature_RNA)

# Filter to one condition
seurat_obj |> filter(phenotype == "induced")

# Count cells per condition and plate
seurat_obj |>
  group_by(phenotype, block) |>
  summarize(n_cells = n())

# ============================================================
# 2. QC and filtering
# ============================================================

# Metadata columns work directly in aes() — no extra extraction needed
seurat_obj |>
  ggplot(aes(nCount_RNA, nFeature_RNA, color = block)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~phenotype) +
  scale_x_log10() + scale_y_log10() +
  theme_bw()

# Adjust thresholds after inspecting the QC plots above
seurat_filtered <- seurat_obj |>
  filter(nFeature_RNA > 1000, nCount_RNA > 1e5)

# ============================================================
# 3. Normalisation, variable features, and scaling
# ============================================================

# Standard Seurat pipeline — pipes work because each function returns
# the modified Seurat object
seurat_filtered <- seurat_filtered |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData()

# ============================================================
# 4. Dimensionality reduction and clustering
# ============================================================

seurat_filtered <- seurat_filtered |>
  RunPCA(npcs = 20) |>
  FindNeighbors(dims = 1:10) |>
  FindClusters(resolution = 0.4) |>
  RunUMAP(dims = 1:10)

# ============================================================
# 5. Tidy exploration of clustering results
# ============================================================

# group_by + summarize: cluster composition by condition
seurat_filtered |>
  group_by(seurat_clusters, phenotype) |>
  summarize(
    n_cells     = n(),
    mean_counts = mean(nCount_RNA)
  )

# Tidy ggplot2: tidySeurat exposes UMAP_1, UMAP_2 directly in aes()
seurat_filtered |>
  ggplot(aes(UMAP_1, UMAP_2, color = phenotype)) +
  geom_point(size = 1.5, alpha = 0.7) +
  theme_bw()

# Equivalent base Seurat
DimPlot(seurat_filtered, reduction = "umap", group.by = "phenotype")

# ============================================================
# 6. Join gene expression into the tidy frame
# ============================================================

# join_features() adds normalised expression columns for selected genes
top2 <- VariableFeatures(seurat_filtered)[1:2]

seurat_filtered |>
  join_features(features = top2, shape = "wide") |>
  ggplot(aes(.data[[top2[1]]], .data[[top2[2]]], color = phenotype)) +
  geom_point(alpha = 0.5) +
  labs(x = "Gene 1 (norm. expr.)", y = "Gene 2 (norm. expr.)") +
  theme_bw()

# Equivalent base Seurat
FeaturePlot(seurat_filtered, features = top2)
