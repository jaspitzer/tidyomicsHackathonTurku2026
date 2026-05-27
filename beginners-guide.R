## Tidyomics Beginner's Guide
## Examples from https://github.com/tidyomics/.github using real Bioconductor datasets

# Install tidyomics suite if needed:
# BiocManager::install("tidyomics")


# ============================================================
# 1. tidySummarizedExperiment — bulk RNA-seq (airway)
#    Comparison: tidy API vs. base R
# ============================================================

library(airway)
library(tidySummarizedExperiment)

data(airway)
airway

# Tidy approach: mean counts per gene across treated samples
airway |>
  filter(dex == "trt") |>
  group_by(gene_id) |>
  summarize(mean_count = mean(counts))

# Equivalent base R
treated <- airway[, airway$dex == "trt"]
data.frame(
  gene_id    = rownames(treated),
  mean_count = rowMeans(assay(treated, "counts"))
)


# ============================================================
# 2. tidybulk — differential expression (airway)
# ============================================================

library(tidybulk)
library(ggplot2)

airway |>
  identify_abundant(factor_of_interest = dex) |>
  scale_abundance() |>
  test_differential_abundance(~ dex + cell) |>
  filter(is_de_edgeR) |>
  ggplot(aes(logFC, -log10(FDR), color = is_de_edgeR)) +
  geom_point() +
  theme_bw()

# PCA to visualise sample separation by treatment
airway |>
  identify_abundant(factor_of_interest = dex) |>
  scale_abundance() |>
  reduce_dimensions(method = "PCA") |>
  pivot_sample() |>
  ggplot(aes(PC1, PC2, color = dex, label = .sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel() +
  theme_bw()


# ============================================================
# 3. tidySingleCellExperiment — single-cell RNA-seq (PBMC 3k)
# ============================================================

library(TENxPBMCData)
library(scater)
library(tidySingleCellExperiment)
library(ggplot2)

pbmc <- TENxPBMCData("pbmc3k")

# QC, normalise, reduce dimensions
pbmc <- addPerCellQCMetrics(pbmc)
pbmc <- pbmc[, pbmc$sum > 500 & pbmc$detected > 200]
pbmc <- logNormCounts(pbmc)
pbmc <- runPCA(pbmc)
pbmc <- runUMAP(pbmc, dimred = "PCA")

# Tidy filtering and visualisation
pbmc |>
  filter(sum > 2000) |>
  ggplot(aes(UMAP1, UMAP2, color = log10(sum))) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(color = "log10(UMI)")


# ============================================================
# 4. plyranges — ATAC-seq peaks overlapping promoters
#    Data: fluentGenomics workflow package
# ============================================================

library(plyranges)
library(readr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load ATAC-seq peaks from fluentGenomics
peaks_file <- system.file("extdata", "ATAC_peak_metadata.txt.gz",
                          package = "fluentGenomics")
peaks_df <- read_tsv(peaks_file, col_types = c("cidciicdc"))

peaks_gr <- peaks_df |>
  as_granges(seqnames = chr) |>
  select(peak_id = gene_id) |>
  set_genome_info(genome = "GRCh38")

peaks_gr

# Build promoter windows from UCSC knownGene (GRCh38)
txdb        <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters_gr <- promoters(genes(txdb), upstream = 2000, downstream = 200)

# Tidy plyranges: count ATAC peaks per gene promoter
peaks_gr |>
  join_overlap_inner(promoters_gr) |>
  group_by(gene_id) |>
  summarize(n_peaks = n())
