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
) |>
  head()


# ============================================================
# 2. tidybulk — differential expression (airway)
# ============================================================

library(tidybulk)
library(ggplot2)

filt <- airway |>
  identify_abundant(factor_of_interest = dex) |>
  scale_abundance(method = "RLE")

de_airway <- filt |>
  test_differential_abundance(method = "DESeq2", .formula = ~cell + dex)

de_airway |>
  filter(padj < .01) |>
  ggplot(aes(log2FoldChange, -log10(padj))) +
  geom_point() +
  theme_bw() + 
  coord_cartesian(xlim=c(-4,4),ylim=c(0,30))

# PCA to visualise sample separation by treatment
filt |>
  reduce_dimensions(method = "PCA") |>
  pivot_sample() |>
  ggplot(aes(PC1, PC2, color = dex, label = .sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel() +
  theme_bw()

# ============================================================
# 3. plyranges — ATAC-seq peaks overlapping promoters
#    Data: fluentGenomics workflow package
# ============================================================

library(plyranges)
library(fluentGenomics)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

data(peaks)
peaks

# Build promoter windows from UCSC knownGene (GRCh38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters_gr <- promoters(genes(txdb), upstream = 2000, downstream = 200)

# Tidy plyranges: count ATAC peaks per gene promoter
peaks |>
  join_overlap_inner(promoters_gr) |>
  group_by(gene_id) |>
  summarize(n_peaks = n())
