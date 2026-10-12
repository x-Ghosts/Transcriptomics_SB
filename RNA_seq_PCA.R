setwd("/home/adam/Desktop/PhD_work/Systems_Biology/SystemsBiology_Brasile/RNA_Seq")

# Libraries :
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "org.Rn.eg.db"))

install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))



library(dplyr)
library(tidyverse)
library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify) #install.packages("ggfortify")
library(stats)
library(sva)


# Loading datasets:

transcripts <- read.csv("transcript_count_rnaseq.csv", header = T, stringsAsFactors = F, sep = ";")
gene <- read.csv("gene_count_rnaseq.csv", header = T, stringsAsFactors = F, sep = ";")

# Blood only:
blood_cols <- grepl("^RNA_B_", colnames(transcripts))
blood_cols_gene <- grepl("^RNA_B_", colnames(gene))
transcripts_blood <- transcripts[,c("tx", "gene_id", colnames(transcripts)[blood_cols])]
gene_count_blood <- gene[,c("gene_id", "gene_name", colnames(gene)[blood_cols])]



######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

# Original gene-count matrix
# ↓
# Correct numeric parsing
# ↓
# Metadata
# SampleID | IID | Breed | Treatment
# ↓
# Check pairing / missing samples
# ↓
# Filtering
# ↓
# Set factor levels
# ↓
# Construct DESeqDataSet + specify design
# ↓
# Normalization
# ↓
# Experimental design
# ↓
# DGE
# ↓
# Multiple testing
# ↓
# Pathway / enrichment analysis

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

######################################################################################################################################################################################
######################################################################################################################################################################################
"Data Processing"
######################################################################################################################################################################################
######################################################################################################################################################################################


# ============================================================
# Original gene-count matrix
# ============================================================

# Identify blood RNA-seq sample columns
blood_sample_cols <- grepl("^RNA_B_", colnames(gene_count_blood))
gene_annotation <- gene_count_blood[, c("gene_id", "gene_name")]
anyDuplicated(gene_annotation$gene_id) #Duplicate genes Checks
gene_count_matrix_raw <- as.matrix(gene_count_blood[, blood_sample_cols, drop = FALSE])
rownames(gene_count_matrix_raw) <- gene_annotation$gene_id

# Basic checks
dim(gene_count_matrix_raw) # 32092 x 32
head(gene_count_matrix_raw[, 1:5])

# ============================================================
# Correct numeric parsing
# ============================================================

gene_count_matrix <- matrix(
  as.numeric(gsub(",", ".", gene_count_matrix_raw, fixed = TRUE)  ),
  nrow = nrow(gene_count_matrix_raw),
  ncol = ncol(gene_count_matrix_raw),
  dimnames = dimnames(gene_count_matrix_raw)
)

# QC checks
typeof(gene_count_matrix) # Expected: "double"
sum(is.na(gene_count_matrix)) # Expected: 0
all(gene_count_matrix >= 0) # Expected: TRUE
range(gene_count_matrix)
quantile(
  gene_count_matrix,
  probs = c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 0.999, 1)
)

# gene <- read.csv("gene_count_rnaseq.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";", dec = ",")  ## FUTURE IMPORT

# ============================================================
# Metadata
# SampleID | IID | Breed | Treatment
# ============================================================

sample_ids <- colnames(gene_count_matrix)
sample_ids_clean <- sub("_+$", "", sample_ids)

sample_info <- strcapture(
  "^RNA_B_([AN])([0-9]+)_([^_]+)$",
  sample_ids_clean,
  proto = list(
    BreedCode = character(),
    AnimalNumber = character(),
    Treatment = character()
  )
)

meta <- data.frame(
  SampleID = sample_ids,
  IID = paste0(sample_info$BreedCode, sample_info$AnimalNumber),
  Breed = ifelse( sample_info$BreedCode == "A", "Angus",
                  ifelse( sample_info$BreedCode == "N", "Nellore", NA
                  )
  ),
  Treatment = sample_info$Treatment,
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$SampleID

# ============================================================
# Check pairing / missing samples
# ============================================================

# Number of samples for each IID x Treatment combination
pair_table <- table(
  IID = meta$IID,
  Treatment = meta$Treatment
)
nrow(pair_table) # 16 Animals = 8 Angus + 8 Nellore
pair_table

expected_treatments <- c("D4", "POST")
incomplete_IIDs <- rownames(pair_table)[!apply( pair_table[, expected_treatments, drop = FALSE], 1, function(x) all(x == 1))]
incomplete_IIDs

# ============================================================
# All variables datasets
# ============================================================
gene_count_matrix
gene_annotation
meta


######################################################################################################################################################################################
######################################################################################################################################################################################
"Data Processing - Filtering and DDS object creation"
######################################################################################################################################################################################
######################################################################################################################################################################################

# ============================================================
# Data filtering :
# ============================================================

smallestGroupSize <- 8
keep <- rowSums(gene_count_matrix >= 10) >= smallestGroupSize

gene_count_matrix_filtered <- gene_count_matrix[keep,,drop = FALSE]

cat("Genes before filtering:", nrow(gene_count_matrix), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "\n")
cat(
  "Percentage retained:",
  round(100 * mean(keep), 2),
  "%\n"
)

# ============================================================
# Factor levels
# ============================================================

head(meta)
meta$Treatment <- factor(
  meta$Treatment,
  levels = c("D4", "POST")
)

meta$Breed <- factor(
  meta$Breed,
  levels = c("Angus", "Nellore")
)

meta$IID <- factor(meta$IID)

# ============================================================
# Creating DDS object
# ============================================================

# Fractional Estimates:
fractional_values <- gene_count_matrix_filtered[gene_count_matrix_filtered %% 1 != 0]
head(fractional_values,30)

# Creation of Two signal datasets: ROUNDED Fractions, and DECIMAL commas removal
gene_count_matrix_rounded <- round(gene_count_matrix_filtered)
gene_count_matrix_dot_removed <- gene_count_matrix_filtered
fraction_index <- gene_count_matrix_dot_removed %% 1 != 0
gene_count_matrix_dot_removed[fraction_index] <- round(gene_count_matrix_dot_removed[fraction_index] * 1000)


# DDS Objects:
dds_rounded <- DESeqDataSetFromMatrix(countData = gene_count_matrix_rounded,
                                      colData = meta,
                                      design = ~ IID + Treatment
)
dds_dot_removed <-DESeqDataSetFromMatrix(countData = gene_count_matrix_dot_removed,
                                         colData = meta,
                                         design = ~ IID + Treatment
)
dds_rounded$Treatment <- relevel(dds_rounded$Treatment, ref = "POST")
dds_dot_removed$Treatment <- relevel(dds_dot_removed$Treatment, ref = "POST")


######################################################################################################################################################################################
######################################################################################################################################################################################
"PCA analysis"
######################################################################################################################################################################################
######################################################################################################################################################################################

shape_column <- "Treatment"
color_column <- "Breed"
label <- TRUE
label_size <- 2

# Variance Stabilizing Transformation
vsd_Round <- vst(dds_rounded, blind = TRUE)
vsd_dot <- vst(dds_dot_removed, blind = TRUE)

t_gene_count_matrix_rounded <- t(assay(vsd_Round))
t_gene_count_matrix_dot <- t(assay(vsd_dot))

meta_Round <- meta[
  rownames(t_gene_count_matrix_rounded), , drop = FALSE
]

meta_dot <- meta[
  rownames(t_gene_count_matrix_dot), , drop = FALSE
]

# stopifnot(
#   identical(
#     rownames(meta_Round),
#     rownames(t_gene_count_matrix_rounded)
#   ),
#   identical(
#     rownames(meta_dot),
#     rownames(t_gene_count_matrix_dot)
#   )
# )

pca_Round <- prcomp(
  t_gene_count_matrix_rounded,
  scale. = FALSE,
  center = TRUE
)

pca_dot <- prcomp(
  t_gene_count_matrix_dot,
  scale. = FALSE,
  center = TRUE
)

pca_plot_Round <- autoplot(
  pca_Round,
  data = meta_Round,
  shape = shape_column,
  colour = color_column,
  label = FALSE
) +
  geom_text_repel(
    aes(label = SampleID),
    size = label_size
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    title = "PCA - Rounded Counts (VST)",
    colour = "Breed",
    shape = "Treatment"
  )


pca_plot_dot <- autoplot(
  pca_dot,
  data = meta_dot,
  shape = shape_column,
  colour = color_column,
  label = FALSE
) +
  geom_text_repel(
    aes(label = SampleID),
    size = label_size
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    title = "PCA - Dot-Removed Counts (VST)",
    colour = "Breed",
    shape = "Treatment"
  )

pca_plot_Round
pca_plot_dot

ggsave(
  filename = "PCA_Rounded.pdf",
  plot = pca_plot_Round,
  width = 8,
  height = 6
)

ggsave(
  filename = "PCA_Dot_Removed.pdf",
  plot = pca_plot_dot,
  width = 8,
  height = 6
)
