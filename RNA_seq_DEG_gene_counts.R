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
"DEG analysis"
######################################################################################################################################################################################
######################################################################################################################################################################################

dds_Round <- DESeq(dds_rounded, parallel = TRUE) 
res_Round <- results(dds_Round)

### "- Rounded Part"
res_Round_df <- as.data.frame(res_Round)
res_Round_df$gene <- rownames(res_Round_df)
padj_significant_res_Round <- res_Round_df[res_Round_df$padj <= 0.05 & !is.na(res_Round_df$padj),]


# ### "- Dot Part"
# res_dot_df <- as.data.frame(res_dot)
dds_dot <- DESeq(dds_dot_removed, parallel = TRUE)
res_dot <- results(dds_dot)
# res_dot_df$gene <- rownames(res_dot_df)
# padj_significant_res_dot <- res_dot_df[res_dot_df$padj <= 0.05,]


######################################################################################################################################################################################
######################################################################################################################################################################################
"Volcano Plot"
######################################################################################################################################################################################
######################################################################################################################################################################################


fdr_cutoff <- 0.05
lfc_cutoff <- 1
n_labels <- 15

# Check available model coefficients
resultsNames(dds_Round)
coef_name <- "Treatment_D4_vs_POST" #stopifnot(coef_name %in% resultsNames(dds_Round))

res_Round_05 <- results(
  dds_Round,
  name = coef_name,
  alpha = fdr_cutoff
)
res_Round_05_df <-as.data.frame(res_Round_05)

res_Round_shrunk <- lfcShrink(
  dds_Round,
  coef = coef_name,
  res = res_Round_05,
  type = "apeglm"
)

volcano_df <- as.data.frame(res_Round_shrunk)
volcano_df$gene <- rownames(volcano_df)

# Exclude unavailable adjusted p-values
volcano_df <- volcano_df[
  !is.na(volcano_df$padj) &
    is.finite(volcano_df$log2FoldChange),
  ,
  drop = FALSE
]

# Calculate -log10(FDR)
volcano_df$neglog10padj <- -log10(
  pmax(volcano_df$padj, .Machine$double.xmin)
)


volcano_df$status <- "Not significant"
# Significant but small fold change
volcano_df$status[
  volcano_df$padj <= fdr_cutoff
] <- "Significant (small FC)"

# Higher expression at D4
volcano_df$status[
  volcano_df$padj <= fdr_cutoff &
    volcano_df$log2FoldChange >= lfc_cutoff
] <- "Up in D4"

# Higher expression at POST
volcano_df$status[
  volcano_df$padj <= fdr_cutoff &
    volcano_df$log2FoldChange <= -lfc_cutoff
] <- "Up in POST"

volcano_df$status <- factor(
  volcano_df$status,
  levels = c(
    "Up in D4",
    "Up in POST",
    "Significant (small FC)",
    "Not significant"
  )
)

volcano_labels <- volcano_df[
  volcano_df$status %in% c("Up in D4", "Up in POST"),
  ,
  drop = FALSE
]

volcano_labels <- volcano_labels[
  order(volcano_labels$padj),
  ,
  drop = FALSE
]
volcano_labels <- head(volcano_labels, n_labels)  #volcano_labels

volcano_plot_Round <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neglog10padj
  )
) +
  
  geom_point(
    aes(colour = status),
    alpha = 0.65,
    size = 1.5
  ) +
  
  scale_colour_manual(
    values = c(
      "Up in D4" = "#D73027",
      "Up in POST" = "#4575B4",
      "Significant (small FC)" = "#E6AB02",
      "Not significant" = "grey75"
    ),
    drop = FALSE
  ) +
  
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  geom_text_repel(
    data = volcano_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4
  ) +
  
  labs(
    title = "Volcano Plot - Angus / Nellore - Blood",
    subtitle = "D4 vs POST | Paired DESeq2 analysis",
    x = "Shrunken log2 Fold Change (D4 vs POST)",
    y = "-log10(adjusted p-value)",
    colour = "Gene classification"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    axis.line = element_line(colour = "black")
  )
print(volcano_plot_Round)

ggsave(
  filename = "Volcano_Rounded_D4_vs_POST.pdf",
  plot = volcano_plot_Round,
  width = 9,
  height = 7
)
