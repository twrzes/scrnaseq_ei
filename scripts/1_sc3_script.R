# Using scater to visualise the data and SC3 to analyse the data

# Loading required libraries
library(SingleCellExperiment)
library(scater)
library(SC3)

# Changing the default R behaviour that all strings are treated as factors when reading the table
options(stringsAsFactors = FALSE)

# Loading gene counts matrix
gene_counts <- read.table("<TO_EDIT>", header = T, sep = "\t", check.names = FALSE)

# Putting gene names as row names
rownames(gene_counts) <- gene_counts[["gene_name"]]

# Removing first column
gene_counts <- subset(gene_counts, select = -gene_name)

# Loading cell annotation
cell_annotations <- read.table("<TO_EDIT>", header = T, sep = "\t", check.names = FALSE)

# Putting cell names as row names
rownames(cell_annotations) <- cell_annotations[["cell_name"]]

# Removing first column
cell_annotations <- subset(cell_annotations, select = -cell_name)

# Creating a SingleCellExperiment object
dataset_sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(gene_counts)),
  colData = cell_annotations
)

# Calculating QC metrics
dataset_sce <- calculateQCMetrics(dataset_sce)

# Adding log2(x) + 1 into SCE object
dataset_sce <- normalize(dataset_sce)

# Running PCA on log2 counts
dataset_sce <- runPCA(dataset_sce)

# Generating a PCA plot
plotPCA(dataset_sce, colour_by = "cell_type1")

# Generating the plot with amount of variance explained by each annotation column
plotExplanatoryVariables(dataset_sce)

# Run tSNE
dataset_sce <- runTSNE(dataset_sce, perplexity=<NUMBER_TO_EDIT>)

# Plot tSNE with perplexity as in previous step
plotTSNE(dataset_sce, colour_by = "cell_type1")

# Running PCA on to identify outliers
dataset_sce <- runPCA(dataset_sce, use_coldata = TRUE, detect_outliers = TRUE, selected_variables = c("pct_counts_in_top_100_features", "total_features_by_counts", "log10_total_features_by_counts"))

# Plot PCA with cell types and outliers
plotPCA(dataset_sce, colour_by = "cell_type1", size_by = "outlier")

# Define feature names in feature_symbol column
rowData(dataset_sce)$feature_symbol <- rownames(dataset_sce)
# Remove features with duplicated names
dataset_sce <- dataset_sce[!duplicated(rowData(dataset_sce)$feature_symbol), ]

# Preparing the data
sce <- sc3_prepare(dataset_sce)

# Calculating the optimal number of clusters
sce <- sc3_estimate_k(sce)

# See the estimated number of clusters
str(metadata(sce)$sc3$k_estimation)

# Calculate distance matrices corresponding to Euclidean, Pearson and Spearman distances
sce <- sc3_calc_dists(sce)

# Distance matrix transformation
sce <- sc3_calc_transfs(sce)

# K-means clustering
sce <- sc3_kmeans(sce, ks = <TO_EDIT>:<TO_EDIT>)

# Clustering solution
sce <- sc3_calc_consens(sce)

# Differentially expressed genes, cell markers and cell outliers
sce <- sc3_calc_biology(sce, ks = <TO_EDIT>:<TO_EDIT>)

# We can write the whole result tables by writing column and row data
write.table(as.data.frame(colData(sce)[,grep("sc3_", colnames(colData(sce)))]), file = "<TO_EDIT>", append = F, quote = F, sep = "\t", row.names = T, col.names = T)
write.table(as.data.frame(rowData(sce)[ , grep("sc3_", colnames(rowData(sce)))]), file = "<TO_EDIT>", append = F, quote = F, sep = "\t", row.names = T, col.names = T)

# Interactive mode - opens a new window in the browser
sc3_interactive(sce)
