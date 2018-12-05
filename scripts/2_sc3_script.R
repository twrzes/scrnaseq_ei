# Using SC3 to analyse the data

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

# Loading required libraries
library(SC3)
library(SingleCellExperiment)
library(scater)

# Create a SingleCellExperiment object
dataset_sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(gene_counts)), 
  colData = cell_annotations
)

# Calculating QC metrics
library(scater)
dataset_sce <- calculateQCMetrics(dataset_sce)

# Adding log2(x) + 1 into SCE object
dataset_sce <- normalize(dataset_sce)

# Define feature names in feature_symbol column
rowData(dataset_sce)$feature_symbol <- rownames(dataset_sce)
# Remove features with duplicated names
dataset_sce <- dataset_sce[!duplicated(rowData(dataset_sce)$feature_symbol), ]



# Calculate the clusters (based on k-means clustering)
sc_three <- sc3(dataset_sce, ks = 2:10, biology = TRUE, k_estimator = TRUE)

# See the estimated number of clusters
str(metadata(sc_three)$sc3$k_estimation)

# Interactive mode - opens a new window in the browser
sc3_interactive(sc_three)

