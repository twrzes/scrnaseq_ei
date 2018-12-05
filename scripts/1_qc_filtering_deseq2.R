# Experiment QC and filtering, followed by DESEq2

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

# Note that row names in cell annotations are in the same order as column names in gene counts - this is absolutely critical as DESEq2 maps gene counts and sample annotations based on order, not on names

# Creating a SingleCellExperiment object
library(SingleCellExperiment)
dataset_sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(gene_counts)), 
  colData = cell_annotations
)

# Calculating QC metrics
library(scater)
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

# Filter out the outliers from the dataset
dataset_sce <- filter(dataset_sce, outlier == FALSE)

# Load DESeq2
library(DESeq2)

# Parallel processing - change number of cores
library(BiocParallel)
# This registers 4 cores in non-Windows computers (Mac, Linux)
register(MulticoreParam(<TO_EDIT>))
# This registers 4 cores in Windows-based computers
register(SnowParam(<TO_EDIT>))

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts(dataset_sce),
                              colData = colData(dataset_sce),
                              design = ~ <TO_EDIT>)

# Differential expression
# Depending of the dataset, this step may take up to few hours, therefore load the RDS file instead (it contains the saved R object with calculated differentially expressed genes)
dds <- DESeq(dds, test = "LRT", reduced = ~ <TO_EDIT>, sfType = "poscounts", useT = TRUE, minmu=1e-6, minReplicatesForReplace = Inf, parallel = TRUE)

# Loading a DESeq2 object
dds <- readRDS("<TO_EDIT>")

# Available results
resultsNames(dds)
# Calculate differential expression for chosen contrasts
res <- results(dds, parallel = TRUE, contrast = c("<TO_EDIT>", "<TO_EDIT>", "<TO_EDIT>"), cooksCutoff=FALSE)

# Writing result to a file
# Converting results to data frame
res_file <- as.data.frame(res)
# Ordering the data frame based on log2 fold-change
res_file <- res_file[order(res_file$log2FoldChange, decreasing = TRUE),]
# Adding row names as a separate column
res_file <- as.data.frame(cbind(gene_name = rownames(res_file), res_file))
# Writing the results to a specified file
write.table(x = res_file, file = "<TO_EDIT>", append = F, quote = F, sep = "\t", row.names = F, col.names = T)
