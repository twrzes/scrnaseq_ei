# Single cell RNA-Seq course
**Authors: Tomasz Wrzesinski, Wilfried Haerty**

**E-mail: [tomasz.wrzesinski@earlham.ac.uk](mailto:tomasz.wrzesinski@earlham.ac.uk)**

**Date: 06-07/12/2018**

**Venue: Earlham Institute**

## Software prerequisites
### Download and Installation of R
The R Language interpreter can be downloaded from: [https://cran.rstudio.com/](https://cran.rstudio.com/)

Specifically:
- **Windows**: [https://cran.rstudio.com/bin/windows/base/R-3.5.1-win.exe](https://cran.rstudio.com/bin/windows/base/R-3.5.1-win.exe)
- **Mac OSX**: [https://cran.rstudio.com/bin/macosx/R-3.5.1.pkg](https://cran.rstudio.com/bin/macosx/R-3.5.1.pkg)
- **Linux**: [https://cran.rstudio.com/bin/linux/](https://cran.rstudio.com/bin/linux/)

### Installation of RStudio
Download the appropriate software version according to your OS: [https://www.rstudio.com/products/rstudio/download/#download](https://www.rstudio.com/products/rstudio/download/#download)

Namely:
- **Windows**: [https://download1.rstudio.org/RStudio-1.1.463.exe](https://download1.rstudio.org/RStudio-1.1.463.exe)
- **Mac OSX**: [https://download1.rstudio.org/RStudio-1.1.463.dmg](https://download1.rstudio.org/RStudio-1.1.463.dmg)
- **Ubuntu**: [https://download1.rstudio.org/rstudio-1.1.463-amd64.deb](https://download1.rstudio.org/rstudio-1.1.463-amd64.deb)

### R libraries
In RStudio enter the commands - accept the installation of dependencies if necessary:
```R
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DESeq2", version = "3.8")
BiocManager::install("SC3", version = "3.8")
BiocManager::install("SingleCellExperiment", version = "3.8")
BiocManager::install("scater", version = "3.8")
BiocManager::install("edgeR", version = "3.8")
BiocManager::install("AneuFinder", version = "3.8")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("Rtsne")
install.packages("mvoutlier")
```

## More information
During preparation of these materials I used materials already available in the Internet:
- Single-cell RNA-Seq course organised by Wellcome Sanger Institute: [https://hemberg-lab.github.io/scRNA.seq.course](https://hemberg-lab.github.io/scRNA.seq.course)
- Single-cell RNA-Seq course organised by Harvard Chan Bioinformatics Core: [https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html](https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)
- An introduction to the scater package: [https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-intro.html](https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-intro.html)
- Data visualisation methods in scater: [https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-dataviz.html](https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-dataviz.html)
- Quality control with scater: [https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-qc.html](https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-qc.html)
- Analyzing RNA-Seq data with DESeq2: [https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- SC3 package manual: [https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html](https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html)
- K-means clustering explanation: [https://bigdata-madesimple.com/possibly-the-simplest-way-to-explain-k-means-algorithm/](https://bigdata-madesimple.com/possibly-the-simplest-way-to-explain-k-means-algorithm/)

## Data repository
You can find all the input files needed for the analysis at [https://drive.google.com/drive/folders/1K_RfWQVglVsjQ5eEKBYhXwrU5XRynDKi?usp=sharing](https://drive.google.com/drive/folders/1K_RfWQVglVsjQ5eEKBYhXwrU5XRynDKi?usp=sharing)

**Important: Download only the `dataset3` directory**

## Datasets
During our course we will be using three datasets, starting from `dataset3` - one focused on expression patterns in developing mouse embryo (`dataset1`) and two concentrated on different cell populations in human pancreas (`dataset2` and `dataset3`).

### Dataset 1 - Deng Q, _et al_., 2014
Link to the publication: [http://science.sciencemag.org/content/343/6167/193](http://science.sciencemag.org/content/343/6167/193)

**Abstract**:
Expression from both alleles is generally observed in analyses of diploid cell populations, but studies addressing allelic expression patterns genome-wide in single cells are lacking. Here, we present global analyses of allelic expression across individual cells of mouse preimplantation embryos of mixed background (CAST/EiJ × C57BL/6J). We discovered abundant (12 to 24%) monoallelic expression of autosomal genes and that expression of the two alleles occurs independently. The monoallelic expression appeared random and dynamic because there was considerable variation among closely related embryonic cells. Similar patterns of monoallelic expression were observed in mature cells. Our allelic expression analysis also demonstrates the de novo inactivation of the paternal X chromosome. We conclude that independent and stochastic allelic transcription generates abundant random monoallelic expression in the mammalian cell.

### Dataset 2 - Muraro MJ, _et al_., 2016
Link to the publication:
[https://www.sciencedirect.com/science/article/pii/S2405471216302927](https://www.sciencedirect.com/science/article/pii/S2405471216302927)

**Abstract**:
To understand organ function, it is important to have an inventory of its cell types and of their corresponding marker genes. This is a particularly challenging task for human tissues like the pancreas, because reliable markers are limited. Hence, transcriptome-wide studies are typically done on pooled islets of Langerhans, obscuring contributions from rare cell types and of potential subpopulations. To overcome this challenge, we developed an automated platform that uses FACS, robotics, and the CEL-Seq2 protocol to obtain the transcriptomes of thousands of single pancreatic cells from deceased organ donors, allowing in silico purification of all main pancreatic cell types. We identify cell type-specific transcription factors and a subpopulation of REG3A-positive acinar cells. We also show that CD24 and TM4SF4 expression can be used to sort live alpha and beta cells with high purity. This resource will be useful for developing a deeper understanding of pancreatic biology and pathophysiology of diabetes mellitus.

### Dataset 3 - Segerstolpe A, _et al_., 2016
Link to the publication:
[https://www.sciencedirect.com/science/article/pii/S1550413116304363](https://www.sciencedirect.com/science/article/pii/S1550413116304363)

**Abstract**:
Hormone-secreting cells within pancreatic islets of Langerhans play important roles in metabolic homeostasis and disease. However, their transcriptional characterization is still incomplete. Here, we sequenced the transcriptomes of thousands of human islet cells from healthy and type 2 diabetic donors. We could define specific genetic programs for each individual endocrine and exocrine cell type, even for rare δ, γ, ε, and stellate cells, and revealed subpopulations of α, β, and acinar cells. Intriguingly, δ cells expressed several important receptors, indicating an unrecognized importance of these cells in integrating paracrine and systemic metabolic signals. Genes previously associated with obesity or diabetes were found to correlate with BMI. Finally, comparing healthy and T2D transcriptomes in a cell-type resolved manner uncovered candidates for future functional studies. Altogether, our analyses demonstrate the utility of the generated single-cell gene expression resource.

## Input files
### Count matrix
The count matrix contains gene IDs as rows and sample names as columns. The filename of this file ends at `*_counts.txt`.
### Annotation file
The annotation file contains sample names as row and all available meta-data (conditions/batches) as columns. The filename of this file ends at `*_ann.txt`.

## 1. Starting the analysis
**Script name:** [https://github.com/twrzes/scrnaseq_ei/blob/master/scripts/1_sc3_script.R](https://github.com/twrzes/scrnaseq_ei/blob/master/scripts/1_sc3_script.R)
### 1.1. Loading the data to R
By default, all character strings in the table are being recognised as factors. To change this behaviour we need to type:
```R
options(stringsAsFactors = FALSE)
```
Next, we need to load the count matrix:
```R
gene_counts <- read.table("<TO_EDIT>", header = T, sep = "\t", check.names = FALSE)
```
**Important: Don't forget to change the path to the counts file in your script before running the line.**

#### TASK:
Look at the `gene_counts` object (you can do that either by clicking on the object under `Environment` tab in RStudio - this may take a lot of time - or by typing `head(gene_counts)` in RStudio which will give you an overview of columns in the file) - assuming that we are interested in differences in different cell types, can you see any potential confounding factors that may be problematic in the analysis of the dataset?

Because all the R libraries require to have gene names as row names, we need to change the first column to row names and remove the column from the data frame. This can be achieved by the following code:
```R
rownames(gene_counts) <- gene_counts[["gene_name"]]
gene_counts <- subset(gene_counts, select = -gene_name)
```
Next, we load the sample annotation file:
```R
cell_annotations <- read.table("<TO_EDIT>", header = T, sep = "\t", check.names = FALSE)
```
**Important: Don't forget to change the path to the annotation file in your script before running the line.**

#### TASK:
Look at the `cell_annotations` object under `Environment` tab in RStudio - assuming that we are interested in differences in different cell types, can you see any potential confounding factors that may be problematic in the analysis of the dataset?

Because all the libraries require to have cell names as row names, we need to change the first column to row names and remove the column from the data frame. This can be achieved by the following code:
```R
rownames(cell_annotations) <- cell_annotations[["cell_name"]]
cell_annotations <- subset(cell_annotations, select = -cell_name)
```

### 1.2. Creating a `SingleCellExperiment` object
**For more information you can refer to `SingleCellExperiment` documentation:** [http://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html](http://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)

The `SingleCellExperiment` class is a light-weight container for single-cell genomics data. It extends the `RangedSummarizedExperiment` class and follows similar conventions, i.e., rows should represent features (genes, transcripts, genomic regions) and columns should represent cells.

We can create an object with our data by using a `SingleCellExperiment` constructor:
```R
dataset_sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(gene_counts)),
  colData = cell_annotations
)
```
where:
- `assays`: list of matrices with unnormalised or normalised counts
- `colData`: DataFrame with cell annotations

#### TASK:
Which assays are available in `dataset_sce` object? (type `assays(dataset_sce)`). Which annotations are available in `dataset_sce` object? (type `colData(dataset_sce)`).

## 2. Expression QC
**For more information you can refer to `scater` documentation:** [http://bioconductor.org/packages/release/bioc/html/scater.html](http://bioconductor.org/packages/release/bioc/html/scater.html)

Once gene expression has been quantified as it is summarized as an **expression matrix** where each row corresponds to a gene (or transcript) and each column corresponds to a single cell. This matrix should be examined to remove poor quality cells which were not detected in either read QC or mapping QC steps. Failure to remove low quality cells at this stage may add technical noise which has the potential to obscure the biological signals of interest in the downstream analysis.

Since there is currently no standard method for performing scRNASeq the expected values for the various QC measures that will be presented here can vary substantially from experiment to experiment. Thus, to perform QC we will be looking for cells which are outliers with respect to the rest of the dataset rather than comparing to independent quality standards. Consequently, care should be taken when comparing quality metrics across datasets collected using different protocols.

All the quality metrics can be automatically calculated with `scater` package:
```R
dataset_sce <- calculateQCMetrics(dataset_sce)
```
We can view the calculated quality metrics with:
```R
colData(dataset_sce)
```

#### TASK:
Which metrics (i.e. columns) in the `colData(dataset_sce)` were added by `calculateQCMetrics` function?

Next, in order to perform a dimensionality reduction, we need to normalise the data - usually by doing a _log2_ transformation:
```R
dataset_sce <- normalize(dataset_sce)
```
### 2.1. Principal Component Analysis
The easiest way to overview the data is by transforming it using the principal component analysis and then visualize the first two principal components.

Principal component analysis (PCA) is a statistical procedure that uses a transformation to convert a set of observations into a set of values of linearly uncorrelated variables called principal components (PCs). The number of principal components is less than or equal to the number of original variables.

Mathematically, the PCs correspond to the eigenvectors of the covariance matrix. The eigenvectors are sorted by eigenvalue so that the first principal component accounts for as much of the variability in the data as possible, and each succeeding component in turn has the highest variance possible under the constraint that it is orthogonal to the preceding components (the figure below is taken from here).

![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/pca.png)

_Image credit: Vladimir Kiselev, Tallulah Andrews, Jennifer Westoby, Davis McCarthy, Maren Büttner and Martin Hemberg, Wellcome Sanger Institute, Cambridge_

Running the PCA can be achieved by running a line of code:
```R
dataset_sce <- runPCA(dataset_sce)
```
To plot the results in RStudio we can type:
```R
plotPCA(dataset_sce, colour_by = "cell_type1")
```
#### TASK:
Run the `runPCA` and `plotPCA` functions few times. Do you see any changes in the resulting plots?


We can investigate the relative importance of different explanatory factors with the  `plotExplanatoryVariables` function. We compute the R² for each factor in `colData(dataset_sce)` when fitting a linear model regressing expression values for each gene against that factor. This is best done on the log-expression values to reduce the effect of the mean on the variance - not a problem because we have already ran `normalize` first.

```R
plotExplanatoryVariables(dataset_sce)
```

#### TASK:
Which variable does explain the most of the variance in the data?

### 2.2. tSNE - t-Distributed Stochastic Neighbor Embedding
An alternative to PCA for visualizing scRNA-Seq data is a tSNE plot. [It has been introduced by van der Maaten and Hinton in 2008](http://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf). tSNE (t-Distributed Stochastic Neighbor Embedding) combines dimensionality reduction (e.g. PCA) although, in contrast with PCA, tSNE is a stochastic algorithm which means running the method multiple times on the same dataset will result in different plots.

Before diving in: if you haven’t encountered t-SNE before, here’s what you need to know about the math behind it. The goal is to take a set of points in a high-dimensional space and find a faithful representation of those points in a lower-dimensional space, typically the 2D plane. The algorithm is non-linear and adapts to the underlying data, performing different transformations on different regions. Those differences can be a major source of confusion.

A second feature of t-SNE is a tuneable parameter, _"perplexity”_, which says (loosely) how to balance attention between local and global aspects of your data. The parameter is, in a sense, a guess about the number of close neighbors each point has. The perplexity value has a complex effect on the resulting pictures. The original paper says, _“The performance of SNE is fairly robust to changes in the perplexity, and typical values are between 5 and 50.”_ But the story is more nuanced than that. Getting the most from t-SNE may mean analyzing multiple plots with different perplexities.

That’s not the end of the complications. The t-SNE algorithm doesn’t always produce similar output on successive runs, for example, and there are additional hyperparameters related to the optimization process.

**A great overview of tSNE pitfalls:** [https://distill.pub/2016/misread-tsne/](https://distill.pub/2016/misread-tsne/)

To prepare tSNE plot on our data, we can run:
```R
dataset_sce <- runTSNE(dataset_sce, perplexity=<NUMBER_TO_EDIT>)
plotTSNE(dataset_sce, colour_by = "cell_type1")
```
**Important: Don't forget to change the perplexity parameter to any number you have in mind.**

#### TASK:
- Run the `runTSNE` and `plotTSNE` functions with different perplexity numbers. How do different perplexity numbers affect your tSNE plot?
- Run the `runTSNE` and `plotTSNE` functions with chosen perplexity number few times. Do you obtain the same or different results (i.e. plots) with the same perplexity numbers?

### 2.3. Identification of outliers - automatic PCA-based method
Another option available in `scater` is to conduct PCA on a set of QC metrics and then use automatic outlier detection to identify potentially problematic cells.

By default, the following metrics are used for PCA-based outlier detection:
- pct_counts_top_100_features
- total_features
- pct_counts_feature_controls
- n_detected_feature_controls
- log10_counts_endogenous_features
- log10_counts_feature_controls

`scater` first creates a matrix where the rows represent cells and the columns represent the different QC metrics. Here, the PCA plot provides a 2D representation of cells ordered by their quality metrics. The outliers are then detected using methods from the `mvoutlier` package.

Since we do not have any ERCC spike-ins or mitochondrial genes in the dataset, we need to tell the function which variables we want to use for outlier detection. This can be done with `selected_variables` option in `runPCA` function:
```R
dataset_sce <- runPCA(dataset_sce,
  use_coldata = TRUE,
  detect_outliers = TRUE,
  selected_variables = c("pct_counts_in_top_100_features", "total_features_by_counts", "log10_total_features_by_counts"))
```

Plotting the PCA plot with outliers depicted as different point sizes:
```R
plotPCA(dataset_sce, colour_by = "cell_type1", size_by = "outlier")
```

### TASK:
Do the data contain any outliers identified with automatic PCA-based method?

## 3. Clustering approaches and problems
Once we have normalized the data we can carry out analyses that are relevant to the biological questions at hand. The exact nature of the analysis depends on the dataset. Nevertheless, there are a few aspects that are useful in a wide range of contexts, one of them being the clustering of scRNA-seq data.

One of the most promising applications of scRNA-seq is _de novo_ discovery and annotation of cell-types based on transcription profiles. Computationally, this is a hard problem as it amounts to **unsupervised clustering**. That is, we need to identify groups of cells based on the similarities of the transcriptomes without any prior knowledge of the labels. Moreover, in most situations we do not even know the number of clusters _a priori_. The problem is made even more challenging due to the high level of noise (both technical and biological) and the large number of dimensions (i.e. genes).

Unsupervised clustering is useful in many different applications and it has been widely studied in machine learning. Some of the most popular approaches are **hierarchical clustering** and **k-means clustering**

### 3.1. Hierarchical clustering
In hierarchical clustering, one can use either a bottom-up or a top-down approach. In the former case, each cell is initially assigned to its own cluster and pairs of clusters are subsequently merged to create a hierarchy:

Raw data:

![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/hierarchical_clustering1.png)

Creating a hierarchical clustering graph:

![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/hierarchical_clustering2.png)

_Image credits: Vladimir Kiselev, Tallulah Andrews, Jennifer Westoby, Davis McCarthy, Maren Büttner and Martin Hemberg, Wellcome Sanger Institute, Cambridge_

With a top-down strategy, one instead starts with all observations in one cluster and then recursively split each cluster to form a hierarchy. One of the advantages of this strategy is that the method is deterministic.

### 3.2. K-means clustering
In k-means clustering, the goal is to partition N cells into k different clusters. In an iterative manner, cluster centers are assigned and each cell is assigned to its nearest cluster:

![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/k-means.png)

_Image credit: Vladimir Kiselev, Tallulah Andrews, Jennifer Westoby, Davis McCarthy, Maren Büttner and Martin Hemberg, Wellcome Sanger Institute, Cambridge_


To run a k-means algorithm, you have to randomly initialize three points called the cluster centroids. In the example we have three cluster centroids, because we want to group a data into three clusters. K-means is an iterative algorithm and it does two steps: 1. Cluster assignment step 2. Move centroid step.

In Cluster assignment step, the algorithm goes through each of the data points and depending on which cluster is closer, whether the red cluster centroid or the blue cluster centroid or the green. It assigns the data points to one of the three cluster centroids.

In move centroid step, K-means moves the centroids to the average of the points in a cluster. In other words, the algorithm calculates the average of all the points in a cluster and moves the centroid to that average location.

### 3.3. Problems with the clustering methods
- What is the number of clusters k?
- What is a cell type?
- Scalability: in the last few years the number of cells in scRNA-seq experiments has grown by several orders of magnitude from hundreds to millions.
- Tools are not user-friendly

## 4. SC3 library - clustering example
**For more information you can refer to `SC3` documentation:** [http://bioconductor.org/packages/release/bioc/html/SC3.html](http://bioconductor.org/packages/release/bioc/html/SC3.html)

Single-Cell Consensus Clustering (`SC3`) is a tool for unsupervised clustering of scRNA-seq data. `SC3` achieves high accuracy and robustness by consistently integrating different clustering solutions through a consensus approach. An interactive graphical implementation makes `SC3` accessible to a wide audience of users. In addition, `SC3` also aids biological interpretation by identifying marker genes, differentially expressed genes and outlier cells. A manuscript describing `SC3` in details is published in [Nature Methods](http://dx.doi.org/10.1038/nmeth.4236).

![image](https://hemberg-lab.github.io/scRNA.seq.course/figures/sc3.png)

_Image credit: Vladimir Kiselev, Tallulah Andrews, Jennifer Westoby, Davis McCarthy, Maren Büttner and Martin Hemberg, Wellcome Sanger Institute, Cambridge_

### 4.1. Creating a column in rowData with feature names
In order to have feature names included in all the analyses (i.e. tables and plots), we need to put their names in the `feature_symbol` column in row data:
```R
rowData(sce)$feature_symbol <- rownames(sce)
```
We also need to remove some duplicated feature names (if there are any):
```R
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
```

### 4.2. Preparing the data
We start with `sc3_prepare`. This method prepares an object of `sce` class for `SC3` clustering. This method also defines all parameters needed for clustering and stores them in the `sc3` slot:
```R
sce <- sc3_prepare(sce)
```

### 4.3. Calculating the optimal number of clusters
When the `sce` object is prepared for clustering, `SC3` can also estimate the optimal number of clusters `k` in the dataset. `SC3` utilizes the Tracy-Widom theory on random matrices to estimate `k`:
```R
sce <- sc3_estimate_k(sce)
```
We can print the predicted number of clusters with:
```R
str(metadata(sce)$sc3$k_estimation)
```

#### TASK:
What is 'the best' number of clusters identified by the algorithm?

### 4.4. Calculating distances
Now we are ready to perform the clustering itself. First `SC3` calculates distances between the cells (i.e. Euclidean, Pearson and Spearman distances):
```R
sce <- sc3_calc_dists(sce)
```

### 4.5. Distance matrix transformation
Next the distance matrices are transformed using PCA and graph Laplacian:
```R
sce <- sc3_calc_transfs(sce)
```

### 4.6. K-means clustering
K-means should then be performed on the transformed distance matrices:
```R
sce <- sc3_kmeans(sce, ks = <TO_EDIT>:<TO_EDIT>)
```
where:
- ks: number of clusters; you can specify a range of `k`, e.g. if you want to specify between 2 to 6 clusters, the command will look like: `sce <- sc3_kmeans(sce, ks = 2:6)`

**Important information: Don't forget to change the number of clusters for k-means clustering**

#### TASK:
Calculate K-means with desired number of clusters (the number should overlap predicted number of clusters).

### 4.7. Clustering solution
In this step `SC3` will provide you with a clustering solution. When calculating consensus for each value of `k`, `SC3` averages the clustering results of k-means using a consensus approach.
```R
sce <- sc3_calc_consens(sce)
```

### 4.8. Differentially expressed genes, cell markers and cell outliers
`SC3` can also calculates DE genes, marker genes and cell outliers based on the calculated consensus clusterings.
```R
sce <- sc3_calc_biology(sce, ks = <TO_EDIT>:<TO_EDIT>)
```
**Important information: Don't forget to change the number of clusters which is equal to the number of clusters in k-means clustering - `ks` parameter**.

### 4.9 Writing result tables by writing column and row data
We can write the whole clustering solutions to the text files:
```R
write.table(as.data.frame(colData(sce)[,grep("sc3_", colnames(colData(sce)))]), file = "<TO_EDIT>", append = F, quote = F, sep = "\t", row.names = T, col.names = T)
write.table(as.data.frame(rowData(sce)[ , grep("sc3_", colnames(rowData(sce)))]), file = "<TO_EDIT>", append = F, quote = F, sep = "\t", row.names = T, col.names = T)
```
**Important information: Don't forget to change the names of the output files**.

#### TASK:
Open the colData and rowData files either in text editor (e.g. Notepad) or in Excel. What types of information do you have in each of these files?

### 4.10. Interactive visualisation
To quickly and easily explore the `SC3` solutions using an interactive Shiny application (opens in your default system Internet browser) use the following method:
```R
sc3_interactive(sce)
```

### TASK:
Play with different `k` parameters in your data. Does the best `k` parameter is the one chosen automatically?

#### 4.10.1. Consensus matrix
The consensus matrix is a N by N matrix, where N is the number of cells in the input dataset. It represents similarity between the cells based on the averaging of clustering results from all combinations of clustering parameters. Similarity 0 (blue) means that the two cells are always assigned to different clusters. In contrast, similarity 1 (red) means that the two cells are always assigned to the same cluster. The consensus matrix is clustered by hierarchical clustering and has a diagonal-block structure. Intuitively, the perfect clustering is achieved when all diagonal blocks are completely red and all off-diagonal elements are completely blue.

#### 4.10.2. Silhouette plot
A silhouette is a quantitative measure of the diagonality of the consensus matrix. An average silhouette width (shown at the bottom left of the silhouette plot) varies from 0 to 1, where 1 represents a perfectly block-diagonal consensus matrix and 0 represents a situation where there is no block-diagonal structure. The best clustering is achieved when the average silhouette width is close to 1.

#### 4.10.3. Expression matrix
The expression panel represents the original input expression matrix (cells in columns and genes in rows) after cell and gene filters. Genes are clustered by kmeans with k = 100 (dendrogram on the left) and the heatmap represents the expression levels of the gene cluster centers after log2-scaling.

#### 4.10.4. Cluster stability
Stability index shows how stable each cluster is accross the selected range of ks. The stability index varies between 0 and 1, where 1 means that the same cluster appears in every solution for different k.

#### 4.10.5. Differentially expressed genes
Differential expression is calculated using the non-parametric Kruskal-Wallis test. A significant p-value indicates that gene expression in at least one cluster stochastically dominates one other cluster. SC3 provides a list of all differentially expressed genes with adjusted p-values < 0.01 and plots gene expression profiles of the 50 genes with the lowest p-values. Note that the calculation of differential expression after clustering can introduce a bias in the distribution of p-values, and thus we advise to use the p-values for ranking the genes only.

#### 4.10.6. Marker genes
To find marker genes, for each gene a binary classifier is constructed based on the mean cluster expression values. The classifier prediction is then calculated using the gene expression ranks. The area under the receiver operating characteristic (ROC) curve is used to quantify the accuracy of the prediction. A p-value is assigned to each gene by using the Wilcoxon signed rank test. By default the genes with the area under the ROC curve (AUROC) > 0.85 and with the p-value < 0.01 are selected and the top 10 marker genes of each cluster are visualized in this heatmap.

## 5. Differential expression testing with DESeq2
**For more information you can refer to `DESeq2` documentation:** [http://bioconductor.org/packages/release/bioc/html/DESeq2.html](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)

Please refer to the [**1. Starting the analysis**](https://github.com/twrzes/scrnaseq_ei#1-starting-the-analysis) and [**2. Expression QC**](https://github.com/twrzes/scrnaseq_ei#2-expression-qc) sections of this tutorial.

**Script name:** [https://github.com/twrzes/scrnaseq_ei/blob/master/scripts/2_deseq2_script.R](https://github.com/twrzes/scrnaseq_ei/blob/master/scripts/2_deseq2_script.R)

In order to run DESeq2 we have to remove outliers from the subsequent analysis:
```R
dataset_sce <- filter(dataset_sce, outlier == FALSE)
```

The final step in our analysis workflow is fitting the raw counts to the negative binomial (NB) model and performing the statistical test for differentially expressed genes. In this step we essentially want to determine whether the mean expression levels of different sample groups are significantly different.

![image](https://hbctraining.github.io/DGE_workshop/img/de_theory.png)

_Image credit: Paul Pavlidis, UBC_

Differential expression analysis with `DESeq2` involves multiple steps as displayed in the flowchart below in blue. Briefly, `DESeq2` will model the raw counts, using normalization factors (size factors) to account for differences in library depth. Then, it will estimate the gene-wise dispersions and shrink these estimates to generate more accurate estimates of dispersion to model the counts. Finally, `DESeq2` will fit the negative binomial model and perform hypothesis testing using the Wald test or Likelihood Ratio Test (which we will be using in our case).

![image](https://hbctraining.github.io/DGE_workshop/img/DESeq2_workflow_2018.png)

_Image credit: Members of the teaching team at the Harvard Chan Bioinformatics Core (HBC)_


**If you want to know more about different steps included in this chart, you can read about them either in DESEq2 documentation or in [this tutorial](https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)**

Prior to performing the differential expression analysis, it is a good idea to know what **sources of variation** are present in your data, either by exploration during the QC and/or prior knowledge. This includes batch effects. Once you know the major sources of variation, you can remove them prior to analysis or control for them in the statistical model by including them in your **design formula**.

### 5.1. Design formula
A design formula tells the statistical software the known sources of variation to control for, as well as, the factor of interest to test for during differential expression testing. For example, if you know that sex is a significant source of variation in your data, then `sex` should be included in your model. **The design formula should have all of the factors in your metadata that account for major sources of variation in your data. The last factor entered in the formula should be the condition of interest**.

For example, suppose you have the following metadata:

![image](https://hbctraining.github.io/DGE_workshop/img/meta_example.png)

If you want to examine the expression differences between treatments, and you know that major sources of variation include sex and age, then your design formula would be:
```R
design <- ~ sex + age + treatment
```

The advantage of `DESeq2` is that we can run the whole analysis in just **2 lines of code!**. The disadvantage is that the analysis takes a considerable amount of time to finish.

### 5.2. Creating a DESeqDataSet object
First, we need to tell `DESeq2` how many cores we want to use for our analysis. If you have UNIX-based laptop (Linux, Mac), then run the following command:
```R
register(MulticoreParam(<TO_EDIT>))
```
If you have a Windows-based laptop, type:
```R
register(SnowParam(<TO_EDIT>))
```
**Important information: Don't forget to edit the number of cores equal to the number of cores you want to use for the analysis.**

We create a DESEq2 object with this line of code:
```R
dds <- DESeqDataSetFromMatrix(countData = counts(dataset_sce),
                              colData = colData(dataset_sce),
                              design = ~ <TO_EDIT>)
```
**Important information: Don't forget to edit the design formula based on the annotation columns available. Reminder: you can see these columns by running the command `colData(dataset_sce)`.**

#### TASK:
Which design formula would you use to see differentially expressed genes between different cell types between different disease conditions?


### 5.3 Running a differential expression with Likelihood Ratio Test (LRT)
`DESeq2` offers two kinds of hypothesis tests: the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero, and the likelihood ratio test (LRT). The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model.

The likelihood ratio test can be performed by specifying `test="LRT"` when using the `DESeq` function, and providing a reduced design formula, e.g. one in which a number of terms from `design(dds)` are removed. The degrees of freedom for the test is obtained from the difference between the number of parameters in the two models.

In the case presented above, the reduced model would look like:
```R
reduced <- ~ sex + age
```
We can run the following code to make a Likelihood Ratio Test between full model (i.e. with variable of interest) and reduced model (i.e. with unwanted sources of variation):
```R
dds <- DESeq(dds, test = "LRT", reduced = ~ <TO_EDIT>, sfType = "poscounts",
  useT = TRUE, minmu=1e-6, minReplicatesForReplace = Inf, parallel = TRUE)
```
where:
- `test = "LRT"`: Denotes the test used (i.e. likelihood ratio test)
- `reduced = `: Specifies the reduced model (i.e. model without variable of interest)
- `sfType = "poscounts"`: Method for size factor estimation. The ``"poscounts"`` estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts. This evolved out of use cases with Paul McMurdie's `phyloseq` package for metagenomic samples
- `useT = TRUE`: Whether to use a t-distribution as a null distribution, for significance testing
- `minmu=1e-6`: Lower bound on the estimated count for fitting gene-wise dispersion
- `minReplicatesForReplace = Inf`: The minimum number of replicates required in order to use replaceOutliers on a sample. Set to `Inf` in order to never replace outliers.
- `parallel = TRUE`: Parallel execution with `BiocParallel`.

**Important information: Don't forget to edit the formula for reduced model.**

#### TASK:
Which reduced model would you use in our case?

### 5.4. Loading a DESEq2 object into environment
The above command will take between 10 minutes and an hour, depending on the dataset. To load a pre-computed `DESeqDataSet` object we can run a following command:

```R
dds <- readRDS("<TO_EDIT>")
```
**Important information: Don't forget to change the path to the `*.rds` file.**

### 5.5. Calculating the results
In order to obtain out results we need to use a function `results`. It extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values.

#### 5.5.1. Contrasts
A contrast is a linear combination of estimated log2 fold changes, which can be used to test if differences between groups are equal to zero. The simplest use case for contrasts is an experimental design containing a factor with three levels, say A, B and C. Contrasts enable the user to generate results for all 3 possible differences: log2 fold change of B vs A, of C vs A, and of C vs B. The contrast argument of results function is used to extract test results of log2 fold changes of interest, for example:
```R
results(dds, contrast=c("condition","C","B"))
```
In this case, the log2 fold-change will be calculated for `C` vs. `B` conditions, specified in the `condition` column in annotation data.

To look what contrasts are available, you can run a following command:
```R
resultsNames(dds)
```
Result names are presented in a following format, e.g.:
```R
> resultsNames(dds)
 [1] "Intercept"                                   "sex_male_vs_female"                         
 [3] "age"                                         "cell_type1_alpha_vs_acinar"                 
 [5] "cell_type1_beta_vs_acinar"                   "cell_type1_co.expression_vs_acinar"         
 [7] "cell_type1_delta_vs_acinar"                  "cell_type1_ductal_vs_acinar"                
 [9] "cell_type1_endothelial_vs_acinar"            "cell_type1_epsilon_vs_acinar"               
[11] "cell_type1_gamma_vs_acinar"                  "cell_type1_mast_vs_acinar"                  
[13] "cell_type1_MHC.class.II_vs_acinar"           "cell_type1_not.applicable_vs_acinar"        
[15] "cell_type1_PSC_vs_acinar"                    "cell_type1_unclassified_vs_acinar"          
[17] "cell_type1_unclassified.endocrine_vs_acinar"
```
In order to compare e.g. `epsilon` samples with `PSC` samples in `cell_type1` column (i.e. log2 fold-change with `epsilon` as numerator and `PSC` as denominator) in the above example, you can specify following contrast:
```R
results(dds, contrast=c("cell_type1","epsilon","PSC"))
```


#### 5.5.2. Obtaining the result matrix
Coming back to our dataset, we can obtain the result table by running the following command:
```R
res <- results(dds, parallel = TRUE, contrast = c("<TO_EDIT>", "<TO_EDIT>", "<TO_EDIT>"), cooksCutoff=FALSE)
```
where:
- `parallel = TRUE`: Parallel execution with `BiocParallel`
- `contrast`: A character vector with three elements for log2 fold-change calculation - full model variable name, numerator and denominator
- `cooksCutoff=FALSE`: The `results` function automatically flags genes which contain a Cook’s distance above a cutoff for samples which have 3 or more replicates. Setting this option to `FALSE` turns this behaviour off.

**Important information: Don't forget to edit the `contrast` character vector.**

#### TASK:
Choose a contrast you would like to use and type in `contrast` parameter in `results` function

We can print the summary of results with `summary` function:
```R
summary(res)
```

### 5.6. Saving results to a file
We can write the table with differentially expressed genes to a file.

First, we need to convert our results to a data frame:
```R
res_file <- as.data.frame(res)
```
Next, we can sort genes by log2 fold-change value:
```R
res_file <- res_file[order(res_file$log2FoldChange, decreasing = TRUE),]
```
We want to add the gene name as a first column of the data frame:
```R
res_file <- as.data.frame(cbind(gene_name = rownames(res_file), res_file))
```
Finally, we can save a data frame to a file:
```R
write.table(res_file, file = "<TO_EDIT>", append = F, quote = F, sep = "\t", row.names = F, col.names = T)
```
**Important information: Don't forget to edit the file name of a file with our results**.

#### TASK:
Open the resulting file either in the text editor (e.g. Notepad) or Excel. What do all the columns mean?
