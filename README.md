# Single cell RNA-Seq course
**Author: Tomasz Wrzesinski**

**Date: 06-07/12/2018**

**Venue: Earlham Institute**

## Software prerequisites


## Data repository
You can find all the input files needed for the analysis at [https://drive.google.com/drive/folders/1K_RfWQVglVsjQ5eEKBYhXwrU5XRynDKi?usp=sharing](https://drive.google.com/drive/folders/1K_RfWQVglVsjQ5eEKBYhXwrU5XRynDKi?usp=sharing)

## Datasets
During our course we will be using three datasets - one focused on expression patterns in developing mouse embryo and two concentrated on different cell populations in human pancreas.

### Deng Q, _et al_., 2014
Link to the publication: [http://science.sciencemag.org/content/343/6167/193](http://science.sciencemag.org/content/343/6167/193)

**Abstract**:
Expression from both alleles is generally observed in analyses of diploid cell populations, but studies addressing allelic expression patterns genome-wide in single cells are lacking. Here, we present global analyses of allelic expression across individual cells of mouse preimplantation embryos of mixed background (CAST/EiJ × C57BL/6J). We discovered abundant (12 to 24%) monoallelic expression of autosomal genes and that expression of the two alleles occurs independently. The monoallelic expression appeared random and dynamic because there was considerable variation among closely related embryonic cells. Similar patterns of monoallelic expression were observed in mature cells. Our allelic expression analysis also demonstrates the de novo inactivation of the paternal X chromosome. We conclude that independent and stochastic allelic transcription generates abundant random monoallelic expression in the mammalian cell.

### Muraro MJ, _et al_., 2016
Link to the publication:
[https://www.sciencedirect.com/science/article/pii/S2405471216302927](https://www.sciencedirect.com/science/article/pii/S2405471216302927)

**Abstract**:
To understand organ function, it is important to have an inventory of its cell types and of their corresponding marker genes. This is a particularly challenging task for human tissues like the pancreas, because reliable markers are limited. Hence, transcriptome-wide studies are typically done on pooled islets of Langerhans, obscuring contributions from rare cell types and of potential subpopulations. To overcome this challenge, we developed an automated platform that uses FACS, robotics, and the CEL-Seq2 protocol to obtain the transcriptomes of thousands of single pancreatic cells from deceased organ donors, allowing in silico purification of all main pancreatic cell types. We identify cell type-specific transcription factors and a subpopulation of REG3A-positive acinar cells. We also show that CD24 and TM4SF4 expression can be used to sort live alpha and beta cells with high purity. This resource will be useful for developing a deeper understanding of pancreatic biology and pathophysiology of diabetes mellitus.

### Segerstolpe A, _et al_., 2016
Link to the publication:
[https://www.sciencedirect.com/science/article/pii/S1550413116304363](https://www.sciencedirect.com/science/article/pii/S1550413116304363)

**Abstract**:
Hormone-secreting cells within pancreatic islets of Langerhans play important roles in metabolic homeostasis and disease. However, their transcriptional characterization is still incomplete. Here, we sequenced the transcriptomes of thousands of human islet cells from healthy and type 2 diabetic donors. We could define specific genetic programs for each individual endocrine and exocrine cell type, even for rare δ, γ, ε, and stellate cells, and revealed subpopulations of α, β, and acinar cells. Intriguingly, δ cells expressed several important receptors, indicating an unrecognized importance of these cells in integrating paracrine and systemic metabolic signals. Genes previously associated with obesity or diabetes were found to correlate with BMI. Finally, comparing healthy and T2D transcriptomes in a cell-type resolved manner uncovered candidates for future functional studies. Altogether, our analyses demonstrate the utility of the generated single-cell gene expression resource.

## Input files
### Count matrix
The count matrix contains gene IDs as rows and sample names as columns.
### Annotation file
The annotation file contains sample names as row and all available meta-data (conditions/batches) as columns.

## DESEq2 analysis
### Loading the data to R
By default, all character strings in the table are being recognised as factors. To change this behaviour we need to type:
```
R
# Changing the default R behaviour that all strings are treated as factors when reading the table
options(stringsAsFactors = FALSE)
```

## Expression QC

**Useful links**




## Differential Expression Testing
# This section is based on:

## More information
