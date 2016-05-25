Leukerbad Spring School Bioinformatics and Population Genomics

---------------------------------------

# Practical: RNA-seq analysis for population genomics

Julien Roux, version 1, May 2016

## Schedule
- [ ] Wednesday 1 June, 16:45 to 17:45: from raw sequencing data to transcript expression levels.
- [x] **Thursday 2 June, 13:45 to 15:30: gene-level clustering and differential expression analysis.**

## Introduction
Today you will pursue the analysis of Bou Sleiman et al. *Drosophila melanogaster* data. Unless you managed to map all samples yesterday (:clap:), use the pre-processed `Kallisto` results located in the `~/data/rnaseq/kallisto/SRR*` folders. A first step of the practical will be to import the data and sum the transcript-level TPM estimates to gene-level TPM expression estimates. This allows to obtain expression levels that are not affected by differential transcript usage or differential splicing. You will then be able to perform some clustering analyses and study differential expression between strains and experimental conditions. 

In the interest of time, I have written fully most of the R commands, that you can copy-paste. I encourage you to read them fully and try to understand what was done. Feel free to try and modify some parameters, or ask the assistants if sonething is not clear.

## Read and format the metadata

Launch R in the console or using Rstudio.
```R
## path to the data
dataFolder <- "~/data/rnaseq/"
```

You will first read and format the experiment metadata using the GEOquery package (<http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html>). The metadata are located in the `GSE59411_series_matrix.txt` file that we pre-downloaded on the GEO page of the experiment. Metadata can also be read from the `GSE59411_family.soft` file, or by passing the experiment name to the GEOquery package (beware, the format of the resulting R object will slightly differ depending on the approach used). The following steps are a bit boring, but try to understand what was done:

```R
## Should not be needed:
## source("http://bioconductor.org/biocLite.R")
## biocLite("GEOquery")
library(GEOquery)
## Read the metadata file
gse <- getGEO(filename=file.path(dataFolder, "GSE59411_series_matrix.txt"))
## Extract the metadata, and only retain some of the columns
samples <- pData(phenoData(gse))[,c(1,2, 8:12, 45)]
## Extract the SRA ID of samples by parsing the FTP URL
samples$library_accession <- unlist(lapply(strsplit(as.character(samples$supplementary_file_2), split="/"), tail, n=1))
## Read the ENA metadata, which includes the run ID (SRR...)
samplesENA <- read.table(file.path(dataFolder, "SRP044339.txt"), h=T, sep="\t")
## Merge both data frames
samples <- merge(samples, samplesENA, by.x="library_accession", by.y="experiment_accession")
## Filter and rearrange the columns
samples <- samples[,c(2, 4:8, 3, 1, 14)]
## Rename the columns
names(samples)[4] <- "treatment"
names(samples)[5] <- "dgrp_line"
names(samples)[6] <- "resistance"
## Simplify the levels of experimental factors
samples$treatment <- gsub("treatment: ", "", samples$treatment)
samples$dgrp_line <- gsub("dgrp line: ", "", samples$dgrp_line)
samples$resistance <- gsub("resistance: ", "", samples$resistance)
## Visualize the formatted metadata table
head(samples)
dim(samples)
```
![Question](round-help-button.png)
According to the metadata table, how many experimental factors of interest do you identify in this experiment?

For simplicity, we will ignore the strain effect (`dgrp_line`), and consider the four susceptible and the four resistant lines as simple biological replicates.

## Import the transcripts expression levels and sum them at the gene level

You will use the Mike Love's `tximport` package to directly parse `Kallisto` results files and sum expression at the gene-level (<https://github.com/mikelove/tximport>). As `tximport` was not yet available on R/Bioconductor 3.2 installed on the virtual machine (it is now available in release 3.3), you will install it from the source code, that is already in the data folder:
```R
install.packages(file.path(dataFolder, "./tximport-master"), repos = NULL, type="source")
library(tximport)

## Path to Kallisto result files
files <- file.path(dataFolder, "kallisto", samples$run_accession, "abundance.tsv")
names(files) <- samples$title
## Check that all files are present in the data folder
all(file.exists(files))
```

To sum up transcript expression levels, `tximport` needs a data.frame with a transcript ID column and a gene ID column. You will retrieve this information from the Ensembl database via the Biomart webservice. Biomart queries require the specification of 1) a dataset to use, 2) a list of filters to select specific Ensembl genes (not needed in our case), 3) a list of attributes to output. It is possible to build the Biomart query from the web interface ([the needed query would be obtain by specifying the following dataset and attributes](<http://www.ensembl.org/biomart/martview/368fdd3310212bc95ef4d904847c1408?VIRTUALSCHEMANAME=default&ATTRIBUTES=dmelanogaster_gene_ensembl.default.feature_page.ensembl_transcript_id|dmelanogaster_gene_ensembl.default.feature_page.ensembl_gene_id&FILTERS=&VISIBLEPANEL=resultspanel)), but there is a Bioconductor package called `biomaRt` that allows to obtain results directly into R:
```R
library(biomaRt)
## Choose D. melanogaster dataset in Ensembl release 84
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host="mar2016.archive.ensembl.org")
## retrieve the mapping between transcript and gene IDs
transcript2gene <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"), mart=mart) 
## Sort by gene ID
transcript2gene <- transcript2gene[order(transcript2gene$ensembl_gene_id), ]
head(transcript2gene)
```
![Tip](elemental-tip.png) 
Tip: specifying the `host` argument in the `biomaRt` query allows to choose which Ensembl release you wish to work with (in this case, `mar2016.archive.ensembl.org` redirects to release 84). This is very helpful to make your code reproducible.

You will now import `Kallisto` results:
```R
txi <- tximport(files, type = "kallisto", tx2gene = transcript2gene, countsFromAbundance="scaledTPM")
## Visualize the created txi object:
lapply(txi, head)
```
![Question](round-help-button.png)
What are abundances? What are "counts" in this particular case? Why are length estimates not integers?

I suggest to have a look at this interesting paper: 

Soneson C, Love M and Robinson M. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*. 2015;4(1521) (<http://f1000research.com/articles/4-1521/v2>). A PDF of the paper is located in the `~/data/papers/` folder.

In this paper, Soneson and colleagues show that gene-level abundance estimates offer advantages over transcript-level analyses. They are more accurate, more stable, more interpretable, and allow better control on false positives in the presence of differential isoform usage. For this purpose, they introduced the "scaledTPM" values, which are obtained by summing the estimated transcript TPMs from Kallisto within genes, and multiplying with the total library size in millions. ScaledTPM values are artificial values, transforming underlying abundance measures to the scale of read counts. This allows to incorporate the information provided by the sequencing depth, and work with RNA-seq differential expression tools that were developed to use read counts, while being independent on the length of the pool of transcripts expressed for each gene in a given condition.

## Data normalization

We will import the gene expression data into a `edgeR` DGE object and obtain normalization factors for each sample. The default normalization method in `edgeR` is called TMM. It was previously shown to perform quite well, while not too stringent. 
```R
library(edgeR)
library(limma)
y <- DGEList(txi$counts)
y <- calcNormFactors(y)
y$samples
```
![Question](round-help-button.png)
What is the meaning of each column in this data frame?

You will now look at the data distribution in the DGE object. The `cpm` function returns counts per million reads. Because we used scaledTPMs instead of real counts, this corresponds  for each sample to the TPM abundances divided by its TMM normalization factor:
```R
## Load a nice color palette of 9 + 8 colors to be used for plots
library(RColorBrewer)
myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
## Plot CPM distribution
plotDensities(cpm(y), col=myPalette[1:16])
plotDensities(cpm(y, log=T), col=myPalette[1:16])
?cpm ## note that a small value is added to the counts if log=T
```
![Question](round-help-button.png)
What do you observe? How are the data distributed? Why are there two modes? What would be a good cutoff to discriminate the two modes?

It is a good practice to filter our the genes that are not expressed, or very lowly expressed. This alleviates the multiple testing burden, and anyway there is very little power to detect differential expression for these genes. It is of course better to filter these genes after data normalization, before differential expression testing. The choice of a criterion to filter genes is arbitrary, but it has to make sense given the distribution of the data. *It is a bad practice to filter genes based on their variance!*
```R
cutoff <- ... ## put here the cutoff expression thta was chosen in the previous question
summary(cpm(y, log=T) > cutoff)
hist(apply(cpm(y, log=T), 1, function(x){ return(sum(x > cutoff)) }))
## We can for example require that each gene has expression above cutoff in at least half of the samples 
selectedGenes <- names(which(apply(cpm(y, log=T), 1, function(x){ return(sum(x > cutoff)) }) >= 8))
length(selectedGenes)
```

We will now rebuild a new DGE object using only selected genes, and renormalize it:
```R
y <- DGEList(txi$counts[selectedGenes, ]) 
y <- calcNormFactors(y)
plotDensities(cpm(y, log=T), col=myPalette[1:16], legend="topright")
```
![Question](round-help-button.png)
What has changed now? Do you like the distribution of expression levels and their normalization?

Let's look at the distribution of samples depending on the 2 experimental factors:
```R
plotDensities(cpm(y, log=T), group=samples$resistance, col=myPalette[1:2])
plotDensities(cpm(y, log=T), group=samples$treatment, col=myPalette[3:4])
```
![Question](round-help-button.png)
Are there systematic differences in the distribution of expression levels that are linked to one or the other experimental factor? Can you guess which experimental factor is likely to have the strongest effect on expression levels in this experiment?

## Data clustering

### PCA
You will perform a principal component analysis (PCA) on the normalized data using the `prcomp` R command. Input data of the PCA need to be ~normally distributed, so you will use the logged-CPM values previously plotted.

```R
pca <- prcomp(t(cpm(y, log=T)), scale = T)
plot(pca)
summary(pca)
loadings <- pca$rotation
scores <- pca$x
```
![Question](round-help-button.png)
What does the `scale=T` argument means? How much variance is explained by the first 2 principal components? The first 4 principal components? What the `scores` and `loadings`?

Plot the samples projected onto the first two principal components.
```R
plot(scores[,1], scores[,2])
## OK, here is a bit nicer plot ;)
plot(scores[,1], scores[,2], 
     xlab=paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained"), ## indicate the % variance explained by PC1
     ylab=paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained"), ## indicate the % variance explained by PC2
     pch=as.numeric(as.factor(samples$treatment))+15, ## points shape according to treatment
     col=myPalette[as.numeric(as.factor(samples$resistance))], ## points color accoriding to susceptibility
     xlim=c(min(scores[,1]), max(scores[,1])) ,
     ylim=c(min(scores[,n+1]), max(scores[,n+1])+(max(scores[,n+1])-min(scores[,n+1]))/4) ## Let a bit of on top of the plot for legend
)
## Plot legends
legend("topleft", legend=levels(as.factor(samples$treatment)), pch=16:17)
legend("topright", legend=levels(as.factor(samples$resistance)), pch=c(16,16), col=myPalette[1:2])
```
![Question](round-help-button.png)
What experimental factor does correlate with PC1? What is special about PC2? It could be useful to add sample names with the following command (the + 5 allows to put the labels next to the points and not on the points)
```R
text(scores[,1], scores[,2] + 5, samples$title)
```
![Tip](elemental-tip.png)
When labels get too messy, it can be nice to only label the interesting points. The command `identify(scores[,1], scores[,2], samples$title)` allows you to click on the points to reveal their labels. You can also influence the way the labels are placed by clicking slightly above/below/left/right of a point. Press escape to exit the clicking mode.

Plot the samples projected onto PC2 and PC3, then PC3 and PC4. Do you observe segregation of the points by any experimental factor? Are these observations consistent with your previous prediction of the experimental factor with the strongest effect on expression levels?

### Heatmap (skip this if timing is tight)
Another way to visualize the data is to plot a heatmap of expression levels, along with a dendrogram obtained by hierarchical clustering of the samples. The `heatmap` R function allows you to perform this, but the `heatmap.2` from the `gplots` package offers more possibilities. There a lot more packages and heatmap functions to explore if you have particular needs. You will first plot a heatmap using a random selection of 100 genes:
```R
library(gplots)
## create a gradient of 100 colors going from light blue to dark blue
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100)
## select randomly 100 genes and extract their expression
selectedExpression <- cpm(y, log=T)[sample(1:length(cpm(y, log=T)[,1]), 100),]
## Plot the heatmap
heatmap.2(selectedExpression, scale="none", col = colors, margins = c(14, 6), trace='none', denscol="white")
```
![Question](round-help-button.png)
What are the rows and the left tree representing? What are the columns and the top tree representing? What does the color intensity mean? Try to replot the heatmap with a new random selection of 100 genes. Are the clustering stable? 

It is difficult to include the expressionation of all genes to create a readable heatmap. An alternative is to calculate the matrix of pairwise correlation coefficients across all samples and plot a heatmap of this matrix. In addition, the `ColSideColors` and `RowSideColors` arguments allow to better visualize the experimental factors of each sample:
```R
allCors <- cor(cpm(y, log=T), method="spearman", use="pairwise.complete.obs")
heatmap.2( allCors, scale="none", col = colors, margins = c(16, 12), trace='none', denscol="white", RowSideColors=myPalette[1:2][as.integer(as.factor(samples$resistance))], ColSideColors=myPalette[3:4][as.integer(as.factor(samples$treatment))])
```
![Question](round-help-button.png)
What are the rows and the columns representing now? What does the color intensity mean? Is the clustering pattern consistent with the PCA? Do you see a manifestatio of the DGRP lines effect on the heatmap?

## Differential expression analysis

### Some background
There are several softwares to test for differential expression. `DESeq2` and `edgeR` are very good, but here you will use the `limma-voom` method, introduced in this paper:

Law C, Chen Y, Shi W, Smyth G. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*. 2014;15(2):R29 (<http://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29>). A PDF of the paper is located in the `~/data/papers/` folder.
 
While `DESeq2` and `edgeR` work directly on count distributions, using a negative binomial modeling framework, `limma-voom` works on the log-transformed CPMs, which are normally distributed data. This is very useful because the theory behind of normal distributions is more tractable. Notably, a large range of statistical methods that were developped for microarray analysis can be used on such normally distributed data. For example, the widely used `limma` package uses an empirical Bayes approach, which borrows information across genes to adjust variance estimates when sample sizes are small. `Limma` is also a very comprehensive package, able to deal with experiments with complex designs.

A problem of the log-transformation is that it does not yield stable variances (i.e., there is heteroscedasticity). The variability of genes with lower expression is lower than those with high expression, This is easily observable when plotting one sample versus a biological replicate: 
```R
## The color argument was set up to plot transparent black points 
plot(cpm(y, log=T)[,5], cpm(y, log=T)[,6], pch=16, col=rgb(0,0,0,0.2))
```
To deal with this problem, `limma-voom` models the mean-variance relationship observed in the data. This trend is incorporated into a precision weight for each individual normalized observation, which are proportional to the expression levels (see <https://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares>). These weights can be used in `limma` during its linear modeling step. `Limma-voom` was shown to perform equally well compared to tools based on the negative-binomial distribution modeling, or even better when the sequencing depths are different across samples. 

An additional interesting possibility offered by the `limma` package is the possibility to adjust the `limma-voom` precision weights to deal with variations in sample quality, frequently encountered in RNA-seq experiments. Removal of high variation samples reduces noise but leads to a loss of power. A compromise is to use all samples, but to down-weight the observations from more variable samples. This is implemented in the `voomWithQualityWeights` function, described in the following paper:

Liu R, et al. Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. *Nucleic Acids Res*. 2015;43(15):e97 (<http://nar.oxfordjournals.org/content/43/15/e97.long>). A PDF of the paper is located in the `~/data/papers/` folder.

![Question](round-help-button.png)
Following your interpretation on the PCA results, do you think this procedure will be beneficial?

### Implementation
```R
## Limma-voom requires the specification of a design matrix. It is simpler to create a single factor made of the combination of the 2 factors of interest: susceptibility and treatment
condition <- factor(paste(samples$treatment, samples$resistance, sep="."))
## Build the design matrix:
design <- model.matrix(~ 0 + condition) 
## The "~ 0 + ..." syntax is optional, but it will later allow an easier specification of the contrasts for differential expression
## Simplify design matrix column names:
colnames(design) <- gsub("condition", "", colnames(design))
## Apply the limma-voom method:
v <- voomWithQualityWeights(y, design, plot=T)
```
![Question](round-help-button.png)
What are the generated plots indicating? Which sample(s) are the most down-weighted?

### Linear modeling and extraction of interesting contrasts

Please refer to the very nice user guide for `limma` for details on the analysis, or to adapt this pipeline to analyze your own data (<http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>)
```R
fit <- lmFit(v)
```

The next step is to specify the contrasts of interest in a contrast matrix. It is build by constructing linear combinations of the design matrix column names. I have specified for you 3 contrast:
* the treatment effect in resistant lines (CvsUinR)
* the treatment effect in susceptible lines
* the treatment effect in general
Following the same logic, please complete the command to add the contrasts for:
* the resistance effect in challenged lines
* the resistance effect in unchallenged lines
* the resistance effect in general
```R
cont.matrix <- makeContrasts(
                             treatmentInR  = Challenged.Resistant - Unchallenged.Resistant,
                             treatmentInS  = Challenged.Susceptible - Unchallenged.Susceptible,
                             treatment     = (Challenged.Resistant + Challenged.Susceptible) - (Unchallenged.Resistant + Unchallenged.Susceptible),
                             resistanceInC = [...],
                             resistanceInU = [...],
                             resistance    = [...],
                             levels=design)
cont.matrix
```
![Question](round-help-button.png)
Which contrasts do you think are the most biologically interesting?

You can now extract the lists of genes differentially expressed for each contrast, at a 10% FDR.
```R
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2, p.value=0.1) 
summary(results)
```
![Question](round-help-button.png)
What do you observe? Is it consistent with the PCA? What is surprising?

You can see how many genes are differentially expressed in common between 2 contrast using Venn diagram:
```R
vennDiagram(results[,c(1,2)])
vennDiagram(results[,c(2,5)])
vennDiagram(results[,c(2,6)])
```
![Question](round-help-button.png)
What can you conclude about the overlap between lists of DE genes. Do these results make sense?

### Extraction of differentially expressed genes
To extract the lists of differentially expressed genes, the `coef` argument is needed. It corresponds to the column number of the contrast matrix. 
```R
## Treatment:
treatmentGenes <- topTable(fit2, coef=3, p.value=0.1, number=Inf, sort.by="P")
## visualize the top 100 genes
selectedExpression <- cpm(y, log=T)[rownames(treatmentGenes)[1:100],]
heatmap.2(selectedExpression, scale="none", col = colors, margins = c(14, 6), trace='none', denscol="white", ColSideColors=myPalette[3:4][as.integer(as.factor(samples$treatment))])
## Resistance:
resistanceGenes <- topTable(fit2, coef=6, p.value=0.1, number=Inf, sort.by="P")
## visualize all DE genes
selectedExpression <- cpm(y, log=T)[rownames(resistanceGenes),]
heatmap.2(selectedExpression, scale="none", col = colors, margins = c(14, 6), trace='none', denscol="white", ColSideColors=myPalette[1:2][as.integer(as.factor(samples$resistance))])
```
![Question](round-help-button.png)
How does the clustering looks like when only resistanceGenes are taken into account? What does it tell you?

## Bonus: characterization of differentially expressed genes
### Manually (long, subjective, ...)
Because there are relatively few genes differentially expressed between resistant and susceptible lines, it is possible to look them up in reference databases (Ensembl, Uniprot, Flybase). 

### GO enrichment test
However for long lists of genes, it is difficult to get an objective view of what are these genes, and to characterize their function. Gene Ontology (GO) enrichment analyses are useful wiht that respect. For each GO category (a group of genes sharing the same function, involved in the same process, or located in the same cellular compartment) it is possible to test whether the proportion of DE genes is higher then expected. We can do this with the `topGO` package.

![Warning](warning.png)
A word of caution for the use of such tools: 
* They can be very sensitive to the set of genes used as background/universe. In this practical, you pre-filtered genes before DE analysis, so the list of genes that were actually tested should be used as universe. 
* Another problem when applying such tests to a list of DE genes is a potential length bias. There will be more RNA-seq reads for long genes, so more power to detect DE for longer genes. This can confound the ontology enrichment tests which will report functional biases associated to longer genes.
* Since the GO categories are not independent, it is debated whether a mutiple testing correction (such as FDR) is appropriate. There is a nice paragraph on that topic in the `topGO` package documentation. 

```R
library(topGO)

## topGO needs a vector with 0 or 1 values depending if a gene is DE or not
geneList <- (rep(0, times=length(rownames(results))))
names(geneList) <- rownames(results)
DEGenes <- row.names(treatmentGenes)
geneList[DEGenes] <- 1
geneList = as.factor(geneList)
summary(geneList)

## We then need a mapping of genes to the GO categories. This can be retrieved from Ensembl using biomaRt, or using the Drosophila melanogaster annotation package in Bioconductor
## If the annotation package is not installed: 
## source("http://bioconductor.org/biocLite.R")
## biocLite("org.Dm.eg.db")
library(org.Dm.eg.db)
## Build the topGO object for biological process ontology
BPdata <- new("topGOdata",
              ontology="BP",
              allGenes = geneList,
              nodeSize = 5,
              annot = annFUN.org,
              mapping = "org.Dm.eg.db",
              ID = "Ensembl")
resultBP <- runTest(BPdata, algorithm = "classic", statistic = "fisher")
myTable <- GenTable(BPdata, classic = resultBP, topNodes=length(BPdata@graph@nodes), numChar=100)
```
If you have time you can run the GO enrichment test on the molecular function (`ontology="MF"`) or the cellular component ontologies (`ontology="CC"`).

![Question](round-help-button.png)
What are the top categories enriched for genes DE with treatment? Does it make sense?

![Tip](elemental-tip.png)
`TopGO` includes the possibility to use several decorrelation algorithms, giving less redundant, and more precise categories in the results. Repeat the analysis with the weight algorithm (`algorithm = "weight"`), and observe the difference in results. 

### topAnat enrichment test
It is also possible to perform a similar ontology enrichment, but on the fly anatomical ontology instead of Gene Ontology. Genes are mapped to a tissue if some expression was detected in this tissue. With the Bgee database team (<http://bgee.org>), I have developped a tool called `BgeeDB` allowing to do this, based on the `topGO` package algorithm. It is available in the latest Bioconductor release, or on GitHub <https://github.com/BgeeDB/BgeeDB_R>. I encourage you to try in addition to the classical GO enrichment tests, it gives very interesting results! A graphical interface is also available at <http://bgee.org/?page=top_anat#/>.

---------------------------------------
<sub>Icons taken from http://www.flaticon.com/</sub>

<!--
## TO DO: some questions like: what is on the x / y axis

## TO DO: thanks amina
## TO DO: replace cpm(y, log=T) by a variable

## TO DO: how to implement code folding/hiding?
          easiest is probably to have 2 versions, one with code, one without
          or change file names to generic file names

* TO DO: prepare short presentation of: 
  * kallisto. Fast + accurate: game changer
  * DTU/DE/DTE. DE confounded by DTU
  * limma-voom on TPM, etc

![Question](round-help-button.png)
![Tip](elemental-tip.png)
![To do](wrench-and-hammer.png)
![Warning](warning.png)

http://www.emoji-cheat-sheet.com/
-->