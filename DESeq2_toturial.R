setwd("C:\\Users\\wentingwang\\Desktop\\NGS\\DE\\")

# ---------------------------- #
# Author: DungChi Wu
# Date: Fri Dec  8 15:02:40 2017
# Description: Brief tutorial for DESeq2, modified from following reference.
# Reference: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#can-i-use-deseq2-to-analyze-a-dataset-without-replicates
# ---------------------------- #

# 0. Setup ----------------------------------------------------------------
# First, install the libraries we'll need.
# If you have them installed, skip this step.

# From cran:
install.packages("RColorBrewer")

# From bioconductor:
source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
biocLite("pasilla")


# 1. Preparation ----------------------------------------------------------
# Before we start, we would need to load the necessary libraries into R by 
# the command: library("name")

library("DESeq2")
library("RColorBrewer")
library("pasilla")


# Preparing the input data: 
#
# DESeq2 needs a "counts matrix".
# The rows are genes while the columns are samples.
#
# In the following example, we will use the pasilla dataset.
# This data set is from an experiment on Drosophila melanogaster cell cultures 
# and investigated the effect of RNAi knock-down of the splicing factor pasilla 
# (Brooks et al. 2011). 
# The detailed transcript of the production of the pasilla data is provided in 
# the vignette of the data package pasilla.

# We read in a count matrix, which we will name cts, 
# and the sample information table, which we will name coldata.

pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]


# Note: the columns of the count matrix and the rows of the column data 
# (information about samples) must be in the same order and consistent. 
# So we need to chop off the "fb" of the row names of coldata and re-arrange 
# one or the other so that they are consistent in terms of sample order.

rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))


# 2. Run DESeq2 -----------------------------------------------------------

# With the count matrix, cts, and the sample information, coldata, 
# we can construct a DESeqDataSet:
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds


# The standard differential expression analysis steps 
# are wrapped into a single function, DESeq. 

# Results tables are generated using the function results, 
# which extracts a results table with 
# log2 fold changes, p values and adjusted p values. 

# Details about the comparison are printed to the console, above the results table. 
# The text, condition treated vs untreated, 
# tells you that the estimates are of 
# the logarithmic fold change log2(treated/untreated).

dds$condition <- relevel(dds$condition, ref = "untreated")
dds <- DESeq(dds)
res <- results(dds)
res

# Note: 
# The estimation steps performed by this function are described in 
# the manual page for ?DESeq and in the Methods section of 
# the DESeq2 publication (Love, Huber, and Anders 2014).

# With no additional arguments to results, 
# the log2 fold change and Wald test p value will be for 
# the last variable in the design formula, and if this is a factor, 
# the comparison will be the last level of this variable over the first level. 
# However, the order of the variables of the design do not matter 
# so long as the user specifies the comparison using the name 
# or contrast arguments of results (described later and in ?results).




# Exploring and exporting results from DESeq2 ------------------------------

# Visualization of DESeq2 results:
# We can first use the function plotMA to shows 
# the log2 fold changes to a given variable over the mean of normalized counts 
# for all the samples in the DESeqDataSet. 
# Points will be colored red if the adjusted p value is less than 0.05. 
# Points which fall out of the window are plotted as 
# open triangles pointing either up or down.

plotMA(res, ylim=c(-5, 5) )

# Or we can look at the distribution of adjusted p-value of results

hist(res$padj, breaks = 100, col = "skyblue", border = "slateblue", main = "")


# Filtering and exporting:
# We can first filter the results by an adjusted p value threshold 0.05, and sorted
# it by the adjusted p value in asecnding order: 

resSig <- subset(res[order(res$padj),], padj < 0.05)
resSig

# We can export the results by the base R function:

write.csv(as.data.frame(resSig), file="condition_treated_results.csv")


# Data visualization ------------------------------------------------------

# To further explore the input data, visualization might be useful.
# For example, we can produce a heatmap by R with the library: pheatmap.
# If you don't have it, please installed first.

library("pheatmap")

# But before visualization, we need to do some transformation.
# Maybe the most obvious choice of transformation is the logarithm.
# y = log2(n + 1)  or  y = log2(n + n0), where n is the raw counts.

# We would try normalization by normTransform function in DESeq2,
# this gives log2(n + 1)

ntd <- normTransform(dds)

# Then we select the 30 highestly expressed gene and prepare the annotation.
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,c("condition","type")])

# Plot heatmaps
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# We can also draw a heatmap with raw counts, which is less informative. 
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



# We can also select all the significant differentially expressed gene by which 
# adjusted p value < 0.05.
selectSig <- rownames(subset(res, padj < 0.05 )) 
sig_mat <- assay(ntd)[selectSig, ]
sig_mat <- sig_mat - rowMeans(sig_mat)
pheatmap(sig_mat, scale = "row", cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)


# Principle Component Analysis: PCA plots
plotPCA(ntd, intgroup=c("condition", "type"), ntop = length(assay(ntd)))


# Plot Sample distance:
colors_rev <- colorRampPalette(rev(brewer.pal(9, "Blues")))( n = 299)
sampleDists <- dist(t(assay(ntd)))
sampleDistsMatrix <- as.matrix(sampleDists)
rownames(sampleDistsMatrix) <- paste(ntd$condition, ntd$type, sep="-")
pheatmap(sampleDistsMatrix,                      # distance matrix    
         display_numbers = TRUE,                 # show the distance number  
         color = colors_rev                      # color
)


     