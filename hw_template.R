# ---------------------------- #
# Author: DungChi Wu
# Date: Sat Dec  9 17:56:41 2017
# Description: 
# ---------------------------- #

# Following is the exercise: 
# You would first need to set up the working directory
setwd("")

# Load necessary library
library("DESeq2")
library("RColorBrewer")
library("ggplot2")


# And then loading the expression value data.
# value data <-read.table()

# Also the metadata for each experiment.
# metadata <-read.table()

# You could using head(rawCts) or View(rawCts) (if you use RStudio) to peek the
# data.



# You may notice that the rowname is very difficult to read, since it contain 
# information about transcript from multiple source. We could did a little trick 
# to keep only the ENSEMBL Transcript id which started with ENST and remove all
# the other name.
# rownames() <- sub("\\|.*", "", rownames())                   # Keep ENSEMBL ID with version
# rownames() <- sub("\\.\\d+\\|.*" , "", rownames(), perl = T) # Keep ENSEMBL ID without version



# Check if experiment (count column name) is contained in the metadata
all(rownames() %in% colnames())
all(rownames() == colnames())


# Because this dataset contains multiple sample from different experiment, 
# the cell types also differ, so we would first see how the sample correlated.
# To make the variation similar across orders of magnitude. We would first do 
# log2 + 1 transformed, and the calculate the pearson correlation r. 
# And then draw the plots.

# Setup some convient color variable for future used.
colors_rev <- colorRampPalette( rev(brewer.pal(9, "Blues")))( n = 299)
colors_for <- colorRampPalette( brewer.pal(9, "Blues"))(n = 299)
colors_rgb <- colorRampPalette( c("green", "black", "red"))(n = 299)


library("pheatmap")
#sampleCor <- cor(log2(1+value data), method = "pearson")
jpeg("pheatmap.jpg",height = 1024,width= 1024)
pheatmap(sampleCor,              # correlation matrix 
         col = colors_for,       # color used in heatmap
         annotation = meta,   # show the information
         display_numbers = TRUE  # show the correlation value
)
dev.off()

# You may pick the two group to plot again




# Exercise ----------------------------------------------------------------
# Based on the above information, the two cell lines A549 and K562 are distinguishable from each other.
# Please select these experiment data and the related
# information from rawCts and rawInfo, and do the following:

# 1. Run the DESeq pipeline on these two cell lines, A549 and K562
#    e.g. using K562 as base level





# 2. Report the number of differentially expressed gene using following criteria: 
#    padj < 1e-5 AND log2FoldChange > 4 or < -4 and save it into a file.




