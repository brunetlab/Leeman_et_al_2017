
library(Biobase)
library(DESeq2)
library("MASS")

load("Data/eset_vsd_replicate_corrected.RData")


## Spearman correlation and average

pdf("Results/dendrogram_spearman_average.pdf")
par(mar=c(6.1, 4.1, 1.1, 1.1))
cor <- cor(exprs(eset), method = "spearman")
d <- as.dist(1 - cor)
dd <- as.dendrogram(hclust(d, method = "average"))
plot(dd)
dev.off()

