### For all plots I use vst normalized values.

library(Biobase)
library(DESeq2)

load("Data/eset_gene_counts_all_samples.RData")
eset <- eset.gene.counts
rm(eset.gene.counts)

## Remove the samples that don't pass quality checks (See Materials and Methods):
out.1 <- which(colnames(eset) == "qNSC_young_C")
out.2 <- which(colnames(eset) == "qNSC_young_E")
out.3 <- which(colnames(eset) == "qNSC_old_C")
out.4 <- which(colnames(eset) == "aNSC_young_C")
eset <- eset[, -c(out.1, out.2, out.3, out.4)]

## Remove samples with less than 1 million counts:
ind.keep <- colSums(exprs(eset)) >= 1e6
samples.keep <- names(ind.keep)[ind.keep]
eset <- eset[, samples.keep]
dim(eset)

## Filter out genes with low coverage:
cov.per.sample <- colSums(exprs(eset))
norm.fact <- rep(cov.per.sample, each = nrow(eset))
cpm <- exprs(eset) / norm.fact * 10^6
ind.keep <- rowSums(cpm > 1) >= 2
## I also tried it with cpm > 1, but this doesn't really help.
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]
dim(eset)

## DESeq only takes integers. So I will have to round the estimated counts:

counts.rounded <- round(exprs(eset))

dds <- DESeqDataSetFromMatrix(countData = counts.rounded,
                              design = ~ label + label:age,
                              colData = pData(eset))
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- varianceStabilizingTransformation(dds)

### vst transformation:
exprs(eset) <- assay(vsd)
save(eset, file = "Data/eset_vsd.RData")
