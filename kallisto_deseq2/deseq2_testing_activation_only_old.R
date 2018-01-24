library(Biobase)
library(DESeq2)

load("Data/eset_gene_counts_all_samples.RData")
eset <- eset.gene.counts
rm(eset.gene.counts)

## Remove the samples that have weird density curves:
out.1 <- which(colnames(eset) == "qNSC_young_C")
out.2 <- which(colnames(eset) == "qNSC_young_E")
out.3 <- which(colnames(eset) == "qNSC_old_C")
out.4 <- which(colnames(eset) == "aNSC_young_C")
eset <- eset[, -c(out.1, out.2, out.3, out.4)]

## Remove Endothelials:
eset <- eset[, pData(eset)$label != "Endo"]

## Remove all young samples:
eset <- eset[, pData(eset)$age != "young"]

pData(eset)$activation[ is.element(pData(eset)$label, c("Ast", "qNSC")) ] <- "quiescent"
pData(eset)$activation[ is.element(pData(eset)$label, c("aNSC", "NPC")) ] <- "activated"
pData(eset)$activation <- relevel(as.factor(pData(eset)$activation), "quiescent")

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
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]
dim(eset)


## DESeq only takes integers. So I will have to round the estimated counts:

counts.rounded <- round(exprs(eset))

dds <- DESeqDataSetFromMatrix(countData = counts.rounded,
                              design = ~ activation + replicate,
                              colData = pData(eset))
dds.deseq <- DESeq(dds)

resultsNames(dds.deseq)

res <- results(dds.deseq, list("activationactivated", "activationquiescent"))
sum(res$padj < 0.05, na.rm = TRUE)


activation.results.deseq <- res

save(activation.results.deseq,
     file = "Data/activation_results_deseq_only_old_with_replicate.RData")
