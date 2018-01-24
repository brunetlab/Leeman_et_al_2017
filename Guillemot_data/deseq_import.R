library(Biobase)
library(DESeq2)

## there are no replicates, just 2 lanes per sample -> I add the counts:
tab1 <- read.table("Data/Gene_counts/counts_aNSC_s1.txt", sep = "\t")
tab2 <- read.table("Data/Gene_counts/counts_aNSC_s2.txt", sep = "\t")
all(rownames(tab1) == rownames(tab2))
## TRUE
tab <- tab1
tab$V2 <- tab$V2 + tab2$V2
write.table(tab, "Data/Gene_counts/Merged/counts_aNSC.txt", quote = FALSE,
            col.names = FALSE, row.names = FALSE, sep = "\t")

aNSC <- read.table("Data/Gene_counts/Merged/counts_aNSC.txt",
                   sep = "\t",
                   row.names = 1,
                   col.names = "aNSC")
qNSC <- read.table("Data/Gene_counts/Merged/counts_qNSC.txt",
                   sep = "\t",
                   row.names = 1,
                   col.names = "qNSC")
all(rownames(aNSC) == rownames(qNSC))
## TRUE

counts <- cbind(qNSC, aNSC)

## Remove lines that are not genes:
ind <- grep("_", rownames(counts))
rownames(counts)[ind]
## [1] "__no_feature"           "__ambiguous"            "__too_low_aQual"
## [4] "__not_aligned"          "__alignment_not_unique"
counts.out <- colSums(counts[ind,])
counts.out/colSums(counts)
##      qNSC      aNSC
## 0.2684599 0.2679903
# 27% of the reads are not in exons
counts <- counts[-ind,]


## Filter out genes with low coverage:
cov.per.sample <- colSums(counts)
norm.fact <- rep(cov.per.sample, each = nrow(counts))
cpm <- counts / norm.fact * 10^6
ind.keep <- rowSums(cpm > 1) >= 1
genes.keep <- names(ind.keep)[ind.keep]
counts <- counts[genes.keep, ]
nrow(counts)
## 12206

dds <- DESeqDataSetFromMatrix(countData = counts,
                              design = ~ label,
                              colData = data.frame(label = colnames(counts)))
dds <- estimateSizeFactors(dds)


dds <- estimateDispersions(dds)
save(dds,
     file = "Data/dds.RData")

vsd <- varianceStabilizingTransformation(dds)

eset <- ExpressionSet(assayData = assay(vsd),
                      phenoData = AnnotatedDataFrame(data.frame(label = colnames(counts))))

save(eset,
     file = "Data/eset.RData")

