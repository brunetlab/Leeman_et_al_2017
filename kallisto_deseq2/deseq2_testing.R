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

## Remove samples with less than 1 million counts:
ind.keep <- colSums(exprs(eset)) >= 1e6
samples.keep <- names(ind.keep)[ind.keep]
eset <- eset[, samples.keep]
dim(eset)
## Features  Samples
##    22410       35

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
                              design = ~ label + label:age + replicate,
                              colData = pData(eset))

dds.deseq <- DESeq(dds)
save(dds.deseq,
     file = "Data/dds_deseq_with_replicate.RData")

pdf("Data/dispersion_estimate_dds_deseq_with_replicate.pdf")
plotDispEsts(dds.deseq)
dev.off()
## Looks good!

resultsNames(dds.deseq)

## Endo
res.endo <- results(dds.deseq, list("labelEndo.ageold", "labelEndo.ageyoung"))
sum(res.endo$padj < 0.05, na.rm = TRUE)

## Ast
res.ast <- results(dds.deseq, list("labelAst.ageold", "labelAst.ageyoung"))
sum(res.ast$padj < 0.05, na.rm = TRUE)

## qNSC
res.qnsc <- results(dds.deseq, list("labelqNSC.ageold", "labelqNSC.ageyoung"))
sum(res.qnsc$padj < 0.05, na.rm = TRUE)

## aNSC
res.ansc <- results(dds.deseq, list("labelaNSC.ageold", "labelaNSC.ageyoung"))
sum(res.ansc$padj < 0.05, na.rm = TRUE)

## NPC
res.npc <- results(dds.deseq, list("labelNPC.ageold", "labelNPC.ageyoung"))
sum(res.npc$padj < 0.05, na.rm = TRUE)

age.results.deseq <- list(Endo = res.endo,
                          Ast = res.ast,
                          qNSC = res.qnsc,
                          aNSC = res.ansc,
                          NPC = res.npc)

save(age.results.deseq,
     file = "Data/age_results_deseq_with_replicate.RData")


age.results.deseq <- lapply(age.results.deseq, as.data.frame)

## write tables (ordered by p-value):
for(i in seq(along = age.results.deseq)){
    age.res.i <- age.results.deseq[[i]]
    age.res.i <- age.res.i[order(age.res.i$padj), ]
    write.table(age.res.i,
                file = paste0("Results/y_vs_o_", names(age.results.deseq)[i], ".txt"),
                sep = "\t",
                col.names = NA,
                quote = FALSE)
}



## This gives me the comparison between young qNSC and young aNSC:
aNSC.vs.qNSC <- results(dds.deseq, contrast = c(
                                       0,  ## intercept
                                       1, 0, 0, 0, -1, ## cell types
                                       0, 0, 0, 0, ## replicates
                                       0, 0, 0, 0, 0, ## effect old
                                       1, 0, 0, 0, -1 ## effect young
))
sum(aNSC.vs.qNSC$padj < 0.05, na.rm = TRUE)

save(aNSC.vs.qNSC,
     file = "Data/aNSC_vs_qNSC_deseq_with_replicate_only_young.RData")


## This gives me the comparison between old qNSC and old aNSC:
aNSC.vs.qNSC <- results(dds.deseq, contrast = c(
                                       0,  ## intercept
                                       1, 0, 0, 0, -1, ## cell types
                                       0, 0, 0, 0, ## replicates
                                       1, 0, 0, 0, -1, ## effect old
                                       0, 0, 0, 0, 0 ## effect young
))
sum(aNSC.vs.qNSC$padj < 0.05, na.rm = TRUE)

save(aNSC.vs.qNSC,
     file = "Data/aNSC_vs_qNSC_deseq_with_replicate_only_old.RData")
aNSC.vs.qNSC <- aNSC.vs.qNSC[order(aNSC.vs.qNSC$padj), ]


## This gives me the comparison between young Ast and young qNSC:
Ast.vs.qNSC <- results(dds.deseq, contrast = c(
                                      0,  ## intercept
                                      0, -1, 0, 0, 1, ## cell types
                                      0, 0, 0, 0, ## replicates
                                      0, 0, 0, 0, 0, ## effect old
                                      0, -1, 0, 0, 1 ## effect young
))
sum(Ast.vs.qNSC$padj < 0.05, na.rm = TRUE)

Ast.vs.qNSC <- Ast.vs.qNSC[order(Ast.vs.qNSC$padj), ]
write.table(Ast.vs.qNSC,
            file = paste0("Results/Ast_vs_qNSC_deseq_with_replicate_only_young.txt"),
            sep = "\t",
            col.names = NA,
            quote = FALSE)


## This gives me the comparison between young aNSC and young NPC:
aNSC.vs.NPC <- results(dds.deseq, contrast = c(
                                      0,  ## intercept
                                      -1, 0, 0, 1, 0, ## cell types
                                      0, 0, 0, 0, ## replicates
                                      0, 0, 0, 0, 0, ## effect old
                                      -1, 0, 0, 1, 0 ## effect young
))
sum(aNSC.vs.NPC$padj < 0.05, na.rm = TRUE)

aNSC.vs.NPC <- aNSC.vs.NPC[order(aNSC.vs.NPC$padj), ]
write.table(aNSC.vs.NPC,
            file = paste0("Results/aNSC_vs_NPC_deseq_with_replicate_only_young.txt"),
            sep = "\t",
            col.names = NA,
            quote = FALSE)
