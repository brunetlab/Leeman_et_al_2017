library(DESeq2)
source("../kallisto_deseq2/Data/KEGG_pathway_enrichment.R")

load("Data/results_qNSC_vs_aNSC.RData")

test.results <- as.data.frame(res)
rownames(test.results) <- toupper(rownames(test.results))

#### Based on test statistic:

set.seed(1234)
kegg.results <- KEGG.pathway.enrichment(gene.sets.file = "../kallisto_deseq2/Data/c2.cp.kegg.v5.0.symbols.gmt",
                                        gene.sets.list = NULL,
                                        test.results = test.results,
                                        statistic = "stat",
                                        num.samples = 10000,
                                        abs = FALSE)

ind.sig <- kegg.results$pathway.results$adj.p.val < 0.05

kegg.results$pathway.results[ind.sig, "test.statistic", drop = FALSE]

save(kegg.results,
     file = "Data/kegg_results.RData")


