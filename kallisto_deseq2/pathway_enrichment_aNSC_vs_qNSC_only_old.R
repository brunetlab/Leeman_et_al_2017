## Test enrichment of KEGG pathways in activation signature (dormant versus proliferating)

source("Scripts/KEGG_pathway_enrichment.R")

## Created in deseq2_testing.R:
load("Data/aNSC_vs_qNSC_deseq_with_replicate_only_old.RData")
test.results <- as.data.frame(aNSC.vs.qNSC)
rownames(test.results) <- toupper(rownames(test.results))

set.seed(1234)
kegg.results <- KEGG.pathway.enrichment(gene.sets.file = "~/Detja/STAR/Data/c2.cp.kegg.v5.0.symbols.gmt",
                                        gene.sets.list = NULL,
                                        test.results = test.results,
                                        statistic = "stat",
                                        num.samples = 10000,
                                        abs = FALSE)
sum(kegg.results$pathway.results$adj.p.val < 0.05)

save(kegg.results,
     file = "Results/kegg_results_aNSC_vs_qNSC_deseq_only_old_with_replicate.RData")
