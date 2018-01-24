## Test enrichment of KEGG pathways in activation signature (dormant versus proliferating)

source("~/My_R_functions/KEGG_pathway_enrichment.R")

## Created in deseq2_testing_activation_only_old.R:
load("Data/activation_results_deseq_only_old_with_replicate.RData")
test.results <- as.data.frame(activation.results.deseq)
rownames(test.results) <- toupper(rownames(test.results))

set.seed(1234)
kegg.results <- KEGG.pathway.enrichment(gene.sets.file = "Data/c2.cp.kegg.v5.0.symbols.gmt",
                                        gene.sets.list = NULL,
                                        test.results = test.results,
                                        statistic = "stat",
                                        num.samples = 10000,
                                        abs = FALSE)
sum(kegg.results$pathway.results$adj.p.val < 0.05)

save(kegg.results,
     file = "Results/kegg_results_activation_deseq_only_old_with_replicate.RData")
