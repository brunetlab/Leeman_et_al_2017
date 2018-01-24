library(DESeq2)

dds.deseq <- DESeq(dds)

res <- results(dds.deseq, list("labelaNSC", "labelqNSC"))

save(res,
     file = "Results/results_qNSC_vs_aNSC.RData")
