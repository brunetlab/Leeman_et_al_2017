## Make barplot of significantly enriched KEGG pathways in activation signature (dormant versus proliferating)

library(Biobase)

## for colors:
load("Data/eset_vsd.RData")

## Created in pathway_enrichment_aNSC_vs_qNSC_only_old.R:
load("Results/kegg_results_aNSC_vs_qNSC_deseq_only_old_with_replicate.RData")
kegg.results <- kegg.results$pathway.results

## Remove "KEGG" and first letter upper case:
kegg.results$pathway <- rownames(kegg.results)
tab.plot <- kegg.results[, c("pathway", "adj.p.val", "test.statistic", "fold.change")]
tab.plot$pathway <- tolower(gsub("KEGG_", "", tab.plot$pathway))
tab.plot$pathway <- gsub("_", " ", tab.plot$pathway)
tab.plot <- tab.plot[nrow(tab.plot):1,]

tab.plot$pathway <- sapply(tab.plot$pathway,
                           function(x){
                               paste0(toupper(substring(x, 1, 1)), substring(x, 2))
                           })
rownames(tab.plot) <- tab.plot$pathway
tab.plot.new <- tab.plot

## Results from the young activated vs quiescent samples for order:
load("Data/pathway_enrichment_activation_only_young_barplot_info.RData")

## Get only the pathways that we showed for young pooled data:
tab.plot.new <- tab.plot.new[tab.plot$pathway,]
tab.plot <- tab.plot.new

## Make barplot:

pdf("Results/pathway_enrichment_aNSC_vs_qNSC_only_old_barplot.pdf", width = 7, height = 8.5)
par(mar=c(4.1, 16, .1, 3.2), xpd = TRUE)
x <- tab.plot$fold.change
col <- character(length = length(x))
col[x < 0] <- pData(eset)$cols.pca[pData(eset)$label == "qNSC"][1]
col[x > 0] <- pData(eset)$cols.pca[pData(eset)$label == "aNSC"][1]
b <- barplot(x,
             horiz = TRUE,
             names.arg = tab.plot$pathway,
             xlab = "Average log2 fold change (aNSC / qNSC)",
             border = NA,
             las = 1,
             col = col,
             axes = FALSE,
             space = 0.3)
stars <- character(length = nrow(tab.plot))
pvals <- tab.plot$adj.p.val
stars[ pvals > 0.05 ] <- ""
stars[ pvals <= 0.05 ] <- "*"
stars[ pvals <= 0.001 ] <- "**"
stars[ pvals <= 0.0001 ] <- "***"
pos <- integer(length = nrow(tab.plot))
pos[ x < 0 ] <- 4
pos[ x > 0 ] <- 2
text(x = x,
     y = b[,1] - 0.3,
     stars,
     pos = pos)
axis(1,
     at = seq(-4, 4, by = 1))
dev.off()
