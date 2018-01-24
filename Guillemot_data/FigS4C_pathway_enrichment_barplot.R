## Show the same pathways as we show for our data.

setwd("~/Detja/Published_data/Guillemot/")
library(Biobase)
## for colors:
load("../kallisto_deseq2/Data/eset_vsd_replicate_corrected.RData")

## Created in pathway_enrichment.R:
load("Data/kegg_results.RData")
kegg.results <- kegg.results$pathway.results

## Remove "KEGG" and first letter upper case:
kegg.results$pathway <- rownames(kegg.results)
tab.plot.public <- kegg.results[, c("pathway", "adj.p.val", "test.statistic", "fold.change")]
tab.plot.public$pathway <- tolower(gsub("KEGG_", "", tab.plot.public$pathway))
tab.plot.public$pathway <- gsub("_", " ", tab.plot.public$pathway)

tab.plot.public$pathway <- sapply(tab.plot.public$pathway,
                                  function(x){
                                      paste0(toupper(substring(x, 1, 1)), substring(x, 2))
                                  })
rownames(tab.plot.public) <- tab.plot.public$pathway

## Results from our data for order of pathways:
load("../kallisto_deseq2/Data/pathway_enrichment_activation_only_young_barplot_info.RData")

## Get only the pathways that we showed for our data:
tab.plot.public <- tab.plot.public[tab.plot$pathway,]
tab.plot <- tab.plot.public

## Make barplot:

pdf("Results/pathway_enrichment_barplot.pdf", width = 7, height = 8.5)
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
     at = seq(-1, 2, by = 1))
dev.off()
