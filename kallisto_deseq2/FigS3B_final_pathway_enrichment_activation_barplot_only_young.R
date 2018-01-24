## Make barplot of significantly enriched KEGG pathways in activation signature (dormant versus proliferating)

## Created in pathway_enrichment_activation_only_young.R:
load("Results/kegg_results_activation_deseq_only_young_with_replicate.RData")

## Categories of the pathways:
cat <- read.table("Data/kegg_results_activation_only_young_categorization.txt",
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  header = TRUE)

## Show only a subset of the pathways in the main manuscript:

pathways.main <- cat$Pathway[which(cat$In.main)]

kegg.results <- kegg.results$pathway.results[pathways.main,]
kegg.results <- merge(kegg.results,
                      cat[,c("Pathway", "Category")],
                      by.x = "row.names",
                      by.y = "Pathway")

## Sort by category and within each category by test statistic:

categories <- c("Protein synthesis and degradation",
                "Metabolism",
                "Mitochondria",
                "Intercellular communication",
                "Proliferation",
                "Other")

tab.plot <- data.frame(pathway = character(),
                       adj.p.val = numeric(),
                       test.statistic = numeric(),
                       categ = character())
for(i in seq(along = categories)){
    cat.i <- categories[i]
    res.i <- kegg.results[kegg.results$Category == cat.i,]
    res.i.ord <- res.i[order(res.i$test.statistic),]
    tab.plot.i <- res.i.ord[, c("Row.names", "adj.p.val", "test.statistic", "Category")]
    colnames(tab.plot.i) <- colnames(tab.plot)
    tab.plot <- rbind(tab.plot,
                      tab.plot.i)
}

## Remove "KEGG" and first letter upper case:
tab.plot$pathway <- tolower(gsub("KEGG_", "", tab.plot$pathway))
tab.plot$pathway <- gsub("_", " ", tab.plot$pathway)
tab.plot <- tab.plot[nrow(tab.plot):1,]

tab.plot$pathway <- sapply(tab.plot$pathway,
                           function(x){
                               paste0(toupper(substring(x, 1, 1)), substring(x, 2))
                           })

## I will need the order of the barplot later for the published (Guillemot) data:
save(tab.plot,
     file = "Data/pathway_enrichment_activation_barplot_info.RData")


## Make barplot:

pdf("Results/final_pathway_enrichment_activation_barplot.pdf", width = 7, height = 9)
par(mar=c(4.1, 16, .1, 3.2), xpd = TRUE)
x <- tab.plot$test.statistic
col <- character(length = length(x))
col[x < 0] <- rgb(99, 167, 109, max = 255, alpha = 150)
col[x > 0] <- rgb(103, 72, 146, max = 255, alpha = 150)
b <- barplot(x,
             horiz = TRUE,
             names.arg = tab.plot$pathway,
             xlab = "Average test statistic",
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
     at = seq(-5, 15, by = 5))
dev.off()
