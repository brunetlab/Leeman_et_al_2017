## Significant genes get the colors of the cell types.

library(Biobase)

load("Data/age_results_deseq_with_replicate.RData")
age.results <- lapply(age.results.deseq, as.data.frame)

## for colors:
load("Data/eset_vsd.RData")
ct.unique <- unique(pData(eset)[order(pData(eset)$label),c("label", "cols.pca", "cols.legend")])
rownames(ct.unique) <- ct.unique$label
ct.unique <- ct.unique[c("Endo", "Ast", "qNSC", "aNSC", "NPC"),]

## Remove rows with NA p-values:
age.results <- lapply(age.results,
                      function(x) {
                          x[!is.na(x$padj), ]
                      })
n <- sapply(age.results, nrow)
names(n) <- names(age.results)

## Order by pvalue:
age.results <- lapply(age.results,
                      function(x) {
                          x[order(x$padj),]
                      })

p <- 0.05
cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
    cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
    ind.sig.i <- age.results[[i]]$padj < p
    cols[[i]][ind.sig.i] <- ct.unique[names(age.results)[i], "cols.pca"]
    xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf("Results/stripplot_deseq_with_replicate_cell_type_colors.pdf", width = 5, height = 3)
par(mar = c(3.1, 4.1, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 5.5),
     ylim = c(-4, 6),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change (old / young)"
     )
abline(h = 0)
abline(h = seq(-4, 6, by = 2)[-3],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
    set.seed(1234)
    points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
           y = rev(age.results[[i]]$log2FoldChange),
           pch = 21,
           col = rev(cols[[i]]),
           bg = rev(cols[[i]]))
}
axis(1,
     at = 1:5,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab)
axis(2,
     las = 1,
     at = seq(-6, 6, 2))
dev.off()
