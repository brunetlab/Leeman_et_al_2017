library(Biobase)
load("Data/eset_vsd_replicate_corrected.RData")

legend <- unique(pData(eset)[,c("label", "age", "pch", "cols.pca", "cols.legend")])
rownames(legend) <- paste(legend$label, legend$age, sep = ".")
legend <- legend[c(2,1,6,5,10,9,8,7,4,3),]

lwd.axes = 1.327
cex = 1.7

pdf("Results/final_pca.pdf", width = 5.2, height = 4.2, bg = "white")
par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 6.6))
pca <- prcomp(t(exprs(eset)), scale = TRUE)
s <- summary(pca)$importance[, 1:6]
plot(pca$x,
     pch = pData(eset)$pch,
     col = pData(eset)$cols.legend,
     bg = pData(eset)$cols.pca,
     cex = cex,
     axes = FALSE,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
     xlim = c(-80, 100),
     ylim = c(-80, 120))
box(lwd = lwd.axes)
axis(1,
     lwd = 0,
     lwd.tick = lwd.axes)
axis(2,
     lwd = 0,
     lwd.tick = lwd.axes)
legend(110, 100,
       col = legend$cols.legend,
       pt.bg = legend$cols.pca,
       pch = legend$pch,
       pt.cex = cex,
       legend = paste(legend$label, legend$age),
       bty = "n")
dev.off()
