### For quality control: Make PCAs with coloring by replicate
library(Biobase)
library(DESeq2)

load("Data/eset_vsd.RData")

legend <- unique(pData(eset)[,c("replicate"), drop = FALSE])
legend$col.replicate <- c("red", "darkgreen", "blue", "darkorange")
rgb.col <- col2rgb(legend$col.replicate)
legend$bg.cols <- apply(rgb.col,
                        2,
                        function(x){
                            rgb(x[1], x[2], x[3], alpha = 200, max = 255)
                        })


rownames(legend) <- paste("FACS day", legend$replicate)

info <- merge(legend, pData(eset),
              by = "replicate",
              sort = FALSE)

info.ind <- match(pData(eset)$original.id, info$original.id)
info.ord <- info[info.ind,]
all( pData(eset)$original.id == info.ord$original.id)
## TRUE

pData(eset) <- info.ord

pdf("Results/pca_replicate.pdf", width = 4, height = 4, bg = "white")
par(xpd=TRUE,
    mfrow = c(2,2))
pca <- prcomp(t(exprs(eset)), scale = TRUE)
s <- summary(pca)$importance[, 1:10]
## make PCA for the first 4 PCs:
for(i in 1:4){
    par(mar=c(4.1, 4.1, 1.1, 1.1))
    pc.x <- paste0("PC", i)
    pc.y <- paste0("PC", i+1)
    plot(pca$x[,pc.x],
         pca$x[,pc.y],
         pch = 21,
         col = pData(eset)$col.replicate,
         bg = pData(eset)$bg.cols,
         cex=1,
         xlab=paste0(pc.x, " (", round(100*s[2, i], digits = 1), "%)"),
         ylab=paste0(pc.y, " (", round(100*s[2, i+1], digits = 1),  "%)"))
}
dev.off()
