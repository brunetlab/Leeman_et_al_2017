### 3d PCAs for each cell type.

library(Biobase)
library(scatterplot3d)

load("Data/eset_vsd_replicate_corrected.RData")

pic.size <- 4


pdf("Results/final_pca3d_Endo.pdf", width = pic.size, height = pic.size - 0.5, bg = "white")
eset.part <- eset[,pData(eset)$label == "Endo"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]
pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:6]
scatterplot3d(x = pca$x[, c("PC1", "PC2", "PC3")],
              pch = pData(eset.part)$pch,
              color = pData(eset.part)$cols.legend,
              bg = pData(eset.part)$cols.pca,
              type = "h",
              box = FALSE,
              angle = 135,
              cex.symbols = 2,
              mar=c(3, 2, 0, 3)+0.1,
              y.margin.add = 0.2,
              xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
              ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
              zlab=paste("PC3 (", round(100*s[2,3], digits = 1),  "%)", sep = ""))
dev.off()


pdf("Results/final_pca3d_Ast.pdf", width = pic.size, height = pic.size - 0.5, bg = "white")
eset.part <- eset[,pData(eset)$label == "Ast"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]
pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:6]
scatterplot3d(x = pca$x[, c("PC1", "PC2", "PC3")],
              pch = pData(eset.part)$pch,
              color = pData(eset.part)$cols.legend,
              bg = pData(eset.part)$cols.pca,
              type = "h",
              box = FALSE,
              angle = 70,
              cex.symbols = 2,
              mar=c(3, 3, 0, 2)+0.1,
              y.margin.add = 0.4,
              xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
              ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
              zlab=paste("PC3 (", round(100*s[2,3], digits = 1),  "%)", sep = ""))
dev.off()


pdf("Results/final_pca3d_qNSC.pdf", width = pic.size, height = pic.size - 0.5, bg = "white")
eset.part <- eset[,pData(eset)$label == "qNSC"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]
pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:5]
scatterplot3d(x = pca$x[, c("PC2", "PC3", "PC1")],
              pch = pData(eset.part)$pch,
              color = pData(eset.part)$cols.legend,
              bg = pData(eset.part)$cols.pca,
              type = "h",
              box = FALSE,
              angle = 135,
              cex.symbols = 2,
              mar=c(3, 2, 0, 3)+0.1,
              y.margin.add = 0.2,
              xlab=paste("PC2 (", round(100*s[2,2], digits = 1), "%)", sep = ""),
              ylab=paste("PC3 (", round(100*s[2,3], digits = 1),  "%)", sep = ""),
              zlab=paste("PC1 (", round(100*s[2,1], digits = 1),  "%)", sep = ""))
dev.off()

pdf("Results/final_pca3d_aNSC.pdf", width = pic.size, height = pic.size - 0.5, bg = "white")
eset.part <- eset[,pData(eset)$label == "aNSC"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]
pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:6]
scatterplot3d(x = pca$x[, c("PC2", "PC3", "PC1")],
              pch = pData(eset.part)$pch,
              color = pData(eset.part)$cols.legend,
              bg = pData(eset.part)$cols.pca,
              type = "h",
              box = FALSE,
              angle = 130,
              cex.symbols = 2,
              mar=c(3, 2, 0, 3)+0.1,
              y.margin.add = 0.2,
              xlab=paste("PC2 (", round(100*s[2,2], digits = 1), "%)", sep = ""),
              ylab=paste("PC3 (", round(100*s[2,3], digits = 1),  "%)", sep = ""),
              zlab=paste("PC1 (", round(100*s[2,1], digits = 1),  "%)", sep = ""))
dev.off()


pdf("Results/final_pca3d_NPC.pdf", width = pic.size, height = pic.size - 0.5, bg = "white")
eset.part <- eset[,pData(eset)$label == "NPC"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]
pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:6]
scatterplot3d(x = pca$x[, c("PC2", "PC1", "PC3")],
              pch = pData(eset.part)$pch,
              color = pData(eset.part)$cols.legend,
              bg = pData(eset.part)$cols.pca,
              type = "h",
              box = FALSE,
              angle = 50,
              cex.symbols = 2,
              mar=c(3, 3, 0, 2)+0.1,
              y.margin.add = 0.4,
              xlab=paste("PC2 (", round(100*s[2,2], digits = 1), "%)", sep = ""),
              ylab=paste("PC1 (", round(100*s[2,1], digits = 1),  "%)", sep = ""),
              zlab=paste("PC3 (", round(100*s[2,3], digits = 1),  "%)", sep = ""))
dev.off()
