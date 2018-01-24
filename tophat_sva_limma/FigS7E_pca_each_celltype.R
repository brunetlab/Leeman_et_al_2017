library(sva)
library(limma)
library(Biobase)
# Created in sv.R:
load("Data/sva_eset_no_age_no_celltype.RData")
eset <- eset.age.celltype

ct.unique <- unique(pData(eset)[order(pData(eset)$label),c("label", "cols.pca", "cols.legend")])
rownames(ct.unique) <- ct.unique$label
ct.unique <- ct.unique[c("Endo", "Ast", "qNSC", "aNSC", "NPC"),]

pic.size <- 2.5

# Endo

eset.part <- eset[,pData(eset)$label == "Endo"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]

pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:6]

pdf("Results/pca_CD31_sv_corrected.pdf", width = pic.size, height = pic.size, bg = "white")
par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 1.1))
plot(pca$x, pch = pData(eset.part)$pch, col = pData(eset.part)$cols.pca, bg = pData(eset.part)$cols.pca, cex=2,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
     xlim = c(-150, 150),
     ylim = c(-150, 150))
dev.off()


# Ast

eset.part <- eset[,pData(eset)$label == "Ast"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]

pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:6]

pdf("Results/pca_G_sv_corrected.pdf", width = pic.size, height = pic.size, bg = "white")
par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 1.1))
plot(pca$x, pch = pData(eset.part)$pch, col = pData(eset.part)$cols.pca, bg = pData(eset.part)$cols.pca, cex=2,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
     xlim = c(-150, 150),
     ylim = c(-150, 150))
dev.off()


# qNSC

eset.part <- eset[,pData(eset)$label == "qNSC"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]

pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:5]

pdf("Results/pca_PG_sv_corrected.pdf", width = pic.size, height = pic.size, bg = "white")
par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 1.1))
plot(pca$x, pch = pData(eset.part)$pch, col = pData(eset.part)$cols.pca, bg = pData(eset.part)$cols.pca, cex=2,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
     xlim = c(-150, 150),
     ylim = c(-150, 150))
dev.off()


# aNSC

eset.part <- eset[,pData(eset)$label == "aNSC"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]

pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:5]

pdf("Results/pca_PGE_sv_corrected.pdf", width = pic.size, height = pic.size, bg = "white")
par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 1.1))
plot(pca$x, pch = pData(eset.part)$pch, col = pData(eset.part)$cols.pca, bg = pData(eset.part)$cols.pca, cex=2,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
     xlim = c(-150, 150),
     ylim = c(-150, 150))
dev.off()


# NPC

eset.part <- eset[,pData(eset)$label == "NPC"]
ind.var <- apply(exprs(eset.part), 1, function(x) ifelse(length(unique(x)) > 1, FALSE, TRUE))
eset.part <- eset.part[!ind.var,]

pca <- prcomp(t(exprs(eset.part)), scale = TRUE)
s <- summary(pca)$importance[, 1:5]

pdf("Results/pca_E_sv_corrected.pdf", width = pic.size, height = pic.size, bg = "white")
par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 1.1))
plot(pca$x, pch = pData(eset.part)$pch, col = pData(eset.part)$cols.pca, bg = pData(eset.part)$cols.pca, cex=2,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""),
     xlim = c(-150, 150),
     ylim = c(-150, 150))
dev.off()
