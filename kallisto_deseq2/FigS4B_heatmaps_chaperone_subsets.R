## Make three (small) heatmaps that are attached underneath each other in Adobe Illustrator: Endoplasmic_reticulum_unfolded_protein_response_GO0030968, TriC_subunits, Prefoldin_subunits

library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(Hmisc)

load("Data/eset_vsd_replicate_corrected.RData")
eset <- eset[, c(21:24, 17:20, 34:35, 31:33,
                 28:30, 25:27, 13:16, 9:12)]
rownames(eset) <- toupper(rownames(eset))

sets <- read.table("Data/Chaperone_subsets_forpaper.txt",
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   header = TRUE)
gene.sets.list <- list(sets[,1], sets[,2], sets[,3])
names(gene.sets.list) <- colnames(sets)
gene.sets.list <- lapply(gene.sets.list,
                         function(x) x[x != ""])

genes.1 <- gene.sets.list$Endoplasmic_reticulum_unfolded_protein_response_GO0030968
genes.2 <- gene.sets.list$TriC_subunits
genes.3 <- gene.sets.list$Prefoldin_subunits

eset.part <- eset[is.element(rownames(eset), genes.1),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)

pdf("Results/heatmap_Endoplasmic_reticulum_unfolded_protein_response.pdf",
    onefile = FALSE,
    width = 5, height = 2)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

eset.part <- eset[is.element(rownames(eset), genes.2),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)

pdf("Results/heatmap_TriC_subunits.pdf",
    onefile = FALSE,
    width = 5, height = 0.5)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

eset.part <- eset[is.element(rownames(eset), genes.3),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)

pdf("Results/heatmap_Prefoldin_subunits.pdf",
    onefile = FALSE,
    width = 5, height = 0.5)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()
