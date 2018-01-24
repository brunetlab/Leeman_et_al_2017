library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(Hmisc)
load("Data/eset.RData")

eset <- eset[, is.element(pData(eset)$cell.type, c("qSC", "aSC"))]
eset <- eset[, pData(eset)$age == "young"]
rownames(eset) <- toupper(rownames(eset))

## Lysosome, Proteasome, GO unfolded protein binding
genes1 <- sort(toupper(read.table("../kallisto_deseq2/Data/Proteostasis_gene_lists/KEGG_Lysosome.txt",
                                  stringsAsFactors = FALSE)[,1]))
genes2 <- sort(toupper(read.table("../kallisto_deseq2/Data/Proteostasis_gene_lists/KEGG_Proteasome.txt",
                                  stringsAsFactors = FALSE)[,1]))
genes3 <- sort(toupper(read.table("../kallisto_deseq2/Data/Proteostasis_gene_lists/Chaperones_Dena_Oct22_2015.txt",
                                  stringsAsFactors = FALSE,
                                  skip = 1)[,1]))

eset.part <- eset[is.element(rownames(eset), genes1),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)
pdf("Results/heatmap_Lysosome.pdf", width = 2, height = 2)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

eset.part <- eset[is.element(rownames(eset), genes2),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)
pdf("Results/heatmap_Proteasome.pdf", width = 2, height = 0.9)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

eset.part <- eset[is.element(rownames(eset), genes3),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)
pdf("Results/heatmap_chaperones.pdf", width = 2, height = 2)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()
