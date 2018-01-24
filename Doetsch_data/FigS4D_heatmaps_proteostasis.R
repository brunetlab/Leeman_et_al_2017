## Make three (small) heatmaps that are attached underneath each other in Adobe Illustrator: proteasome, lysosome, unfolded protein binding.

library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(Hmisc)

load("Data/eset.RData")

## Lysosome, Proteasome, GO unfolded protein binding
genes1 <- sort(toupper(read.table("../kallisto_deseq2/Data/Proteostasis_gene_lists/KEGG_Lysosome.txt",
                                  stringsAsFactors = FALSE)[,1]))

genes2 <- sort(toupper(read.table("../kallisto_deseq2/Data/Proteostasis_gene_lists/KEGG_Proteasome.txt",
                                  stringsAsFactors = FALSE)[,1]))

genes3 <- sort(toupper(read.table("../kallisto_deseq2/Data/Proteostasis_gene_lists/Chaperones_Dena_Oct22_2015.txt",
                                  stringsAsFactors = FALSE,
                                  skip = 1)[,1]))


eset.part <- eset[is.element(toupper(fData(eset)$gene), genes1),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)

pdf("Results/heatmap_Lysosome.pdf", width = 3, height = 2)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

eset.part <- eset[is.element(toupper(fData(eset)$gene), genes2),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)

pdf("Results/heatmap_Proteasome.pdf", width = 3, height = 0.9)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

eset.part <- eset[is.element(toupper(fData(eset)$gene), genes3),]
data <- exprs(eset.part)
colnames(data) <- paste(pData(eset.part)$label, pData(eset)$age)

pdf("Results/heatmap_chaperones.pdf", width = 3, height = 2)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()
