library(edgeR)
library(Biobase)
setwd("~/Detja/Kallisto")

load("Data/age_results_deseq_with_replicate.RData")
age.results.deseq <- lapply(age.results.deseq, as.data.frame)

## For colors:
load("Data/eset_vsd.RData")

cols <- c(unique(pData(eset)$cols.legend[pData(eset)$label == "Endo"]),
          unique(pData(eset)$cols.legend[pData(eset)$label == "Ast"]),
          unique(pData(eset)$cols.legend[pData(eset)$label == "qNSC"]),
          unique(pData(eset)$cols.legend[pData(eset)$label == "aNSC"]),
          unique(pData(eset)$cols.legend[pData(eset)$label == "NPC"]))
names(cols) <- c("Endo", "Ast", "qNSC", "aNSC", "NPC")
lwd = 2

pdf("Results/final_difference_distributions.pdf", width = 3, height = 6)
par(mar = c(3.1, 4.1, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 5.5),
     ylim = c(0, 1.75),
     axes = FALSE,
     xlab = "",
     ylab = "Absolute log2 fold change (old / young)"
     )
abline(h = seq(0, 1.75, by = 0.25),
       lty = "dotted",
       col = "grey")
boxplot(list(Endo = abs(age.results.deseq$Endo$log2FoldChange),
             Ast = abs(age.results.deseq$Ast$log2FoldChange),
             qNSC = abs(age.results.deseq$qNSC$log2FoldChange),
             aNSC = abs(age.results.deseq$aNSC$log2FoldChange),
             NPC = abs(age.results.deseq$NPC$log2FoldChange)),
        col = cols,
        notch = TRUE,
        horizontal = FALSE,
        outline = FALSE,
        boxwex = 0.4,
        axes = FALSE,
        add = TRUE)
axis(1,
     at = 1:5,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = names(cols))
axis(2,
     las = 1,
     at = seq(0, 1.75, by = 0.25))
dev.off()
