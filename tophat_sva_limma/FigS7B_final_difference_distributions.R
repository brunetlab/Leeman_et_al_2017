library(edgeR)
library(Biobase)

## Created in sva.R
age.results <- list(Endo = read.table("Results/sva_y_vs_o_CD31.csv", header = TRUE),
                    Ast = read.table("Results/sva_y_vs_o_G.csv", header = TRUE),
                    qNSC = read.table("Results/sva_y_vs_o_PG.csv", header = TRUE),
                    aNSC = read.table("Results/sva_y_vs_o_PGE.csv", header = TRUE),
                    NPC = read.table("Results/sva_y_vs_o_E.csv", header = TRUE))

## For colors:
load("Data/log2fpkm_eset.RData")

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
     ylim = c(0, 2),
     axes = FALSE,
     xlab = "",
     ylab = "Absolute log2 fold change (old / young)"
     )
abline(h = seq(0, 2, by = 0.25),
       lty = "dotted",
       col = "grey")
boxplot(list(Endo = abs(age.results$Endo$logFC),
             Ast = abs(age.results$Ast$logFC),
             qNSC = abs(age.results$qNSC$logFC),
             aNSC = abs(age.results$aNSC$logFC),
             NPC = abs(age.results$NPC$logFC)),
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
     at = seq(0, 2, by = 0.25))
dev.off()
