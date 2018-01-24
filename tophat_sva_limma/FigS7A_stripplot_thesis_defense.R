## Dena used my script and changed colors.
## I want that the significant genes have the colors of the cell types (with colors chosen by Dena for projector).
library(Biobase)

## Created in sva.R:
age.results <- list(Endo = read.table("Results/sva_y_vs_o_CD31.csv", header = TRUE),
                    Ast = read.table("Results/sva_y_vs_o_G.csv", header = TRUE),
                    qNSC = read.table("Results/sva_y_vs_o_PG.csv", header = TRUE),
                    aNSC = read.table("Results/sva_y_vs_o_PGE.csv", header = TRUE),
                    NPC = read.table("Results/sva_y_vs_o_E.csv", header = TRUE))

## for colors:
load("Data/log2fpkm_eset.RData")
ct.unique <- unique(pData(eset)[order(pData(eset)$label),c("label", "cols.pca", "cols.legend")])
#ct.unique cols.pca (transparent)   #### TO SEE COLORS ####
rownames(ct.unique) <- ct.unique$label
ct.unique <- ct.unique[c("Endo", "Ast", "qNSC", "aNSC", "NPC"),]

ct.unique$cols.pca
# original colors   "#FF634796" "#9ACD3296" "#008B8B96" "#0000CD96" "#8B478996"

col2rgb("#9ACD3296", alpha=T)

#############################################################################
ct.unique #cols.pca is the (transparent) column  #### TO SEE COLORS ####
# color picking
col2rgb("red") # to find out a color's RGB values from its name
rgb(1, 1, 1, alpha = 0.5) # gives it a hash name
rgb(255, 0, 0, maxColorValue = 255, alpha = 100)
#############################################################################

#Endo
col2rgb("tomato2") # to find out a color's RGB values from its name
ct.unique$cols.pca[1]<-rgb(230, 102, 76, maxColorValue = 255, alpha = 150)

#Ast
#col2rgb("greenyellow") # to find out a color's RGB values from its name
#rgb(154, 205, 50, maxColorValue = 255, alpha = 110) # gives it a hash name
ct.unique$cols.pca[2]<-rgb(150, 200, 50, maxColorValue = 255, alpha = 150)

#qNSC
col2rgb("darkcyan") # to find out a color's RGB values from its name
#rgb(0, 139, 139, maxColorValue = 255, alpha = 150) # gives it a hash name
ct.unique$cols.pca[3]<-rgb(10, 169, 149, maxColorValue = 255, alpha = 150)


## Remove rows with NA p-values:
age.results <- lapply(age.results,
                      function(x) {
                          x[!is.na(x$adj.P.Val), ]
                      })
n <- sapply(age.results, nrow)
names(n) <- names(age.results)

## Order by pvalue:
age.results <- lapply(age.results,
                      function(x) {
                          x[order(x$adj.P.Val),]
                      })

p <- 0.05
cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
    cols[[i]] <- rep(rgb(163, 163, 163, maxColorValue = 255, alpha = 50), n[i]) # grey60
    ind.sig.i <- age.results[[i]]$adj.P.Val < p
    cols[[i]][ind.sig.i] <- ct.unique[names(age.results)[i], "cols.pca"]
    xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)


pdf("Results/stripplot_thesis_defense.pdf", width = 5, height = 3)
par(mar = c(3.1, 4.1, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 5.5),
     ylim = c(-5, 7),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change (old / young)"
     )
abline(h = 0)
abline(h = seq(-5, 7, by = 1)[-6],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
    set.seed(1234)
    points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
           y = rev(age.results[[i]]$logFC),
           cex = 0.9,
           pch = 21,
           col = NA,
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
     at = seq(-5, 7, by = 1))
dev.off()
