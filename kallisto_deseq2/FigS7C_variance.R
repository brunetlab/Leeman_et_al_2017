library(Biobase)
library(edgeR)

load("Data/eset_vsd.RData")

cts <- c("Endo", "Ast", "qNSC", "aNSC", "NPC")
cv.list <- list()
sd.list <- list()
for(i in 1:length(cts)){
    ct <- cts[i]
    ind.ct <- pData(eset)$label == ct
    ind.old <- pData(eset)$age == "old"
    sd.old <- apply(exprs(eset)[, which(ind.old & ind.ct)], 1, sd)
    cv.old <- sd.old / rowMeans(exprs(eset)[, which(ind.old & ind.ct)])
    sd.young <-  apply(exprs(eset)[, which(!ind.old & ind.ct)], 1, sd)
    cv.young <- sd.young / rowMeans(exprs(eset)[, which(!ind.old & ind.ct)])
    cv.list[[i]] <- data.frame(cv.young = cv.young, cv.old = cv.old)
    sd.list[[i]] <- data.frame(sd.young = sd.young, sd.old = sd.old)
}
names(cv.list) <- cts
names(sd.list) <- cts

cv.df <- do.call("cbind", cv.list)
sd.df <- do.call("cbind", sd.list)
cols <- sapply(cts, function(x, y=pData(eset)){
    c(young = y$cols.pca[(y$label == x) & (y$age == "young")][1],
      old = y$cols.pca[(y$label == x) & (y$age == "old")][1])
})


load("Data/eset_vsd_replicate_corrected.RData")

cts <- c("Endo", "Ast", "qNSC", "aNSC", "NPC")
cv.list <- list()
sd.list <- list()
for(i in 1:length(cts)){
    ct <- cts[i]
    ind.ct <- pData(eset)$label == ct
    ind.old <- pData(eset)$age == "old"
    sd.old <- apply(exprs(eset)[, which(ind.old & ind.ct)], 1, sd)
    cv.old <- sd.old / rowMeans(exprs(eset)[, which(ind.old & ind.ct)])
    sd.young <-  apply(exprs(eset)[, which(!ind.old & ind.ct)], 1, sd)
    cv.young <- sd.young / rowMeans(exprs(eset)[, which(!ind.old & ind.ct)])
    cv.list[[i]] <- data.frame(cv.young = cv.young, cv.old = cv.old)
    sd.list[[i]] <- data.frame(sd.young = sd.young, sd.old = sd.old)
}
names(cv.list) <- cts
names(sd.list) <- cts

cv.df <- do.call("cbind", cv.list)
sd.df <- do.call("cbind", sd.list)
cols <- sapply(cts, function(x, y=pData(eset)){
    c(young = y$cols.pca[(y$label == x) & (y$age == "young")][1],
      old = y$cols.pca[(y$label == x) & (y$age == "old")][1])
})


pdf("Results/variance_across_cell_types_cv_replicate_corrected.pdf", width = 4, height = 4)
par(mar=c(2.1, 4.1, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 5.5),
     ylim = c(0, 0.3),
     axes = FALSE,
     xlab = "",
     ylab = "CV per gene"
     )
abline(h = seq(0, 0.3, by = 0.05),
       lty = "dotted",
       col = "grey")
boxplot(cv.df[seq(1, ncol(cv.df), 2)],
        at = 1:5 - 0.2,
        boxwex = 0.25,
        xaxt = "n",
        outline = FALSE,
        col = cols["young",],
        notch = TRUE,
        axes = FALSE,
        add = TRUE
)
boxplot(cv.df[seq(2, ncol(cv.df), 2)],
        at = 1:5 + 0.2,
        boxwex = 0.25,
        xaxt = "n",
        outline = FALSE,
        axes = FALSE,
        col = cols["old",],
        notch = TRUE,
        add = TRUE)
axis(1,
     at = 1:5,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = cts)
axis(2,
     las = 1,
     at = seq(0, 0.3, by = 0.05))
dev.off()
