tab <- read.table("Rossi_data/Data/GEXC_Normalized_Signal_Intensity.csv", header = TRUE, sep =",")
rownames(tab) <- tab[,1]

exprs.m <- tab[,-c(1,2)]

exprs.means <- by(exprs.m, tab$Gene.Symbol, function(x){
    colMeans(x)
})
exprs.means.m <- do.call("rbind", exprs.means)
exprs.means.m <- exprs.means.m[-which(rownames(exprs.means.m) == "---"),]

eset <- ExpressionSet(assayData = exprs.means.m)

s.info <- strsplit(colnames(eset), "_")
age <- sapply(s.info, function(x) x[1])
age[age == "o"] <- "old"
age[age == "y"] <- "young"
age <- as.factor(age)
age <- relevel(age, "young")
replicate <- as.integer(sapply(s.info, function(x) rev(x)[1]))
cell.type <- sapply(s.info, function(x){
    y = x[-1]
    z = rev(rev(y)[-1])
    return(paste(z, collapse = "_"))
})

fData(eset) <- data.frame(gene = rownames(eset))
pData(eset) <- data.frame(sample = colnames(eset),
                          cell.type = as.factor(cell.type),
                          age = age,
                          replicate = replicate)
pData(eset)$label <- as.character(pData(eset)$cell.type)
pData(eset)$label[pData(eset)$cell.type == "MPP_Flk2_neg"] <- "ST_HSC"
pData(eset)$label[pData(eset)$cell.type == "HSC"] <- "LT_HSC"
pData(eset)$label[pData(eset)$cell.type == "MPP_Flk2_pos"] <- "MPP"
pData(eset)$label <- as.factor(pData(eset)$label)
colnames(eset) <- paste(pData(eset)$label, pData(eset)$age, pData(eset)$replicate, sep = "_")

eset.rossi <- eset
load("Data/log2fpkm_eset.RData")

pData(eset.rossi)$equivalent[pData(eset.rossi)$label == "LT_HSC"] <- "qNSC"
pData(eset.rossi)$equivalent[pData(eset.rossi)$label == "ST_HSC"] <- "aNSC"
pData(eset.rossi)$equivalent[pData(eset.rossi)$label == "MPP"] <- "NPC"

col.info <- pData(eset)[,c("label", "age", "color", "alpha", "pch", "lty.s", "cols", "cols.pca", "cols.legend")]
colnames(col.info)[1] <- "equivalent"
col.info <- unique(col.info)

me <- merge(pData(eset.rossi), col.info, all.x = TRUE, all.y = FALSE)
ind.ord <- match(pData(eset.rossi)$sample, as.character(me$sample))
me <- me[ind.ord,]
rownames(me) <- me$sample
pData(eset.rossi) <- me
eset <- eset.rossi

pData(eset)$cols.pca[pData(eset)$label == "MEP"] <- "red"
pData(eset)$cols.pca[pData(eset)$label == "GMP"] <- "cornflowerblue"
pData(eset)$cols.pca[pData(eset)$label == "CLP"] <- "darkgreen"
pData(eset)$cols.pca[pData(eset)$label == "BLP"] <- "pink"
pData(eset)$cols.pca[pData(eset)$label == "pre_pro.B"] <- "orange"
pData(eset)$cols.pca[pData(eset)$label == "CMP"] <- "brown"
pData(eset)$pch[pData(eset)$age == "young"] <- 21
pData(eset)$pch[pData(eset)$age == "old"] <- 24

save(eset, file = "Data/eset.RData")
