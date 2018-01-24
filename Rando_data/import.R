library(GEOquery)
library(Biobase)
eset <- getGEO("GSE47177")
eset <- eset[[1]]

t <- pData(eset)[,c("title",
                    "characteristics_ch1",
                    "characteristics_ch1.1",
                    "characteristics_ch1.2",
                    "characteristics_ch1.3",
                    "characteristics_ch1.4",
                    "description")]
pData(eset)$cell.type <- NA
pData(eset)$age <- NA
pData(eset)[pData(eset)$description == "Gene expression data from VCAM+ quiescent satellite cells isolated from hindlimb muscles of uninjured 2-3-month old mice", c("cell.type", "age")] <- rep(c("qSC", "young"), each = 3)
pData(eset)[pData(eset)$description == "Gene expression data from VCAM+ cells isolated from hindlimb muscles of 2-3-month old mice 60 hours after BaCl2 injury", c("cell.type", "age")] <- rep(c("aSC", "young"), each = 3)
pData(eset)[pData(eset)$description == "Gene expression data from VCAM+ quiescent satellite cells isolated from hindlimb muscles of uninjured 22-24-month old mice", c("cell.type", "age")] <- rep(c("qSC", "old"), each = 3)
pData(eset)[pData(eset)$description == "Gene expression data from VCAM+ cells isolated from hindlimb muscles of 22-24-month old mice 60 hours after BaCl2 injury", c("cell.type", "age")] <- rep(c("aSC", "old"), each = 3)
pData(eset)[!is.na(pData(eset)$age), c("cell.type","age", "description")]


pData(eset)$age <- relevel(as.factor(pData(eset)$age), ref="young")
pData(eset)$cell.type <- relevel(as.factor(pData(eset)$cell.type), ref="qSC")

pData(eset)$pch <- rep(21, ncol(eset)) # circle
pData(eset)$pch[pData(eset)$age == "old"] <- 24 # triangle

pData(eset)$color[pData(eset)$cell.type == "qSC"] <- "darkcyan"
pData(eset)$color[pData(eset)$cell.type == "aSC"] <- "mediumblue"

cols.pca <- character(length=ncol(eset))
for(i in 1:length(cols.pca)){
    col.i <- col2rgb(pData(eset)$color[i])
    cols.pca[i] <- rgb(col.i[1,], col.i[2,], col.i[3,], alpha = 150, maxColorValue = 255)
}
pData(eset)$cols.pca <- cols.pca
cols.legend <- character(length=ncol(eset))
for(i in 1:length(cols.pca)){
    col.i <- col2rgb(pData(eset)$color[i])
    cols.legend[i] <- rgb(col.i[1,], col.i[2,], col.i[3,], alpha = 255, maxColorValue = 255)
}
pData(eset)$cols.legend <- cols.legend

ind.gene <- fData(eset)$gene_assignment != "---"

eset <- eset[ind.gene,]

ga <- as.character(fData(eset)$gene_assignment)
ga1 <- strsplit(ga, " // ")
fData(eset)$gene.symbol <- sapply(ga1, function(x) x[2])
length(unique(fData(eset)$gene.symbol))

eset <- eset[!is.na(fData(eset)$gene.symbol),]


length(unique(toupper(fData(eset)$gene.symbol))) == length(unique(fData(eset)$gene.symbol))
## TRUE

## median per gene:
genes <- unique(as.character(fData(eset)$gene.symbol))
exprs.new <- matrix(numeric(), nrow = length(genes), ncol = ncol(eset))
dimnames(exprs.new) <- list(genes, colnames(eset))
for( i in seq(along = genes)){
    gene.i <- genes[i]
    exprs.i <- exprs(eset)[fData(eset)$gene.symbol == gene.i,, drop = FALSE]
    if(nrow(exprs.i) > 1){
        exprs.i.new <- apply(exprs.i, 2, median)
        exprs.new[i,] <- exprs.i.new
    } else {
        exprs.new[i,] <- exprs.i
    }
}

dim(exprs.new)


eset.new <- ExpressionSet(assayData = exprs.new,
                          phenoData = phenoData(eset))
fData(eset.new)$gene <- rownames(eset.new)
eset <- eset.new


save(eset, file = "Data/eset.RData")
