library(GEOquery)
eset <- getGEO("GSE54653")
eset <- eset[[1]]

pData(eset) <- pData(eset)[,c("title", "geo_accession")]
pData(eset)$cell.type <- c(rep("qNSC", 3), rep("aNSC", 3))
pData(eset)$replicate <- rep(1:3, 2)

sampleNames(eset) <- paste(pData(eset)$cell.type, pData(eset)$replicate, sep = "_")

pData(eset)$col[pData(eset)$cell.type == "qNSC"] <- "blue"
pData(eset)$col[pData(eset)$cell.type == "aNSC"] <- "red"

pData(eset)$cell.type <- relevel(as.factor(pData(eset)$cell.type), ref = "qNSC")

gene.info <- as.character(fData(eset)[, "Gene Symbol"])
rows.out <- grep( " /// ", gene.info)
eset <- eset[-rows.out,]
rows.out <- which(fData(eset)[, "Gene Symbol"] == "")
eset <- eset[-rows.out,]
fData(eset)$gene <- fData(eset)[, "Gene Symbol"]

## median per gene:
genes <- unique(as.character(fData(eset)$gene))
exprs.new <- matrix(numeric(), nrow = length(genes), ncol = ncol(eset))
dimnames(exprs.new) <- list(genes, colnames(eset))
for( i in seq(along = genes)){
    gene.i <- genes[i]
    exprs.i <- exprs(eset)[fData(eset)$gene == gene.i,, drop = FALSE]
    if(nrow(exprs.i) > 1){
        exprs.i.new <- apply(exprs.i, 2, median)
        exprs.new[i,] <- exprs.i.new
    } else {
        exprs.new[i,] <- exprs.i
    }
}

eset.new <- ExpressionSet(assayData = exprs.new,
                          phenoData = phenoData(eset))
fData(eset.new)$gene <- rownames(eset.new)
eset <- eset.new

save(eset, file = "Data/eset.RData")
