## Created in edgeR_import.R
load("Data/rc.RData")
counts.out <- rc$counts["__no_feature",] + rc$counts["__ambiguous",] + rc$counts["__alignment_not_unique",]
counts.out/colSums(rc$counts)
# around 40% for each sample not in exons

out.1 <- which(colnames(rc) == "PG_young_C")
out.2 <- which(colnames(rc) == "PG_young_E")
out.3 <- which(colnames(rc) == "PG_old_C")
out.4 <- which(colnames(rc) == "PGE_young_C")
genes.out <- which(is.element(rownames(rc), c("__no_feature", "__ambiguous", "__alignment_not_unique")))
rc <- rc[-genes.out, -out.4]

rc$samples$age <- relevel(as.factor(rc$samples$age), ref="young")

design <-  model.matrix(~ cell.type + cell.type:age, data = rc$samples)

# filter out genes with low coverage:
keep <- rowSums(cpm(rc) > 1) >= 2
rc <- rc[keep,]

rownames(rc) <- tolower(rownames(rc))

rc$samples$lib.size <- colSums(rc$counts)

rc.norm <- calcNormFactors(rc)
rc.com <- estimateGLMCommonDisp(rc.norm, design)
rc.trend <- estimateGLMTrendedDisp(rc.com, design)
rc.tag <- estimateGLMTagwiseDisp(rc.trend, design)

load("Data/exonic_gene_sizes.RData")
names(exonic.gene.sizes) <- tolower(names(exonic.gene.sizes))
ind <- match(rownames(rc.norm), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
rc.norm$genes <- data.frame(gene.symbol=rownames(rc.norm),
                            exonic.size=exonic.gene.sizes.ord)

fpkm <- rpkm(rc.norm, gene.length=rc.norm$genes$exonic.size,
             normalized.lib.sizes=TRUE, log=FALSE)
log2fpkm <- log2( fpkm + 1)


library(Biobase)
eset <- ExpressionSet(assayData = log2fpkm,
                      phenoData = AnnotatedDataFrame(rc$samples))
pData(eset)$label[pData(eset)$cell.type == "CD31"] <- "Endo"
pData(eset)$label[pData(eset)$cell.type == "E"] <- "NPC"
pData(eset)$label[pData(eset)$cell.type == "G"] <- "Ast"
pData(eset)$label[pData(eset)$cell.type == "PG"] <- "qNSC"
pData(eset)$label[pData(eset)$cell.type == "PGE"] <- "aNSC"

pData(eset)$color[pData(eset)$label == "Endo"] <- "tomato"
pData(eset)$color[pData(eset)$label == "Ast"] <- "yellowgreen"
pData(eset)$color[pData(eset)$label == "qNSC"] <- "darkcyan"
pData(eset)$color[pData(eset)$label == "aNSC"] <- "mediumblue"
pData(eset)$color[pData(eset)$label == "NPC"] <- "orchid4"

pData(eset)$alpha[pData(eset)$age == "young"] <- 150
pData(eset)$alpha[pData(eset)$age == "old"] <- 255

pData(eset)$pch <- rep(21, ncol(eset)) # circle
pData(eset)$pch[pData(eset)$age == "old"] <- 24 # triangle
pData(eset)$lty.s <- rep(1, ncol(eset))
pData(eset)$lty.s[pData(eset)$age == "old"] <- 2 # dashed

cols <- character(length=ncol(eset))
for(i in 1:length(cols)){
    col.i <- col2rgb(pData(eset)$color[i])
    cols[i] <- rgb(col.i[1,], col.i[2,], col.i[3,], alpha = pData(eset)$alpha[i], maxColorValue = 255)
}
pData(eset)$cols <- cols
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

save(eset, file = "Data/log2fpkm_eset.RData")


