## Bring all counts from all samples plus annotation of the transcripts into one expressionSet object

library(Biobase)

dirs <- dir("Results")
length(dirs)
## 40
## Do all abundance files exist?
files <- list.files("Results",
                    pattern = ".tsv",
                    recursive = TRUE)
length(files)
## 40
## Yes! All samples could be quantified :)

aligned <- sapply(strsplit(files, "/"), function(x) x[1])
dirs <- aligned

n <- length(dirs)
system("wc -l Results/CD31_old_C/abundance.tsv")
## 54595 Results/CD31_old_C/abundance.tsv
## first row is the header
m <- 54595 - 1
count.matrix <- matrix(numeric(length = n * m),
                       ncol = n,
                       nrow = m)
colnames(count.matrix) <- dirs
tpm.matrix <- matrix(numeric(length = n * m),
                     ncol = n,
                     nrow = m)
colnames(tpm.matrix) <- dirs

for( i in 1:n ){
    data.i <- read.table(paste0("Results/",
                                dirs[i],
                                "/abundance.tsv"),
                         colClasses = c("NULL", "NULL","NULL", "numeric", "numeric"),
                         sep = "\t",
                         header = TRUE)
    count.matrix[, i] <- data.i$est_counts
    tpm.matrix[, i] <- data.i$tpm
    rm(data.i)
    cat(i, "\n")
}


save(count.matrix, tpm.matrix,
     file = "Data/kallisto_results_matrices.RData")



#### Executable: 

load("Data/kallisto_results_matrices.RData")

counts.per.sample <- apply(count.matrix, 2, sum)
summary(counts.per.sample)

anno.tab <- read.table(paste0("Results/CD31_old_C/abundance.tsv"),
                       colClasses = c("character", "integer","NULL", "NULL", "NULL"),
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = FALSE)

anno.list <- strsplit(anno.tab[,1], "|", fixed = TRUE)


## Do all transcripts have a Ensembl transcript ID that appear on first position?
check <- sapply(anno.list, function(x) grep("ENSMUST", x))
all(check == 1)
## TRUE

## Do all transcripts have a Ensembl gene model ID that appears on second position?
check <- sapply(anno.list, function(x) grep("ENSMUSG", x))
all(check == 2)
## TRUE

## Do all transcripts have a VEGA gene model ID that appears on third position?
check <- sapply(anno.list, function(x) grep("OTTMUSG", x))
all(check == 3)
## NA
## That is not great. I have to grep everything.

grep.anno <- function(pattern) {
    sapply(anno.list,
           function(x) {
               ind <- grep(pattern, x)
               if(length(ind) == 1) {
                   x[ind]
               } else {
                   NA
               }
           })
}

## To extract the symbol and the symbol together with transcript number:
symbols <- sapply(anno.list,
                  function(x) {
                      y <- x[x != "-"]
                      excl <- paste(c("ENSMUST", "ENSMUSG", "OTTMUST", "OTTMUSG", "UTR5", "CDS", "UTR3"),
                                    collapse = "|")
                      ## contains symbol transcript with number, symbol transcript and length:
                      y <- y[!grepl(excl, y)]
                      ## gives NAs back if a character is not convertable to a string:
                      y.int <- strtoi(y)
                      y.symbols <- y[is.na(y.int)]
                      if(length(y.symbols) > 0) {
                          return(y.symbols)
                      } else {
                          return(c(NA, NA))
                      }
                  })

anno.df <- data.frame(ensembl.transcript = grep.anno("ENSMUST"),
                      ensembl.gene.model = grep.anno("ENSMUSG"),
                      symbol.transcript.num = symbols[1,],
                      symbol = symbols[2,],
                      length = anno.tab$length,
                      vega.transcript = grep.anno("OTTMUST"),
                      vega.gene.model = grep.anno("OTTMUSG"),
                      utr5 = grep.anno("UTR5"),
                      cds = grep.anno("CDS"),
                      utr3 = grep.anno("UTR3"),
                      stringsAsFactors = FALSE)
rownames(anno.df) <- anno.df$ensembl.transcript

identical(colnames(count.matrix), colnames(tpm.matrix))
## TRUE

cnames <- strsplit(colnames(count.matrix), "_")
sample.anno <- data.frame(cell.type = as.factor(sapply(cnames, function(x) x[1])),
                      age = as.factor(sapply(cnames, function(x) x[2])),
                      replicate = as.factor(sapply(cnames, function(x) x[3])))
sample.anno$label[sample.anno$cell.type == "CD31"] <- "Endo"
sample.anno$label[sample.anno$cell.type == "G"] <- "Ast"
sample.anno$label[sample.anno$cell.type == "PG"] <- "qNSC"
sample.anno$label[sample.anno$cell.type == "PGE"] <- "aNSC"
sample.anno$label[sample.anno$cell.type == "E"] <- "NPC"
sample.anno$label <- as.factor(sample.anno$label)
sample.anno$age <- relevel(sample.anno$age, "young")

sample.anno$color[sample.anno$label == "Endo"] <- "tomato"
sample.anno$color[sample.anno$label == "Ast"] <- "yellowgreen"
sample.anno$color[sample.anno$label == "qNSC"] <- "darkcyan"
sample.anno$color[sample.anno$label == "aNSC"] <- "mediumblue"
sample.anno$color[sample.anno$label == "NPC"] <- "orchid4"

sample.anno$alpha[sample.anno$age == "young"] <- 150
sample.anno$alpha[sample.anno$age == "old"] <- 255

sample.anno$pch <- rep(21, ncol(count.matrix)) # circle
sample.anno$pch[sample.anno$age == "old"] <- 24 # triangle
sample.anno$lty.s <- rep(1, ncol(count.matrix))
sample.anno$lty.s[sample.anno$age == "old"] <- 2 # dashed


cols.pca <- character(length=ncol(count.matrix))
for(i in 1:length(cols.pca)){
    col.i <- col2rgb(sample.anno$color[i])
    cols.pca[i] <- rgb(col.i[1,], col.i[2,], col.i[3,], alpha = 150, maxColorValue = 255)
}
sample.anno$cols.pca <- cols.pca
cols.legend <- character(length=ncol(count.matrix))
for(i in 1:length(cols.pca)){
    col.i <- col2rgb(sample.anno$color[i])
    cols.legend[i] <- rgb(col.i[1,], col.i[2,], col.i[3,], alpha = 255, maxColorValue = 255)
}
sample.anno$cols.legend <- cols.legend

sample.anno$original.id <- colnames(count.matrix)
rownames(sample.anno) <- paste(sample.anno$label, sample.anno$age, sample.anno$replicate, sep = "_")


### Add lane information:

lane.info <- read.table("Data/lane_pooling_BGI.txt",
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        header = TRUE)
lane.info$cell.type <- sapply(strsplit(lane.info$sample.id, "_"), function(x) x[3])
lane.info$age <- sapply(strsplit(lane.info$sample.id, "_"), function(x) x[2])
lane.info$replicate <- sapply(strsplit(lane.info$sample.id, "_"), function(x) x[1])

info <- merge(lane.info, sample.anno,
              sort = FALSE)

info.ind <- match(sample.anno$original.id, info$original.id)
info.ord <- info[info.ind,]
all( sample.anno$original.id == info.ord$original.id)
## TRUE
rownames(info.ord) <- rownames(sample.anno)

sample.anno <- info.ord


## create eset object:

rownames(count.matrix) <- rownames(anno.df)
eset.counts <- ExpressionSet(assayData = count.matrix)
pData(eset.counts) <- sample.anno
fData(eset.counts) <- anno.df
colnames(eset.counts) <- rownames(sample.anno)

save(eset.counts, file = "Data/eset_counts_all_samples.RData")


## Collapse transcripts to genes by just summing the reads of all transcripts per gene:

genes <- unique(fData(eset.counts)$symbol)

exprs.new <- matrix(nrow = length(genes),
                    ncol = ncol(eset.counts),
                    dimnames = list(genes, colnames(eset.counts)))

for(i in seq(along = genes)){
    gene.i <- genes[i]
    transcripts.i <- fData(eset.counts)$ensembl.transcript[fData(eset.counts)$symbol == gene.i]
    exprs.i <- colSums(exprs(eset.counts)[transcripts.i, , drop = FALSE])
    exprs.new[gene.i,] <- exprs.i
}


eset.gene.counts <- ExpressionSet(assayData = exprs.new)
pData(eset.gene.counts) <- pData(eset.counts)

save(eset.gene.counts,
     file = "Data/eset_gene_counts_all_samples.RData")
