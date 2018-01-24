## This function tests for enrichment of KEGG pathways in a gene list.
## The null hypothesis is per gene set is: "The genes in this gene set show the same pattern of association with the phenotype compared with the rest of the genes."
## As a test statistic per gene set I will use the average (absolute) test statistic / standardized correlation coefficient etc.
## See Ackermann & Strimmer 2009.

## INPUT:

## gene.sets.file: file name; default: Data/c2.cp.kegg.v5.0.symbols.gmt (or any other file in .gmt format, available from MSigDB: http://www.broadinstitute.org/gsea/downloads.jsp), should be all UPPERCASE letters!

## gene.sets.list: alternative to gene.sets.file. This is a list where names of the list elements are the gene set names. Each list element is a character vector with the genes in this gene set. Should be all UPPERCASE letters!

## num.genes: minimum number of genes a pathway should have; pathways with fewer than num.genes are removed from the analysis.

## test.results: data.frame with test results, row.names: gene names (UPPERCASE letters!)

## statistic: corresponds to the column in test.results that contains the gene-wise test statistics that should be used, default: "t" (t-test statistic)

## fold.change: corresponds to the column in test.results that contains the gene-wise fold changes (the average fold changes can be used plotting enrichments)

## num.samples: how often should gene sets be sampled in order to get a p-value? default: 1,000

## abs: If abs = FALSE, the direction of the statistic plays a role (similar to GSEA), if TRUE, the direction of the statistic does not play a role and the absolute values of the gene-wise statistics are used.

## HEADS UP: If p values are reported as 0, that just means that they are actually p < 1/num.samples and should be reported like this.

KEGG.pathway.enrichment <- function(gene.sets.file = "Data/c2.cp.kegg.v5.0.symbols.gmt",
                                    gene.sets.list = NULL,
                                    num.genes = 20,
                                    test.results,
                                    statistic = "t",
                                    fold.change = "log2FoldChange",
                                    num.samples = 1000,
                                    abs = FALSE) {
    if(is.null(gene.sets.list)){
        ## Import gmt file and create list of gene sets:
        gene.sets <- read.table(gene.sets.file,
                                sep = "\n",
                                quote = "",
                                stringsAsFactors = FALSE)
        gene.sets.l <- apply(gene.sets, 1,
                             function(x){
                                 y <- strsplit(x, "\t")
                                 return(as.character(unlist(y)))
                             })
        names(gene.sets.l) <- sapply(gene.sets.l, function(x) x[1])
        ## Remove empty strings:
        gene.sets.l <- lapply(gene.sets.l,
                              function(x) {
                                  ind.empty <- x == ""
                                  return(x[!ind.empty])
                              })

    } else {
        gene.sets.l <- gene.sets.list
    }
    ## Remove genes in gene lists that were not tested:
    gene.sets.l.red <- lapply(gene.sets.l, function(x){
        intersect(x, rownames(test.results))
    })
    ## Remove gene lists with fewer than num.genes:
    ind.too.small <- sapply(gene.sets.l.red, length) < num.genes
    gene.sets.l.red <- gene.sets.l.red[!ind.too.small]
    ## Calculate average of squared statistics per gene set and return all genes that are more extreme than the average, also calculate average of fold changes:
    test.statistics <- lapply(gene.sets.l.red,
                              function(x){
                                  gene.statistics <- test.results[x, statistic, drop = FALSE]
                                  gene.fold.changes <- test.results[x, fold.change, drop = FALSE]
                                  if(abs){
                                      avg.fold.change <- mean(abs(gene.fold.changes[, fold.change]))
                                      m <- mean(abs(gene.statistics[, statistic]))
                                      ind.driver.genes <- abs(gene.statistics[, statistic]) > m
                                      driver.genes <- rownames(gene.statistics)[ind.driver.genes]
                                  }
                                  if(!abs){
                                      avg.fold.change <- mean(gene.fold.changes[, fold.change])
                                      m <- mean(gene.statistics[, statistic])
                                      sign <- sign(m)
                                      if(sign == 1) {
                                          ind.driver.genes <- gene.statistics[, statistic] > m
                                      }
                                      if(sign == -1) {
                                          ind.driver.genes <- gene.statistics[, statistic] < m
                                      }
                                      if(sign == 0) {
                                          ind.driver.genes <- rep(FALSE, nrow(gene.statistics))
                                      }
                                      driver.genes <- rownames(gene.statistics)[ind.driver.genes]
                                  }
                                  return(list(statistic = m,
                                              fold.change = avg.fold.change,
                                              driver.genes = driver.genes))
                              }
                              )
    gene.set.lengths <- sapply(gene.sets.l.red, length)
    ## For each gene set: Sample a gene set of the same size, calculate averages of (squared) statistics and calculate p-value
    p.values <- numeric(length = length(test.statistics))
    names(p.values) <- names(gene.set.lengths)
    for(i in seq(along = gene.set.lengths)){
        gene.statistics.i <- sample(test.results[, statistic],
                                    size = gene.set.lengths[i] * num.samples,
                                    replace = TRUE)
        ## per sample, one column:
        gene.statistics.i.matrix <- matrix(gene.statistics.i, ncol = num.samples)
        if(abs){
            test.statistic.i <- colMeans(abs(gene.statistics.i.matrix))
        }
        if(!abs){
            test.statistic.i <- colMeans(gene.statistics.i.matrix)
        }
        ## How often is sampled gene set statistic equal or greater than the real gene set statistic? -> the percentage is the p-value
        p.values[i] <- sum(abs(test.statistic.i) >= abs(test.statistics[[i]]$statistic)) / num.samples
    }
    results.df <- data.frame(num.genes = gene.set.lengths,
                             p.val = p.values,
                             adj.p.val = p.adjust(p.values, method = "BH"),
                             test.statistic = sapply(test.statistics, function(x) x$statistic),
                             fold.change = sapply(test.statistics, function(x) x$fold.change))
    driver.genes <- lapply(test.statistics,
                           function(x){
                               y <- test.results[x$driver.genes, ]
                               y <- y[order(abs(y[, statistic]), decreasing = TRUE),]
                           })
    ind.ord <- order(results.df$p.val)
    results.df <- results.df[ind.ord, ]
    driver.genes <- driver.genes[ind.ord]
    return(list(pathway.results = results.df,
                driver.genes = driver.genes))
}

