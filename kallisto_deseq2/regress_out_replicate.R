########  I decided to go with the vst normalized data where I regress out the replicate.

library(Biobase)
library(limma)
load("Data/eset_vsd.RData")
full.model <- model.matrix(~ label + label:age + replicate, data = pData(eset)) # all variables
fit <- lmFit(eset, full.model)
fit.eb <- eBayes(fit)

### Regress out replicate:

mod <- coefficients(fit)[,-c(1:5, 9:13)] %*% t(fit$design[,-c(1:5, 9:13)])
eset.corrected <- eset
exprs(eset.corrected) <- exprs(eset) - mod
exprs(eset) <- exprs(eset.corrected)

save(eset, file = "Data/eset_vsd_replicate_corrected.RData")


### Write out replicate corrected vst values to add it as a supplementary table:

write.table(exprs(eset),
            file = "Data/vsd_replicate_corrected.txt",
            sep = "\t",
            col.names = NA,
            quote = FALSE)
