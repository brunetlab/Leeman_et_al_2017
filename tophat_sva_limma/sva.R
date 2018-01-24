library(sva)
library(limma)
library(Biobase)

## Created in edgeR_normalize.R:
load("Data/log2fpkm_eset.RData")

full.model <- model.matrix(~ cell.type + cell.type:age, data = pData(eset)) # all variables
null.model <- model.matrix(~ 1, data = pData(eset)) # only adjustment variables
svobj <- sva(exprs(eset),
             mod = full.model,
             mod0 = null.model)

full.model.sv <- cbind(full.model, svobj$sv)

fit <- lmFit(eset, full.model.sv)
fit.eb <- eBayes(fit)



# effect old in CD31
tt.CD31 <- topTable(fit.eb, coef = "cell.typeCD31:ageold", number = Inf)
write.table(data.frame(gene=rownames(tt.CD31), tt.CD31),
            file = "Results/sva_y_vs_o_CD31.csv",
            sep = "\t", row.names=FALSE)


# effect old in E, NPC
tt.E <- topTable(fit.eb, coef = "cell.typeE:ageold", number = Inf)
write.table(data.frame(gene=rownames(tt.E), tt.E),
            file = "Results/sva_y_vs_o_E.csv",
            sep = "\t", row.names=FALSE)

# effect old in G, Ast
tt.G <- topTable(fit.eb, coef = "cell.typeG:ageold", number = Inf)
write.table(data.frame(gene=rownames(tt.G), tt.G),
            file = "Results/sva_y_vs_o_G.csv",
            sep = "\t", row.names=FALSE)

# effect old in PG, qNSC
tt.PG <- topTable(fit.eb, coef = "cell.typePG:ageold", number=Inf)
write.table(data.frame(gene=rownames(tt.PG), tt.PG),
            file = "Results/sva_y_vs_o_PGC.csv",
            sep = "\t", row.names=FALSE)

# effect old in PGE, aNSC
tt.PGE <- topTable(fit.eb, coef = "cell.typePGE:ageold", number = Inf)
write.table(data.frame(gene=rownames(tt.PGE), tt.PGE),
            file = "Results/sva_y_vs_o_PGE.csv",
            sep = "\t", row.names=FALSE)

## Get clean data for PCA:
# Regress out surrogate variables (keep variance caused by age and cell type, which are columns 2 to 10, and intercept, which is column 1):
mod.age.celltype <- coefficients(fit)[,-c(1:10)] %*% t(fit$design[,-c(1:10)])
eset.age.celltype <- eset
exprs(eset.age.celltype) <- exprs(eset) - mod.age.celltype
save(eset.age.celltype, file = "Data/sva_eset_no_age_no_celltype.RData")
