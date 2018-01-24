library(edgeR)

files <- list.files("/srv/gsfs0/projects/brunet/Detja/Data/Gene_counts",
                    recursive = FALSE, pattern = "counts",
                    full.names=TRUE)

labels <- gsub("/srv/gsfs0/projects/brunet/Detja/Data/Gene_counts/counts_",
               	"", files)
labels <- gsub(".txt", "", labels)

rc <- readDGE(files, labels = labels)

out1 <- which(colnames(rc) == "PGE_old_E")
rc <- rc[,-out1]

rc$samples$cell.type <- as.factor(sapply(strsplit(colnames(rc), "_"), function(x) x[1]))
rc$samples$age <- as.factor(sapply(strsplit(colnames(rc), "_"), function(x) x[2]))
rc$samples$replicate <-	as.factor(sapply(strsplit(colnames(rc), "_"), function(x) x[3]))

save(rc, file = "Data/rc.RData")

