library(data.table)
setwd("/g/data/zk16/cc3704/replication_timing/human")

thyroid <- fread("ENCFF682VUY_thyroid_ATAC_hg38.bed.gz")
thyroid[,7] <- scale(x=thyroid[,7], center=T, scale=T)
write.table(x = thyroid[,c(1:3,7)], file= "ENCFF682VUY_thyroid_ATAC_hg38_scaled.bed"
            , col.names = F, row.names = F, sep = '\t', quote = F)
