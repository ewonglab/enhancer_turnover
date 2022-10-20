
library(data.table)

#read and scale ATAC-seq signal

testis <- fread("ENCFF289LWU_ATAC_testis_hg38.bed")
testis[,7] <- scale(x=testis[,7], center=T, scale=T)
ovary <- fread("ENCFF586QDS_ovary_hg38.bed.gz")
ovary[,7] <- scale(x=ovary[,7], center=T, scale=T)
write.table(x= testis[,c(1:3,7)], file = "ENCFF289LWU_ATAC_testis_hg38_scaled.bed", col.names=F, row.names=F, sep="\t", quote=F)
write.table(x= ovary[,c(1:3,7)], file = "ENCFF586QDS_ATAC_ovary_hg38_scaled.bed", col.names=F, row.names=F, sep="\t", quote=F)



