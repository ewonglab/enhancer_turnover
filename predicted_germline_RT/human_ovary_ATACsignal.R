Replicon for replication time simulation - from human testis and ovary ATAC-seq

Mouse test


cd /g/data/zk16/cc3704/replication_timing
export PATH=/g/data/zk16/xzhang/bin:$PATH
mouse_path="/g/data/zk16/cc3704/replication_timing/mouse/testis_DNAseq_mm10_norm/chr1/"
human_path="/g/data/zk16/cc3704/replication_timing/human/ENCFF289LWU_ATAC_testis_hg19_IPLS/chr1/"
replicon -threads 8 ${mouse_path}chr1_ENCFF728ZKQ_testis_DNAseq_mm10_small_norm chr1_mm10_norm

Human

library(data.table)
testis <- fread("ENCFF289LWU_ATAC_testis_hg38.bed")
testis[,7] <- scale(x=testis[,7], center=T, scale=T)
ovary <- fread("ENCFF586QDS_ovary_hg38.bed.gz")
ovary[,7] <- scale(x=ovary[,7], center=T, scale=T)
write.table(x= testis[,c(1:3,7)], file = "ENCFF289LWU_ATAC_testis_hg38_scaled.bed", col.names=F, row.names=F, sep="\t", quote=F)
write.table(x= ovary[,c(1:3,7)], file = "ENCFF586QDS_ATAC_ovary_hg38_scaled.bed", col.names=F, row.names=F, sep="\t", quote=F)



