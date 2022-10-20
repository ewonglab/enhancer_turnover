library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library(data.table)
library(reshape2)


setwd("/g/data/zk16/veronika/projects/repli_time/AML/ATAC/peaks/")
SU484 <- fread("aml_SU484_1kb.bed", header = F)
SU501 <- fread("aml_SU501_1kb.bed", header = F)
SU654 <- fread("aml_SU654_1kb.bed", header = F)

setwd("/g/data/zk16/cc3704/replication_timing/human/erythroid")

SU484_gains <- SU484[SU484$V4=="gain",]
SU484_unch <- SU484[SU484$V4=="unchanged",]
SU484_lose <- SU484[SU484$V4=="lose",]

SU501_gains <- SU501[SU501$V4=="gain",]
SU501_unch <- SU501[SU501$V4=="unchanged",]
SU501_lose <- SU501[SU501$V4=="lose",]

SU654_gains <- SU654[SU654$V4=="gain",]
SU654_unch <- SU654[SU654$V4=="unchanged",]
SU654_lose <- SU654[SU654$V4=="lose",]


#make gr objects
SU484_gains.gr <- with(SU484_gains, GRanges(V1 , IRanges( V2+1, V3 )))
SU484_unch.gr <- with(SU484_unch, GRanges(V1 , IRanges( V2+1, V3 )))
SU484_lose.gr <- with(SU484_lose, GRanges(V1 , IRanges( V2+1, V3 )))

SU501_gains.gr <- with(SU501_gains, GRanges(V1 , IRanges( V2+1, V3 )))
SU501_unch.gr <- with(SU501_unch, GRanges(V1 , IRanges( V2+1, V3 )))
SU501_lose.gr <- with(SU501_lose, GRanges(V1 , IRanges( V2+1, V3 )))

SU654_gains.gr <- with(SU654_gains, GRanges(V1 , IRanges( V2+1, V3 )))
SU654_unch.gr <- with(SU654_unch, GRanges(V1 , IRanges( V2+1, V3 )))
SU654_lose.gr <- with(SU654_lose, GRanges(V1 , IRanges( V2+1, V3 )))


#get union of gains, loses and unchanged enhancers

tumor_peaks <- GenomicRanges::union(SU484_gains.gr, c(SU501_gains.gr, SU654_gains.gr))
normal_peaks <- GenomicRanges::union(SU484_lose.gr, c(SU501_lose.gr, SU654_lose.gr))
common_peaks <- GenomicRanges::union(SU484_unch.gr, c(SU501_unch.gr, SU654_unch.gr))

tumor_peaks.df <- data.frame(chr=seqnames(tumor_peaks), start=start(tumor_peaks)-1
                             , end=end(tumor_peaks), stringsAsFactors = F)
normal_peaks.df <- data.frame(chr=seqnames(normal_peaks), start=start(normal_peaks)-1
                              , end=end(normal_peaks), stringsAsFactors = F)
common_peaks.df <- data.frame(chr=seqnames(common_peaks), start=start(common_peaks)-1
                              , end=end(common_peaks), stringsAsFactors = F)

write.table(x = tumor_peaks.df, file = "gains_AML_union_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
write.table(x = normal_peaks.df, file = "loses_AML_union_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
write.table(x = common_peaks.df, file = "unchanged_AML_union_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
