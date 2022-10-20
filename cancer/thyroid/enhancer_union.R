library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library(data.table)
library(reshape2)


setwd("/g/data/zk16/veronika/projects/repli_time/thyroid/ATAC/peaks/")
hg19_txdb <- makeTxDbFromGFF(file="/g/data/zk16/cc3704/human_data/hg19.ensGene.gtf"
                             , format="gtf",organism="Homo sapiens")
hg19_prom <- GenomicFeatures::promoters(x = hg19_txdb
                                        , upstream = 1000
                                        , downstream = 1000)
#Loading macs3 output
chrOrder<-paste0("chr", c(1:22,"X","Y"))

#xsl files from macs are 1-based
dir = "/g/data/zk16/veronika/projects/repli_time/thyroid/ATAC/peaks/"
sample <- c("F001_C1_T1", "F002_C1_N1", "F017_C7_T1", "F018_C7_N1", "F020_C8_T1", "F021_C8_N1")

macsPeaks <- list()
macsPeaks.gr <- GRangesList()
macsPeaks_distal.gr <- GRangesList()
for (i in 1:length(sample)){
  print(sample[i])
  macsPeaks[[i]] <- read.table(paste0(dir, sample[i], "_peaks.xls"), header = T)
  print(dim(macsPeaks[[i]]))
  macsPeaks[[i]]$chr <- paste0("chr", macsPeaks[[i]]$chr)
  macsPeaks[[i]] <- macsPeaks[[i]][macsPeaks[[i]]$chr %in% chrOrder,]
  print(dim(macsPeaks[[i]]))
  macsPeaks.gr[[i]] <- makeGRangesFromDataFrame(macsPeaks[i], keep.extra.columns = T)
  x <- as.data.frame(findOverlaps(macsPeaks.gr[[i]], hg19_prom))
  print(length(unique(x$queryHits)))
  macsPeaks_distal.gr[[i]] <- macsPeaks.gr[[i]][-unique(x$queryHits)]
}

setwd("/g/data/zk16/cc3704/replication_timing/human/thyroid")
F001_C1_T1.gr <- macsPeaks_distal.gr[[1]]
F002_C1_N1.gr <- macsPeaks_distal.gr[[2]]
F017_C7_T1.gr <- macsPeaks_distal.gr[[3]]
F018_C7_N1.gr <- macsPeaks_distal.gr[[4]]
F020_C8_T1.gr <- macsPeaks_distal.gr[[5]]
F021_C8_N1.gr <- macsPeaks_distal.gr[[6]]

#Defining peaks
C1_olaps <- findOverlaps(F001_C1_T1.gr, F002_C1_N1.gr)
C1_tumor <- subsetByOverlaps(F001_C1_T1.gr, F002_C1_N1.gr, invert = TRUE)
C1_normal <- subsetByOverlaps(F002_C1_N1.gr, F001_C1_T1.gr, invert = TRUE)
C1_common_1 <- F001_C1_T1.gr[queryHits(C1_olaps)]
C1_common_2 <- F002_C1_N1.gr[subjectHits(C1_olaps)]
C1_common <- GenomicRanges::union(C1_common_1, C1_common_2)

C7_olaps <- findOverlaps(F017_C7_T1.gr, F018_C7_N1.gr)
C7_tumor <- subsetByOverlaps(F017_C7_T1.gr, F018_C7_N1.gr, invert = TRUE)
C7_normal <- subsetByOverlaps(F018_C7_N1.gr, F017_C7_T1.gr, invert = TRUE)
C7_common_1 <- F017_C7_T1.gr[queryHits(C7_olaps)]
C7_common_2 <- F018_C7_N1.gr[subjectHits(C7_olaps)]
C7_common <- GenomicRanges::union(C7_common_1, C7_common_2)

C8_olaps <- findOverlaps(F020_C8_T1.gr, F021_C8_N1.gr)
C8_tumor <- subsetByOverlaps(F020_C8_T1.gr, F021_C8_N1.gr, invert = TRUE)
C8_normal <- subsetByOverlaps(F021_C8_N1.gr, F020_C8_T1.gr, invert = TRUE)
C8_common_1 <- F020_C8_T1.gr[queryHits(C8_olaps)]
C8_common_2 <- F021_C8_N1.gr[subjectHits(C8_olaps)]
C8_common <- GenomicRanges::union(C8_common_1, C8_common_2)

tumor_peaks <- GenomicRanges::union(C1_tumor, c(C7_tumor, C8_tumor))
normal_peaks <- GenomicRanges::union(C1_normal, c(C7_normal, C8_normal))
common_peaks <- GenomicRanges::union(C1_common, c(C7_common, C8_common))

tumor_peaks.df <- data.frame(chr=seqnames(tumor_peaks), start=start(tumor_peaks)-1
                             , end=end(tumor_peaks), stringsAsFactors = F)
normal_peaks.df <- data.frame(chr=seqnames(normal_peaks), start=start(normal_peaks)-1
                              , end=end(normal_peaks), stringsAsFactors = F)
common_peaks.df <- data.frame(chr=seqnames(common_peaks), start=start(common_peaks)-1
                              , end=end(common_peaks), stringsAsFactors = F)

write.table(x = tumor_peaks.df, file = "gains_thyroid_union_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
write.table(x = normal_peaks.df, file = "loses_thyroid_union_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
write.table(x = common_peaks.df, file = "unchanged_thyroid_union_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
