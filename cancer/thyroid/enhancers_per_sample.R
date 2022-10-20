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

#Defining enhancers
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

#C1 enhancers
C1_tumor.df <- data.frame(chr=seqnames(C1_tumor), start=start(C1_tumor)-1
                             , end=end(C1_tumor), stringsAsFactors = F)
C1_normal.df <- data.frame(chr=seqnames(C1_normal), start=start(C1_normal)-1
                              , end=end(C1_normal), stringsAsFactors = F)
C1_common.df <- data.frame(chr=seqnames(C1_common), start=start(C1_common)-1
                              , end=end(C1_common), stringsAsFactors = F)

#C7 enhancers
C7_tumor.df <- data.frame(chr=seqnames(C7_tumor), start=start(C7_tumor)-1
                             , end=end(C7_tumor), stringsAsFactors = F)
C7_normal.df <- data.frame(chr=seqnames(C7_normal), start=start(C7_normal)-1
                              , end=end(C7_normal), stringsAsFactors = F)
C7_common.df <- data.frame(chr=seqnames(C7_common), start=start(C7_common)-1
                              , end=end(C7_common), stringsAsFactors = F)

#C8 enhancers
C8_tumor.df <- data.frame(chr=seqnames(C8_tumor), start=start(C8_tumor)-1
                             , end=end(C8_tumor), stringsAsFactors = F)
C8_normal.df <- data.frame(chr=seqnames(C8_normal), start=start(C8_normal)-1
                              , end=end(C8_normal), stringsAsFactors = F)
C8_common.df <- data.frame(chr=seqnames(C8_common), start=start(C8_common)-1
                              , end=end(C8_common), stringsAsFactors = F)

#add enhancer type
C1_tumor.df$type <- "gain"
C1_normal.df$type <- "lose"
C1_common.df$type <- "unchanged"

C7_tumor.df$type <- "gain"
C7_normal.df$type <- "lose"
C7_common.df$type <- "unchanged"

C8_tumor.df$type <- "gain"
C8_normal.df$type <- "lose"
C8_common.df$type <- "unchanged"

c1_all <- rbind(C1_tumor.df, C1_normal.df, C1_common.df)
c7_all <- rbind(C7_tumor.df, C7_normal.df, C7_common.df)
c8_all <- rbind(C8_tumor.df, C8_normal.df, C8_common.df)


write.table(x = c1_all, file = "thyroid_C1_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
write.table(x = c7_all, file = "thyroid_C7_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
write.table(x = c8_all, file = "thyroid_C8_1kb.bed"
            , quote = F, col.names = F, row.names = F, sep = '\t')
