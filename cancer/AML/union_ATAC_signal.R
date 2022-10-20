#Prediction of pre-leukemic HSC from narrowPeak ATAC-seq files

library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library(data.table)
library(reshape2)


#Loading macs3 output
chrOrder<-paste0("chr", c(1:22,"X","Y"))


dir = "/g/data/zk16/veronika/projects/repli_time/AML/ATAC/peaks/"
sample <- c("SU484_blast", "SU484_pHSC", "SU501_blast"
            , "SU501_pHSC", "SU654_blast", "SU654_pHSC")
macsPeaks <- list()
macsPeaks.gr <- GRangesList()
for (i in 1:length(sample)){
  print(sample[i])
  #narrowPeaks are 0-based coordinates
  macsPeaks[[i]] <- read.table(paste0(dir, sample[i], "_peaks.narrowPeak"), header = F)
  # print(dim(macsPeaks[[i]]))
  macsPeaks[[i]]$V1 <- paste0("chr", macsPeaks[[i]]$V1)
  macsPeaks[[i]] <- macsPeaks[[i]][macsPeaks[[i]]$V1 %in% chrOrder,]
  # print(dim(macsPeaks[[i]]))
  # print(class(macsPeaks[[i]]))
  # print(head(macsPeaks[[i]]))
  macsPeaks.gr[[i]] <- with(macsPeaks[[i]], GRanges(V1 , IRanges( V2+1, V3 )))
}

setwd("/g/data/zk16/cc3704/replication_timing/human/AML")
SU484_blast.gr <- macsPeaks.gr[[1]]
SU484_pHSC.gr <- macsPeaks.gr[[2]]
SU501_blast.gr <- macsPeaks.gr[[3]]
SU501_pHSC.gr <- macsPeaks.gr[[4]]
SU654_blast.gr <- macsPeaks.gr[[5]]
SU654_pHSC.gr <- macsPeaks.gr[[6]]

SU484_normal <- subsetByOverlaps(SU484_pHSC.gr, SU484_blast.gr, invert = TRUE)
SU501_normal <- subsetByOverlaps(SU501_pHSC.gr, SU501_blast.gr, invert = TRUE)
SU654_normal <- subsetByOverlaps(SU654_pHSC.gr, SU654_blast.gr, invert = TRUE)

SU484_olaps <- findOverlaps(SU484_blast.gr, SU484_pHSC.gr)
SU484_common_1 <- SU484_blast.gr[queryHits(SU484_olaps)]
SU484_common_2 <- SU484_pHSC.gr[subjectHits(SU484_olaps)]
SU484_common <- GenomicRanges::union(SU484_common_1, SU484_common_2)

SU501_olaps <- findOverlaps(SU501_blast.gr, SU501_pHSC.gr)
SU501_common_1 <- SU501_blast.gr[queryHits(SU501_olaps)]
SU501_common_2 <- SU501_pHSC.gr[subjectHits(SU501_olaps)]
SU501_common <- GenomicRanges::union(SU501_common_1, SU501_common_2)

SU654_olaps <- findOverlaps(SU654_blast.gr, SU654_pHSC.gr)
SU654_common_1 <- SU654_blast.gr[queryHits(SU654_olaps)]
SU654_common_2 <- SU654_pHSC.gr[subjectHits(SU654_olaps)]
SU654_common <- GenomicRanges::union(SU654_common_1, SU654_common_2)

common_peaks <- GenomicRanges::union(SU484_common, c(SU501_common, SU654_common))
normal_peaks <- GenomicRanges::union(SU484_normal, c(SU501_normal, SU654_normal))

#Overlap union of enhancers with MACS3 output and calculate average signal for every region

names(macsPeaks) <- sample
macsPeaks.df <- do.call("rbind", macsPeaks)
dim(macsPeaks.df)
dim(unique(macsPeaks.df[,c("V1", "V2", "V3")]))#465372      3

#macs output is  0-based (narrowPeak)
macsPeaks.df.gr <- with(macsPeaks.df, GRanges(V1 , IRanges( V2+1, V3 )))
macsPeaks.df$sample <- rownames(macsPeaks.df)

x <- as.data.frame(findOverlaps(common_peaks, macsPeaks.df.gr))
common_signal <- cbind(as.data.frame(common_peaks)[x$queryHits,]
                       , macsPeaks.df[x$subjectHits, c("V7", "sample")])

#start position back to 0-based
common_signal$start <- common_signal$start-1
common_signal$sample <- gsub("\\..*", "", common_signal$sample)

# Keep only losses and unchanged enhancers  

common_signal <- common_signal[common_signal$sample %in%
                                 c("SU484_pHSC", "SU501_pHSC", "SU654_pHSC" ),]

# calculate mean signal
common_signal$sample <- NULL
common_signal <- aggregate(V7 ~ ., common_signal[,c("seqnames", "start", "end", "width"
                                                   , "V7")], mean)#66575     5

x <- as.data.frame(findOverlaps(normal_peaks, macsPeaks.df.gr))
normal_signal <- cbind(as.data.frame(normal_peaks)[x$queryHits,]
                       , macsPeaks.df[x$subjectHits, c("V7", "sample")])

normal_signal$sample <- gsub("\\..*", "", normal_signal$sample)
unique(normal_signal$sample)
normal_signal <- normal_signal[normal_signal$sample %in%
                                 c("SU484_pHSC", "SU501_pHSC", "SU654_pHSC" ),]
normal_signal$sample <- NULL
normal_signal <- aggregate(V7 ~ ., normal_signal[,c("seqnames", "start", "end", "width"
                                                   , "V7")], mean)#57491     5

peaks <- rbind(common_signal, normal_signal)
peaks$width <- NULL
peaks <- aggregate(V7 ~., peaks, mean)

write.table(x = peaks, file = "AML_pHSC_mean_signal.bed", quote = F, col.names = F
            , row.names = F, sep = '\t')

#scale signal values
peaks$V7 <- scale(x = peaks$V7, center = T, scale = T)
peaks$start <- format(peaks$start, scientific = F)
peaks$end <- format(peaks$end, scientific = F)

write.table(x = peaks, file = "AML_pHSC_mean_signal.scaled.bed", quote = F, col.names = F
            , row.names = F, sep = '\t')
