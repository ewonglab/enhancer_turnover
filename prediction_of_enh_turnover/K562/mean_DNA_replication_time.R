library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

continuous_var_per_bin <- function(bins, peaks){
  #make sure variable is in column number 4
  colnames(peaks) <- c("V1", "V2", "V3", "my_var")
  #add ids
  peaks$peak_id <- with(peaks, paste(V1 , paste( V2, V3, sep = "_"), sep = ":"))
  bins$bin_id <- with(bins, paste(V1 , paste( V2, V3, sep = "_"), sep = ":"))
  #create gr objects
  peaks.gr <- with(peaks, GRanges(V1 , IRanges( V2+1, V3 )))
  bins.gr <- with(bins, GRanges(V1 , IRanges( V2+1, V3 )))
  #cound number of peaks per genomic bin
  x <- as.data.frame(findOverlaps(bins.gr, peaks.gr))
  bins_ov_peaks <- cbind(bins[x$queryHits, ], peaks[x$subjectHits, "my_var", drop =F])
  #aggregate continuous variable to calculate a single value per genomic bin
  peaks_with_var <- aggregate(my_var~bin_id, bins_ov_peaks, mean)
  return(peaks_with_var)}

setwd("/g/data/zk16/cc3704/replication_timing/human/enh_gene/k562")
args <- commandArgs(trailingOnly=TRUE)
genome_bins <- fread(args[1])

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]

#REPLICATION TIME DATA SETS

#K562 DNA replication time
k562_rt <- fread("GSM923448_hg19_wgEncodeUwRepliSeqK562WaveSignalRep1.bedGraph")

#H9 ESC DNA replication time
h9_rt <- fread("/g/data/zk16/cc3704/replication_timing/human/RT_H9_ESC_Ext29405702_hg19.bedgraph")

# replication time predicted from ovary ATAC-seq
replicon_ovary_files <- list.files(path = "/g/data/zk16/cc3704/replication_timing/human/ovary_hg19_IPLS"
                                   , pattern = "(.*)timing$", full.names = T)
                                   
# replication time predicted from testis ATAC-seq
replicon_testis_files <- list.files(path = "/g/data/zk16/cc3704/replication_timing/human/testis_hg19_IPLS"
                                    , pattern = "(.*)timing$", full.names = T)
                                    
replicon_ovary <- lapply(replicon_ovary_files, fread)
replicon_testis <- lapply(replicon_testis_files, fread)
names(replicon_ovary) <- basename(replicon_ovary_files)
names(replicon_testis) <- basename(replicon_testis_files)

for(i in 1:length(replicon_ovary)){
  replicon_ovary[[i]]$chr <- names(replicon_ovary[i])}

for(i in 1:length(replicon_testis)){
  replicon_testis[[i]]$chr <- names(replicon_testis[i])}

replicon_ovary <- do.call("rbind", replicon_ovary)
replicon_testis <- do.call("rbind", replicon_testis)

# add coordinates to predicted ovary and testis DNA replication time
replicon_ovary$chr <- gsub("\\.timing", "", replicon_ovary$chr)
replicon_ovary$end <- replicon_ovary$V1 + 500
replicon_ovary <- replicon_ovary[,c("chr", "V1", "end", "V2")]

replicon_testis$chr <- gsub("\\.timing", "", replicon_testis$chr)
replicon_testis$end <- replicon_testis$V1 + 500
replicon_testis <- replicon_testis[,c("chr", "V1", "end", "V2")]

# correct direction of predicted DNA replication time
replicon_ovary$V2 <- (replicon_ovary$V2) * (-1)
replicon_testis$V2 <- (replicon_testis$V2) * (-1)

bins_with_k562_rt <- continuous_var_per_bin(genome_bins, k562_rt)
bins_with_h9_rt <- continuous_var_per_bin(genome_bins, h9_rt)
bins_with_ovary_rt <- continuous_var_per_bin(genome_bins, replicon_ovary)
bins_with_testis_rt <- continuous_var_per_bin(genome_bins, replicon_testis)

colnames(bins_with_k562_rt) <- c("bin_id", "k562_rt")
colnames(bins_with_ovary_rt) <- c("bin_id", "ovary_rt")
colnames(bins_with_testis_rt) <- c("bin_id", "testis_rt")
colnames(bins_with_h9_rt) <- c("bin_id", "H9_rt")

out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
write.table(x = bins_with_k562_rt, file = paste(out_name, "k562_rt.txt", sep = "_"), quote = F
, row.names=F, sep = '\t')
write.table(x = bins_with_ovary_rt, file = paste(out_name, "ovary_rt.txt", sep = "_"), quote = F
, row.names=F, sep = '\t')
write.table(x = bins_with_testis_rt, file = paste(out_name, "testis_rt.txt", sep = "_"), quote = F
, row.names=F, sep = '\t')
write.table(x = bins_with_h9_rt, file = paste(out_name, "H9_rt.txt", sep = "_"), quote = F, row.names=F
, sep = '\t')
