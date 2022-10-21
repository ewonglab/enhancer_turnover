library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

count_peaks_per_bin <- function(bins, peaks){
  #add ids
  peaks$peak_id <- with(peaks, paste(V1 , paste( V2, V3, sep = "_"), sep =":"))
  bins$bin_id <- with(bins, paste(V1 , paste( V2, V3, sep = "_"), sep =":"))
  #create gr objects
  peaks.gr <- with(peaks, GRanges(V1 , IRanges( V2+1, V3 )))
  bins.gr <- with(bins, GRanges(V1 , IRanges( V2+1, V3 )))
  #cound number of peaks per genomic bin
  x <- as.data.frame(findOverlaps(bins.gr, peaks.gr))
  bins_ov_peaks <- cbind(bins[x$queryHits, "bin_id", drop =F], peaks[x$subjectHits, "peak_id", drop =F])
  n_peaks <- as.data.frame(table(bins_ov_peaks$bin_id))
  return(n_peaks)}

proportion_peaks <- function(bins, peaks){
  #add ids
  colnames(bins) <- c("bin_seqnames", "bin_start", "bin_end")
  colnames(peaks) <- c("peak_seqnames", "peak_start", "peak_end")
  bins <- as.data.frame(bins)
  peaks <- as.data.frame(peaks)
  bins.gr <- with(bins, GRanges(bin_seqnames , IRanges( bin_start+1, bin_end )))
  peaks.gr <- with(peaks, GRanges(peak_seqnames , IRanges( peak_start+1, peak_end )))
  peaks.gr <- GenomicRanges::reduce(peaks.gr)
  peaks <- data.frame(peak_seqnames=seqnames(peaks.gr), peak_start=start(peaks.gr)-1
                         , peak_end=end(peaks.gr), stringsAsFactors = F)
  
  x <- as.data.frame(findOverlaps(bins.gr, peaks.gr))
  bins$bin_id <- with(bins, paste(bin_seqnames , paste( bin_start, bin_end, sep = "_"), sep =":"))
  bins_ov_peak <- cbind(bins[x$queryHits,], peaks[x$subjectHits,])
  bins_ov_peak$adj_start <- ifelse(bins_ov_peak$peak_start < bins_ov_peak$bin_start
                                    , bins_ov_peak$bin_start, bins_ov_peak$peak_start)
  bins_ov_peak$adj_end   <- ifelse(bins_ov_peak$peak_end > bins_ov_peak$bin_end
                                    , bins_ov_peak$bin_end, bins_ov_peak$peak_end)
  bins_ov_peak$inter_width <- with(bins_ov_peak, adj_end - adj_start)
  bins_ov_peak <- as.data.frame(bins_ov_peak)
  bins_ov_peak$peak_seqnames <- NULL
  bins_ov_peak$peak_start <- NULL
  bins_ov_peak$peak_end <- NULL
  bins_ov_peak$adj_start  <- NULL
  bins_ov_peak$adj_end <- NULL
  
  bins_ov_peak <- aggregate(inter_width~., bins_ov_peak, sum)
  bins_ov_peak$bin_width <-with(bins_ov_peak, bin_end-bin_start)
  bins_ov_peak$prop <- with(bins_ov_peak, inter_width/bin_width)
  return(bins_ov_peak[,c("bin_id", "prop")])}

setwd("/g/data/zk16/cc3704/replication_timing/human/enh_gene/liver")
args <- commandArgs(trailingOnly=TRUE)
genome_bins <- fread(paste0("../k562/", args[1]))

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]


#HISTONE MARKS
H3K27ac <- fread("ENCFF837ELM_H3K27ac_hg19.bed")
H3K4me2 <- fread("ENCFF557QJF_H3K4me2_hg19.bed.gz")
H4K20me1 <- fread("ENCFF244UXH_H4K20me1_hg19.bed.gz")
H3K9ac <- fread("ENCFF128PTM_H3K9ac_hg19.bed.gz")
H3K4me1 <- fread("ENCFF608RNK_H3K4me1_hg19.bed")
H3K27me3 <- fread("ENCFF877YOY_H3K27me3_hg19.bed")
H3K36me3 <- fread("ENCFF734MXG_H3K36me3_hg19.bed")
H2AFZ <- fread("ENCFF987QKQ_H2AFZ_hg19.bed.gz")
H3K4me3 <- fread("ENCFF830CCA_H3K4me3_hg19.bed")
H3K9me2 <- fread("ENCFF273GGC_H3K9me2_hg19.bed.gz")
H3K9me3 <- fread("ENCFF681HRQ_H3K9me3_hg19.bed")
H3K79me2 <- fread("ENCFF999FNX_H3K79me2_hg19.bed.gz")

# make a list with histone marks datasets
histones <- list(unique(H3K27ac[,1:3]), unique(H3K4me2[,1:3]), unique(H4K20me1[,1:3])
                 , unique(H3K9ac[,1:3]), unique(H3K4me1[,1:3]), unique(H3K27me3[,1:3])
                 , unique(H3K36me3[,1:3]), unique(H2AFZ[,1:3]), unique(H3K4me3[,1:3])
                 , unique(H3K9me2[,1:3]), unique(H3K9me3[,1:3]), unique(H3K79me2[,1:3]))

names(histones) <- c("H3K27ac", "H3K4me2", "H4K20me1"
                     , "H3K9ac", "H3K4me1", "H3K27me3"
                     , "H3K36me3", "H2AFZ", "H3K4me3"
                     , "H3K9me2", "H3K9me3", "H3K79me2")


bins_with_histones_N <- lapply(histones, function(x) count_peaks_per_bin(genome_bins, x))
names(bins_with_histones_N) <- c("H3K27ac", "H3K4me2", "H4K20me1"
                                 , "H3K9ac", "H3K4me1", "H3K27me3"
                                 , "H3K36me3", "H2AFZ", "H3K4me3"
                                 , "H3K9me2", "H3K9me3", "H3K79me2")


bins_with_histones <- lapply(histones, function(x) proportion_peaks(genome_bins, x))
names(bins_with_histones) <- c("H3K27ac", "H3K4me2", "H4K20me1"
                               , "H3K9ac", "H3K4me1", "H3K27me3"
                               , "H3K36me3", "H2AFZ", "H3K4me3"
                               , "H3K9me2", "H3K9me3", "H3K79me2")

for(i in 1:length(bins_with_histones_N)){
  colnames(bins_with_histones_N[[i]]) <- c("bin_id", paste0("N_", names(bins_with_histones_N[i])))}


for(i in 1:length(bins_with_histones)){
  colnames(bins_with_histones[[i]]) <- c("bin_id", paste0("prop_", names(bins_with_histones[i])))}

out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
for(i in 1:length(bins_with_histones_N)){
  write.table(x = bins_with_histones_N[[i]]
  , file = paste(out_name, paste0(names(bins_with_histones_N[i]), "_N.txt"), sep = "_"), quote = F
  , row.names=F, sep = '\t')}

for(i in 1:length(bins_with_histones)){
  write.table(x = bins_with_histones[[i]], file = paste(out_name, paste0(names(bins_with_histones[i]), ".txt")
  , sep = "_"), quote = F, row.names=F, sep = '\t')}
