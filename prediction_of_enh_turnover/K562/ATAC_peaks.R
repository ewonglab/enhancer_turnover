library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

# Function to count the number of ATAC-seq peaks overlapping every genomic window
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

# Function to calculate the proportion of overlap with ATAC-seq peaks per genomic window

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
  # print(head(bins_ov_peak))
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


setwd("/g/data/zk16/cc3704/replication_timing/human/enh_gene/k562")
args <- commandArgs(trailingOnly=TRUE)
genome_bins <- fread(args[1])

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]

# reading ATAC-seq peaks   
atac <- fread("ENCFF333TAT_ATAC_hg19.bed")

atac <- unique(atac[,1:3])
bins_with_atac_N <- count_peaks_per_bin(genome_bins, atac)
bins_with_atac <- proportion_peaks(genome_bins, atac)
colnames(bins_with_atac_N) <- c("bin_id", "N_ATAC")
colnames(bins_with_atac) <- c("bin_id", "prop_ATAC")

out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
write.table(x = bins_with_atac_N, file = paste(out_name, "N_ATAC.txt", sep = "_"), quote = F, row.names=F, sep = '\t')
write.table(x = bins_with_atac, file = paste(out_name, "prop_ATAC.txt", sep = "_"), quote = F, row.names=F, sep = '\t')
