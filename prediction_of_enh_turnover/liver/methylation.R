library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

count_methCs <- function(bins, peaks){

  #peak id as it is in the bed file V2 and V3 meaning the position of C
  peaks$peak_id <- with(peaks, paste(V1 , paste( V2, V3, sep = "_"), sep =":"))
  bins$bin_id <- with(bins, paste(V1 , paste( V2, V3, sep = "_"), sep =":"))
  #create gr objects
  #including CpG positions
  peaks.gr <- with(peaks, GRanges(V1 , IRanges( V3, V3 +1)))
  bins.gr <- with(bins, GRanges(V1 , IRanges( V2+1, V3 )))
  x <- as.data.frame(findOverlaps(bins.gr, peaks.gr))
  bins_ov_peaks <- cbind(bins[x$queryHits, "bin_id", drop =F], peaks[x$subjectHits, "peak_id", drop =F])
  n_peaks <- as.data.frame(table(bins_ov_peaks$bin_id))
  return(n_peaks)}


mean_meth_per_bin <- function(bins, peaks){
  colnames(peaks) <- c("V1", "V2", "V3", "my_var")
  peaks$peak_id <- with(peaks, paste(V1 , paste( V2, V3, sep = "_"), sep = ":"))
  bins$bin_id <- with(bins, paste(V1 , paste( V2, V3, sep = "_"), sep = ":"))
  #create gr objects
  peaks.gr <- with(peaks, GRanges(V1 , IRanges( V3, V3+1 )))
  bins.gr <- with(bins, GRanges(V1 , IRanges( V2+1, V3 )))
  x <- as.data.frame(findOverlaps(bins.gr, peaks.gr))
  bins_ov_peaks <- cbind(bins[x$queryHits, ], peaks[x$subjectHits, "my_var", drop =F])
  #aggregate continuous variable to calculate a single value per genomic bin
  peaks_with_var <- aggregate(my_var~bin_id, bins_ov_peaks, mean)
  return(peaks_with_var)}


setwd("/g/data/zk16/cc3704/replication_timing/human/enh_gene/liver")
args <- commandArgs(trailingOnly=TRUE)
genome_bins <- fread(paste0("../k562/", args[1]))

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]

#METHYLATION
WGBS <- fread("ENCFF080ETP_WGBS_CpG_hg19.bed")
WGBS$coverage <- gsub(pattern = "_.*", replacement = "", WGBS$V7)
WGBS$perc_meth_reads <- gsub(pattern = ".*_", replacement = "", WGBS$V7)
WGBS$coverage <- as.numeric(WGBS$coverage)
WGBS$perc_meth_reads <- as.numeric(WGBS$perc_meth_reads)
WGBS <- WGBS[WGBS$coverage >= 10,]
WGBS <- WGBS[WGBS$perc_meth_reads >= 15,]

# keep only sites that map correctly to CpG sites in Hg19
WGBS_CpG.gr <- with(WGBS, GRanges(V1 , IRanges( V3, V3+1 ), strand = V6))

#getting sequences
WGBS_CpG_seq <- (Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, WGBS_CpG.gr))
WGBS_CpG_seq.df <- as.data.frame(WGBS_CpG_seq)
tmp <- as.data.frame(table(WGBS_CpG_seq.df$x))

#remove positions that do not map to CpGs in hg19
if(length(setdiff(tmp$Var1, "C"))>0){
  WGBS <- WGBS[-which(WGBS_CpG_seq.df$x !="CG"),]
}

bins_with_metCpG_N <- count_methCs(genome_bins, unique(WGBS[,1:3]))
bins_with_mean_methCpG <- mean_meth_per_bin(genome_bins, unique(WGBS[,c(1:3,9)]))
colnames(bins_with_metCpG_N) <- c("bin_id", "N_metCpG")

colnames(bins_with_mean_methCpG) <- c("bin_id", "mean_methylation_at_MethCpG")

out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
write.table(x = bins_with_metCpG_N, file = paste(out_name, "N_metCpG.txt", sep = "_"), quote = F
, row.names=F, sep = '\t')
write.table(x = bins_with_mean_methCpG, file = paste(out_name, "mean_meth_at_MethCpG.txt", sep = "_")
, quote = F, row.names=F, sep = '\t')
