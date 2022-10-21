library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

cal_mut_rate <- function(x){
  #I don't need to sum 1 if coordinates are 0-based
  x$width <- with(x, V3-V2)
  x$mut_rate <- with(x, V4/width)
  return(x[,c("V1", "V2", "V3", "mut_rate")])
}

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
#keep only main chromosomes
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]


#read mutations data and calculate local mutation rates
print("Reading mutations data")
anc_100kb <- fread("/g/data/zk16/cc3704/replication_timing/human/num_anc_variants_100kb_human_genome.bed")
rare_100kb <- fread("/g/data/zk16/cc3704/replication_timing/human/num_rare_variants_100kb_human_genome.bed")
anc_1mb <- fread("/g/data/zk16/cc3704/replication_timing/human/num_anc_variants_1mb_human_genome.bed")
rare_1mb <- fread("/g/data/zk16/cc3704/replication_timing/human/num_rare_variants_1mb_human_genome.bed")

anc_100kb <- cal_mut_rate(anc_100kb)
rare_100kb <- cal_mut_rate(rare_100kb)
anc_1mb <- cal_mut_rate(anc_1mb)
rare_1mb <- cal_mut_rate(rare_1mb)

bins_with_anc_100kb <- continuous_var_per_bin(genome_bins, anc_100kb)
bins_with_rare_100kb <- continuous_var_per_bin(genome_bins, rare_100kb)
bins_with_anc_1mb <- continuous_var_per_bin(genome_bins, anc_1mb)
bins_with_rare_1mb <- continuous_var_per_bin(genome_bins, rare_1mb)


colnames(bins_with_anc_100kb) <- c("bin_id", "anc_100kb_mutRate")
colnames(bins_with_anc_1mb) <- c("bin_id", "anc_1mb_mutRate")
colnames(bins_with_rare_100kb) <- c("bin_id", "rare_100kb_mutRate")
colnames(bins_with_rare_1mb) <- c("bin_id", "rare_1mb_mutRate")

out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
write.table(x = bins_with_anc_100kb, file = paste(out_name, "anc_100kb_mutRate.txt", sep = "_")
, quote = F, row.names=F, sep = '\t')

write.table(x = bins_with_anc_1mb, file = paste(out_name, "anc_1mb_mutRate.txt", sep = "_")
, quote = F, row.names=F, sep = '\t')

write.table(x = bins_with_rare_100kb, file = paste(out_name, "rare_100kb_mutRate.txt", sep = "_")
, quote = F, row.names=F, sep = '\t')

write.table(x = bins_with_rare_1mb, file = paste(out_name, "rare_1mb_mutRate.txt", sep = "_")
, quote = F, row.names=F, sep = '\t')
