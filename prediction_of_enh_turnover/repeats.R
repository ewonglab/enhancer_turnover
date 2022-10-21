library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

# Function to read repeatmasker output
read_rm_tab <- function(rm_file){
  data <- read.table(rm_file, skip = 3, header = FALSE, sep = "\t")
  data$V1 <- gsub("\\s+", " ", str_trim(data$V1))
  rm_df <- as.data.frame(str_split_fixed(data$V1, " ", 15))  
  colnames(rm_df) <- c("SW_score", "perc_div.", "perc_del.", "perc_ins."
                       , "query_sequence", "query_start", "query_end"
                       , "query_left", "strand", "matching_repeat"
                       , "repeat_class/family", "repeat_start"
                       , "repeat_end", "repeat_left", "ID")
  return(rm_df)}

cal_mut_rate <- function(x){
  x$width <- with(x, V3-V2)
  x$mut_rate <- with(x, V4/width)
  return(x[,c("V1", "V2", "V3", "mut_rate")])
}

## Function to calculate the proportion of overlap with repetitive elements per genomic window
proportion_TE <- function(bins, repeats){#ignore strand of repeats
  #add ids
  colnames(bins) <- c("bin_seqnames", "bin_start", "bin_end")
  bins <- as.data.frame(bins)
  colnames(repeats) <- c("te_seqnames", "te_start", "te_end")
  repeats <- as.data.frame(repeats)
  repeats$te_start <- as.integer(repeats$te_start)
  repeats$te_end <- as.integer(repeats$te_end)
  
  bins.gr <- with(bins, GRanges(bin_seqnames , IRanges( bin_start+1, bin_end )))
  repeats.gr <- with(repeats, GRanges(te_seqnames , IRanges( te_start+1, te_end )))
  repeats.gr <- GenomicRanges::reduce(repeats.gr)
  repeats <- data.frame(te_seqnames=seqnames(repeats.gr), te_start=start(repeats.gr)-1
                      , te_end=end(repeats.gr), stringsAsFactors = F)
  
  x <- as.data.frame(findOverlaps(bins.gr, repeats.gr))
  bins$bin_id <- with(bins, paste(bin_seqnames , paste( bin_start, bin_end, sep = "_"), sep =":"))
  bins_ov_te <- cbind(bins[x$queryHits,], repeats[x$subjectHits,])
  # print(head(bins_ov_peak))
  bins_ov_te$adj_start <- ifelse(bins_ov_te$te_start < bins_ov_te$bin_start
                                   , bins_ov_te$bin_start, bins_ov_te$te_start)
  bins_ov_te$adj_end   <- ifelse(bins_ov_te$te_end > bins_ov_te$bin_end
                                   , bins_ov_te$bin_end, bins_ov_te$te_end)
  bins_ov_te$inter_width <- with(bins_ov_te, adj_end - adj_start)
  bins_ov_te <- as.data.frame(bins_ov_te)
  bins_ov_te$te_seqnames <- NULL
  bins_ov_te$te_start <- NULL
  bins_ov_te$te_end <- NULL
  bins_ov_te$adj_start  <- NULL
  bins_ov_te$adj_end <- NULL
  
  bins_ov_te <- aggregate(inter_width~., bins_ov_te, sum)
  bins_ov_te$bin_width <-with(bins_ov_te, bin_end-bin_start)
  bins_ov_te$prop <- with(bins_ov_te, inter_width/bin_width)
  return(bins_ov_te[,c("bin_id", "prop")])}


count_TE_per_bin <- function(bins, repeats){
  #add ids
  bins$bin_id <- with(bins, paste(V1 , paste( V2, V3, sep = "_"), sep =":"))
  repeats$query_start <- as.integer(repeats$query_start)
  repeats$query_end <- as.integer(repeats$query_end)
  repeats$repeat_id <- with(repeats, paste(query_sequence
                                           , paste( query_start, query_end, sep = "_")
                                           , sep =":"))
  #create gr objects
  bins.gr <- with(bins, GRanges(V1 , IRanges( V2+1, V3 )))
  repeats.gr <- with(repeats, GRanges(query_sequence, IRanges( query_start+1, query_end )))
  #cound number of peaks per genomic bin
  x <- as.data.frame(findOverlaps(bins.gr, repeats.gr))
  bins_ov_peaks <- cbind(bins[x$queryHits, "bin_id", drop =F]
                         , repeats[x$subjectHits, "repeat_id", drop =F])
  n_repeats <- as.data.frame(table(bins_ov_peaks$bin_id))
  return(n_repeats)}

setwd("/g/data/zk16/cc3704/replication_timing/human/enh_gene/k562")
args <- commandArgs(trailingOnly=TRUE)
genome_bins <- fread(args[1])

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]

print("Reading repeats")
out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
humanonly_f <- list.files(path = "/g/data/zk16/cc3704/replication_timing/human/repeats/humanonly_TE_hg19"
                          , pattern = "(.*)out", full.names = T)
for(i in 1:length(humanonly_f)){
  humanonly_te <- read_rm_tab(humanonly_f[i])
  bins_with_humanonlyTE_N <- count_TE_per_bin(genome_bins, humanonly_te)
  write.table(x = bins_with_humanonlyTE_N, file = paste(out_name, paste0("N_humanonly.TE_", basename(humanonly_f[i]), "txt"), sep = "_")
              , quote = F, row.names=F, col.names = F, sep = '\t')
  bins_with_humanonlyTE <- proportion_TE(genome_bins, humanonly_te[,c("query_sequence", "query_start", "query_end")])
  write.table(x = bins_with_humanonlyTE, file = paste(out_name, paste0("prop_humanonly.TE_", basename(humanonly_f[i]), "txt"), sep = "_"), quote = F, row.names=F, sep = '\t')
}

ancestral_f <- list.files(path = "/g/data/zk16/cc3704/replication_timing/human/repeats/ancestral_TE_hg19"
                          , pattern = "(.*)out", full.names = T)
for(i in 1:length(ancestral_f)){
  ancestral_te <- read_rm_tab(ancestral_f[i])
  bins_with_ancTE_N <- count_TE_per_bin(genome_bins, ancestral_te)
  write.table(x = bins_with_ancTE_N, file = paste(out_name, paste0("N_ancestral.TE_", basename(ancestral_f[i]), "txt"), sep = "_")
              , quote = F, row.names=F, col.names = F, sep = '\t')
  bins_with_ancTE <- proportion_TE(genome_bins, ancestral_te[,c("query_sequence", "query_start", "query_end")])
  write.table(x = bins_with_ancTE, file = paste(out_name, paste0("prop_ancestral.TE_", basename(ancestral_f[i]), "txt"), sep = "_"), quote = F, row.names=F, sep = '\t')
}

