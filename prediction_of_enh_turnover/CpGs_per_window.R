library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

# Function to count CpGs per genomic window
count_CpGs_per_bin <- function(bins, CpG_sites=cpgr){
  #add ids
  bins$bin_id <- with(bins, paste(V1 , paste( V2, V3, sep = "_"), sep =":"))
  #create gr objects
  bins.gr <- with(bins, GRanges(V1 , IRanges( V2+1, V3 )))
  #cound number of peaks per genomic bin
  x <- as.data.frame(findOverlaps(bins.gr, CpG_sites))
  bins_ov_cpgs <- cbind(bins[x$queryHits, "bin_id", drop =F]
                        , as.data.frame(CpG_sites)[x$subjectHits, ])
  n_CpGs <- as.data.frame(table(bins_ov_cpgs$bin_id))
  return(n_CpGs)}

setwd("/g/data/zk16/cc3704/replication_timing/human/enh_gene/k562")
args <- commandArgs(trailingOnly=TRUE)
genome_bins <- fread(args[1])

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]

# Getting hg19 CpGs

hg19 <- BSgenome.Hsapiens.UCSC.hg19
chrs <- names(Hsapiens)[1:24]

#getting start positions of CGs
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))


#Number of CpGs per bin

bins_with_CpG_N <- count_CpGs_per_bin(genome_bins, cpgr)
colnames(bins_with_CpG_N) <- c("bin_id", "N_CpG")

out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
write.table(x = bins_with_CpG_N, file = paste(out_name, "N_CpG.txt", sep = "_"), quote = F, row.names=F, sep = '\t')
