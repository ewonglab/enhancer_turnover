library(data.table)
library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(stringr)
library("bedtoolsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

proportion_genes <- function(bins, genes){
  #add ids
  strand(genes) <- "*"
  genes <- GenomicRanges::reduce(genes)
  colnames(bins) <- c("seqnames", "start", "end")
  bins <- as.data.frame(bins)
  #create gr objects
  genes.df <- data.frame(seqnames=seqnames(genes), start=start(genes)-1
                         , end=end(genes), stringsAsFactors = F)#I will ignore strand
  bins_intersect <-  bedtoolsr::bt.intersect(bins, genes.df)
  colnames(bins_intersect) <- colnames(bins)
  #I don't need to sum 1 if coordinates are 0 based
  bins_intersect$width <- with(bins_intersect, end-start)
  bins_intersect.gr <- with(bins_intersect, GRanges(seqnames , IRanges( start+1, end )))
  bins_intersect.gr <- GenomicRanges::reduce(bins_intersect.gr)
  bins_intersect <- data.frame(seqnames=seqnames(bins_intersect.gr)
                               , start=start(bins_intersect.gr)-1
                         , end=end(bins_intersect.gr), stringsAsFactors = F)#I 
  bins.gr <- with(bins, GRanges(seqnames , IRanges( start+1, end )))
  bins$bin_id <- with(bins, paste(seqnames , paste( start, end, sep = "_"), sep =":"))
  #cound number of peaks per genomic bin
  x <- as.data.frame(findOverlaps(bins.gr, bins_intersect.gr))
  bins_ov_genes <- cbind(bins[x$queryHits, ]
                         , bins_intersect[x$subjectHits,])
  colnames(bins_ov_genes) <- c("bin_seqnames", "bin_start", "bin_end"
                  , "bin_id", "inter_seqnames", "inter_start", "inter_end")
  bins_ov_genes$adj_start <- ifelse(bins_ov_genes$inter_start < bins_ov_genes$bin_start
                                    , bins_ov_genes$bin_start, bins_ov_genes$inter_start)
  bins_ov_genes$adj_end   <- ifelse(bins_ov_genes$inter_end > bins_ov_genes$bin_end
                                    , bins_ov_genes$bin_end, bins_ov_genes$inter_end)
  bins_ov_genes$inter_width <- with(bins_ov_genes, adj_end - adj_start)
  bins_ov_genes <- as.data.frame(bins_ov_genes)
  bins_ov_genes$inter_seqnames <- NULL 
  bins_ov_genes$inter_start <- NULL
  bins_ov_genes$inter_end <- NULL
  bins_ov_genes$adj_start  <- NULL
  bins_ov_genes$adj_end <- NULL
  # print(head(bins_ov_genes))

  bins_ov_genes <- aggregate(inter_width~., bins_ov_genes, sum)
  bins_ov_genes$bin_width <-with(bins_ov_genes, bin_end-bin_start)
  bins_ov_genes$prop <- with(bins_ov_genes, inter_width/bin_width)
  return(bins_ov_genes[,c("bin_id", "prop")])
  }

proportion_genes2 <- function(bins, genes){
  #add ids
  strand(genes) <- "*"
  genes <- GenomicRanges::reduce(genes)
  colnames(bins) <- c("bin_seqnames", "bin_start", "bin_end")
  bins <- as.data.frame(bins)
  #create gr objects
  genes.df <- data.frame(gene_seqnames=seqnames(genes), gene_start=start(genes)-1
                         , gene_end=end(genes), stringsAsFactors = F)#I will ignore strand
  # bins_intersect <-  bedtoolsr::bt.intersect(bins, genes.df)
  bins.gr <- with(bins, GRanges(bin_seqnames , IRanges( bin_start+1, bin_end )))
  x <- as.data.frame(findOverlaps(bins.gr, genes))
  bins$bin_id <- with(bins, paste(bin_seqnames , paste( bin_start, bin_end, sep = "_"), sep =":"))
  bins_ov_genes <- cbind(bins[x$queryHits,], genes.df[x$subjectHits,])
  print(head(bins_ov_genes))
  bins_ov_genes$adj_start <- ifelse(bins_ov_genes$gene_start < bins_ov_genes$bin_start
                                    , bins_ov_genes$bin_start, bins_ov_genes$gene_start)
  bins_ov_genes$adj_end   <- ifelse(bins_ov_genes$gene_end > bins_ov_genes$bin_end
                                    , bins_ov_genes$bin_end, bins_ov_genes$gene_end)
  bins_ov_genes$inter_width <- with(bins_ov_genes, adj_end - adj_start)
  bins_ov_genes <- as.data.frame(bins_ov_genes)
  bins_ov_genes$gene_seqnames <- NULL
  bins_ov_genes$gene_start <- NULL
  bins_ov_genes$gene_end <- NULL
  bins_ov_genes$adj_start  <- NULL
  bins_ov_genes$adj_end <- NULL
  # print(head(bins_ov_genes))
  
  bins_ov_genes <- aggregate(inter_width~., bins_ov_genes, sum)
  bins_ov_genes$bin_width <-with(bins_ov_genes, bin_end-bin_start)
  bins_ov_genes$prop <- with(bins_ov_genes, inter_width/bin_width)
  return(bins_ov_genes[,c("bin_id", "prop")])
}

count_genes_per_bin <- function(bins, genes){
  #add ids
  bins$bin_id <- with(bins, paste(V1 , paste( V2, V3, sep = "_"), sep =":"))
  #create gr objects
  bins.gr <- with(bins, GRanges(V1 , IRanges( V2+1, V3 )))
  #cound number of peaks per genomic bin
  x <- as.data.frame(findOverlaps(bins.gr, genes))
  bins_ov_peaks <- cbind(bins[x$queryHits, "bin_id", drop =F]
                         , as.data.frame(genes)[x$subjectHits, "gene_id", drop =F])
  n_genes <- as.data.frame(table(bins_ov_peaks$bin_id))
  return(n_genes)}

setwd("/g/data/zk16/cc3704/replication_timing/human/enh_gene/k562")
args <- commandArgs(trailingOnly=TRUE)
genome_bins <- fread(args[1])

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX"
          , "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14"
          , "chr15", "chr16", "chr17", "chr18", "chr20", "chrY"
          , "chr19", "chr22", "chr21")
genome_bins <- genome_bins[genome_bins$V1 %in% chrs,]

print("reading genome annotation")
hg19_txdb <- makeTxDbFromGFF(file="/g/data/zk16/cc3704/human_data/hg19.ensGene.gtf"
                             , format="gtf",organism="Homo sapiens")
genes <- GenomicFeatures::genes(hg19_txdb)

# Count number of genes overlapping every genomic window
bins_with_genes_N <- count_genes_per_bin(genome_bins, genes)

# Calculate the proportion of overlap with genes per genomic window
bins_with_genes <- proportion_genes(genome_bins, genes)

colnames(bins_with_genes_N) <- c("bin_id", "N_genes")
colnames(bins_with_genes) <- c("bin_id", "prop_genes")

out_name <- gsub(pattern = "\\.bed", replacement = "", x = args[1])
write.table(x = bins_with_genes_N, file = paste(out_name, "N_genes.txt", sep = "_"), quote = F, row.names=F, sep = '\t')
write.table(x = bins_with_genes, file = paste(out_name, "prop_genes.txt", sep = "_"), quote = F, row.names=F, sep = '\t')
