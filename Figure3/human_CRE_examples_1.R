# This script was used to get an initial subset of candidate human-specific enhancers and non-functional regions in the mouse genome for the sequence alignment examples in Figure 3.
# The candidate human enhancers were filtered to keep those that overlap either Hnf4a or CEBPA binding sites (ChIP-seq) 
# The sequences of the selected human candidate enhancers and non-functional regions in the mouse genome were saved for further analysis.

library(data.table)
library("GenomeInfoDb")
library("GenomicFeatures")
library("AnnotationFilter")
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library(Biostrings)
library("seqinr")

# reading human candidate enhancers and aligned non-functional elements in mouse genome
human_mm10 <- read.table(file = "Hsap_Enhancers_not_Mmus_alignMmus_mm10.bed"
                         , header = F, stringsAsFactors = F, sep ='\t')# non-functional in mouse genome
human_hg38 <- read.table(file = "Hsap_Enhancers_not_Mmus_alignMmus_hg38.bed"
                         , header = F, stringsAsFactors = F, sep ='\t')# human candidate enhancers

colnames(human_mm10)[1:3] <- c("chr_mm10", "start_mm10", "end_mm10")
colnames(human_hg38)[1:3] <- c("chr_hg38", "start_hg38", "end_hg38")

# merge pairs of candidate enhancers and non-functional aligned region
human_enh <- merge(human_hg38, human_mm10, by="V4")# 

mouse_non_fun_gr <- with(unique(human_enh[,c("chr_mm10", "start_mm10", "end_mm10")])
                         , GRanges(chr_mm10, IRanges(start_mm10+1, end_mm10)))#
human_cre_gr <- with(unique(human_enh[,c("chr_hg38", "start_hg38", "end_hg38")])
                     , GRanges(chr_hg38, IRanges(start_hg38+1, end_hg38)))#

# genomes
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
mm10 <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10

# reading chipseq datasets
chipseq_path <- "/g/data/zk16/cc3704/replication_timing/da_model/raw_data/"
Hnf4a_human <- read.table(file = paste0(chipseq_path, "hg38/Hnf4a/Hnf4a_hg38.bed")
                          , header = F, stringsAsFactors = F, sep ='\t')
CEBPA_human <- read.table(file = paste0(chipseq_path, "hg38/CEBPA/CEBPA_hg38.bed")
                          , header = F, stringsAsFactors = F, sep ='\t')

# GenomicRanges objects
CEBPA_human_gr <- with(CEBPA_human, GRanges(V1, IRanges(V2+1, V3)))
Hnf4a_human_gr <- with(Hnf4a_human, GRanges(V1, IRanges(V2+1, V3)))
human_enh_gr <- with(human_enh, GRanges(chr_hg38, IRanges(start_hg38+1, end_hg38)))

## overlapping human enhancers with Chipseq datasets
x <- as.data.frame(findOverlaps(human_enh_gr, CEBPA_human_gr))
human_ov_CEBPA <- cbind(human_enh[x$queryHits,], CEBPA_human[x$subjectHits,])#
x <- as.data.frame(findOverlaps(human_enh_gr, Hnf4a_human_gr))
human_ov_Hnf4a <- cbind(human_enh[x$queryHits,], Hnf4a_human[x$subjectHits,])#

# saving sequences of human candidate enhancers that overlap Chip-seq TF binding sites and non-functional regions 
tmp <- getSeq(hg38, with(unique(human_ov_CEBPA[,c("chr_hg38", "start_hg38", "end_hg38")])
                         , GRanges(chr_hg38, IRanges(start_hg38+1, end_hg38))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(human_ov_CEBPA[,c("chr_hg38", "start_hg38", "end_hg38")])
                           , paste(chr_hg38, paste(start_hg38, end_hg38, sep ='-'), sep = ":"))
            , file.out = "hg38_enh_ov_CEBPA.fa")

tmp <- getSeq(mm10, with(unique(human_ov_CEBPA[,c("chr_mm10", "start_mm10", "end_mm10")])
                         , GRanges(chr_mm10, IRanges(start_mm10+1, end_mm10))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(human_ov_CEBPA[,c("chr_mm10", "start_mm10", "end_mm10")])
                           , paste(chr_mm10, paste(start_mm10, end_mm10, sep ='-'), sep = ":"))
            , file.out = "mm10_non_functional_ov_hg38_CEBPA.fa")

tmp <- getSeq(hg38, with(unique(human_ov_Hnf4a[,c("chr_hg38", "start_hg38", "end_hg38")])
                         , GRanges(chr_hg38, IRanges(start_hg38+1, end_hg38))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(human_ov_Hnf4a[,c("chr_hg38", "start_hg38", "end_hg38")])
                           , paste(chr_hg38, paste(start_hg38, end_hg38, sep ='-'), sep = ":"))
            , file.out = "hg38_enh_ov_Hnf4a.fa")


tmp <- getSeq(mm10, with(unique(human_ov_Hnf4a[,c("chr_mm10", "start_mm10", "end_mm10")])
                         , GRanges(chr_mm10, IRanges(start_mm10+1, end_mm10))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(human_ov_Hnf4a[,c("chr_mm10", "start_mm10", "end_mm10")])
                           , paste(chr_mm10, paste(start_mm10, end_mm10, sep ='-'), sep = ":"))
            , file.out = "mm10_non_functional_ov_hg38_Hnf4a.fa")
