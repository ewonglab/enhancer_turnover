# Script used to subset mouse candidate enhancers to those that overlap CEBPA or HNF4 binding sites (ChIP-seq)
# The sequences of the filtered candidate enhancers and aligned non-functional regions in the human region were saved for further analysis

library(data.table)
library("GenomeInfoDb")
library("GenomicFeatures")
library("AnnotationFilter")
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library(Biostrings))
library("seqinr")

# Reading enhancers and aligned non-functional elements
mouse_mm10 <- read.table(file = "Mmus_Enhancers_not_Hsap_alignHsap_mm10.bed"
                         , header = F, stringsAsFactors = F, sep ='\t')# mouse candidate enhancers
mouse_hg38 <- read.table(file = "Mmus_Enhancers_not_Hsap_alignHsap_hg38.bed"
                         , header = F, stringsAsFactors = F, sep ='\t')# aligned non-functional hg38

colnames(mouse_mm10)[1:3] <- c("chr_mm10", "start_mm10", "end_mm10")
colnames(mouse_hg38)[1:3] <- c("chr_hg38", "start_hg38", "end_hg38")

mouse_enh <- merge(mouse_mm10, mouse_hg38, by="V4")

# genomes
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
mm10 <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10


# reading chipseq datasets

Hnf4a_mouse <- read.table(file = Hnf4a_mm10.bed, header = F, stringsAsFactors = F, sep ='\t')
CEBPA_mouse <- read.table(file = CEBPA_mm10.bed, header = F, stringsAsFactors = F, sep ='\t')

# GenomicRanges objects
CEBPA_mouse_gr <- with(CEBPA_mouse, GRanges(V1, IRanges(V2+1, V3)))
Hnf4a_mouse_gr <- with(Hnf4a_mouse, GRanges(V1, IRanges(V2+1, V3)))

## Overlapping mousee candidate enhancers with Chipseq datasets
mouse_enh_gr <- with(mouse_enh, GRanges(chr_mm10, IRanges(start_mm10+1, end_mm10)))
x <- as.data.frame(findOverlaps(mouse_enh_gr, CEBPA_mouse_gr))
mouse_ov_CEBPA <- cbind(mouse_enh[x$queryHits,], CEBPA_mouse[x$subjectHits,])
x <- as.data.frame(findOverlaps(mouse_enh_gr, Hnf4a_mouse_gr))
mouse_ov_Hnf4a <- cbind(mouse_enh[x$queryHits,], Hnf4a_mouse[x$subjectHits,])

# Saving candidate enhancers (mm10) and non-functional aligned regions (hg38)

tmp <- getSeq(hg38, with(unique(mouse_ov_CEBPA[,c("chr_hg38", "start_hg38", "end_hg38")])
                         , GRanges(chr_hg38, IRanges(start_hg38+1, end_hg38))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(mouse_ov_CEBPA[,c("chr_hg38", "start_hg38", "end_hg38")])
                           , paste(chr_hg38, paste(start_hg38, end_hg38, sep ='-'), sep = ":"))
            , file.out = "hg38_non_functional_ov_mm10_CEBPA.fa")

tmp <- getSeq(mm10, with(unique(mouse_ov_CEBPA[,c("chr_mm10", "start_mm10", "end_mm10")])
                         , GRanges(chr_mm10, IRanges(start_mm10+1, end_mm10))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(mouse_ov_CEBPA[,c("chr_mm10", "start_mm10", "end_mm10")])
                           , paste(chr_mm10, paste(start_mm10, end_mm10, sep ='-'), sep = ":"))
            , file.out = "mm10_enh_ov_CEBPA.fa")


tmp <- getSeq(hg38, with(unique(mouse_ov_Hnf4a[,c("chr_hg38", "start_hg38", "end_hg38")])
                         , GRanges(chr_hg38, IRanges(start_hg38+1, end_hg38))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(mouse_ov_Hnf4a[,c("chr_hg38", "start_hg38", "end_hg38")])
                           , paste(chr_hg38, paste(start_hg38, end_hg38, sep ='-'), sep = ":"))
            , file.out = "hg38_non_functional_ov_mm10_Hnf4a.fa")

tmp <- getSeq(mm10, with(unique(mouse_ov_Hnf4a[,c("chr_mm10", "start_mm10", "end_mm10")])
                         , GRanges(chr_mm10, IRanges(start_mm10+1, end_mm10))))
write.fasta(sequences = as.list(as.character(tmp))
            , names = with(unique(mouse_ov_Hnf4a[,c("chr_mm10", "start_mm10", "end_mm10")])
                           , paste(chr_mm10, paste(start_mm10, end_mm10, sep ='-'), sep = ":"))
            , file.out = "mm10_enh_ov_Hnf4a.fa")

