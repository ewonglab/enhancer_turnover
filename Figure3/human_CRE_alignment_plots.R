# Script used to plot segments of the selected human candidate enhancers and mouse non-functional regions 
# The segments show in the plots include the CEBPA or HNF4A motif locations identified by FIMO

library("ggmsa")
library("Biostrings")

# Reading alignments for selected human candidate enhancers and mouse non--functional regions

chr1_enh <- "align5.fasta"
chr2_enh <- "align6.fasta"

chr1_seq <- readDNAStringSet(chr1_enh)
chr2_seq <- readDNAStringSet(chr2_enh)


# Function to identify TF binding motif locations in a gapped sequence (alignment) given the motif location in the non-gapped sequence (from FIMO).

find_motif_location <- function(gapped_sequence, original_motif_start, original_motif_end) {
  ungapped_sequence <- gsub("-", "", gapped_sequence)  # Remove gaps from gapped sequence
  
  if (original_motif_start > nchar(ungapped_sequence) || original_motif_end > nchar(ungapped_sequence) || original_motif_start > original_motif_end) {
    return(NULL)  # Invalid motif location
  }
  
  gapped_motif_start <- 0
  gapped_motif_end <- 0
  gapped_index <- 0
  
  for (i in 1:nchar(gapped_sequence)) {
    char <- substr(gapped_sequence, i, i)
    if (char != '-') {
      gapped_index <- gapped_index + 1
      if (gapped_index == original_motif_start) {
        gapped_motif_start <- i
      }
      if (gapped_index == original_motif_end) {
        gapped_motif_end <- i
        break
      }
    }
  }
  
  return(list(start = gapped_motif_start, end = gapped_motif_end))
}


# find CEBPA motif location in the aligned sequence using the location in the non-gapped sequence (from FIMO output) 
# Cebpa_human[1,]
# motif_id motif_alt_id        sequence_name start stop strand   score
# 1 M04045_2.00        CEBPA chr1:2493678-2495554   779  790      - 13.9598
# p.value q.value matched_sequence
# 1 7.94e-07       1     GCTTGCGCAACC

find_motif_location(as.character(chr1_seq[1]), 779, 790)# 2 hits both strands same position
# $start
# [1] 801
#
# $end
# [1] 812

# I am plotting a region that includes the location of the CEBPA motif
pdf("chr1_2493678_2495554_hg38_clustalo.pdf")
ggmsa(chr1_enh, start = 790, end = 840, char_width = 0.5, seq_name = T, color="Chemistry_NT") +
  geom_seqlogo() + geom_msaBar()
dev.off()

# Getting the location of the HNF4A motif in the aligned sequence
# Hnf4a_human[Hnf4a_human$sequence_name == "chr2:120515719-120518417",]
# motif_id motif_alt_id            sequence_name start stop strand   score
# 4 M03357_2.00        HNF4A chr2:120515719-120518417  1155 1170      - 19.6121
# p.value q.value matched_sequence
# 4 1.06e-07   0.208 ATTGGACTTTGAACTG

find_motif_location(as.character(chr2_seq[1]), 1155, 1170)
# $start
# [1] 1503
#
# $end
# [1] 1519

pdf("chr2_120515719_120518417_hg38_clustalo.pdf")
ggmsa(chr2_enh, start = 1490, end = 1540, char_width = 0.5, seq_name = T, color="Chemistry_NT") +
  geom_seqlogo() + geom_msaBar()
dev.off()


