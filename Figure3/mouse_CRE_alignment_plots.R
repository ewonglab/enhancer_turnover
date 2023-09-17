library("ggmsa")
library("Biostrings")

# Alignments from clustalo
chr14_enh <- "align2.fasta"
chr9_enh <- "align4.fasta"

chr14_seq <- readDNAStringSet(chr14_enh)
chr9_seq <- readDNAStringSet(chr9_enh)

# Function to transform motif locations in non-gapped sequence to motif location in alignment (gapped sequence)

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

# Cebpa_mouse[2, ] # motif location from FIMO output
# motif_id motif_alt_id           sequence_name start stop strand   score
# 2 M00993_2.00        Cebpa chr14:30860809-30862925  1037 1045      + 13.9515
# p.value q.value matched_sequence
# 2 6.62e-06       1        GTTGCGCAA

# find motif location in the aligned sequence
find_motif_location(as.character(chr14_seq[1]), 1037, 1045)
# $start
# [1] 1243
#
# $end
# [1] 1251
# there is another instance in the - strand

find_motif_location(as.character(chr14_seq[1]), 1038, 1046)
# $start
# [1] 1244
#
# $end
# [1] 1252


pdf("chr14_30860809_30862925_clustalo.pdf")
ggmsa(chr14_enh, start = 1210, end = 1260, char_width = 0.5, seq_name = T, color="Chemistry_NT") +
  geom_seqlogo() + geom_msaBar()
dev.off()

# Hnf4a_mouse[1:2, ]
# motif_id motif_alt_id             sequence_name start stop strand   score
# 1 M00182_2.00        Hnf4a chr11:106056870-106058963   216  225      - 8.44048
# 2 M00182_2.00        Hnf4a    chr9:63200369-63202508   386  395      + 8.44048
# p.value q.value matched_sequence
# 1 1.12e-06   0.454       AGGGGCCACC
# 2 1.12e-06   0.454       AGGGTCCACC

find_motif_location(as.character(chr9_seq[1]), 386, 395)
# $start
# [1] 420
#
# $end
# [1] 429

find_motif_location(as.character(chr9_seq[1]), 1602, 1611)
# $start
# [1] 1709
# $end
# [1] 1718

pdf("chr9_63200369_63202508_clustalo.pdf")
ggmsa(chr9_enh, start = 400, end = 450, char_width = 0.5, seq_name = T, color="Chemistry_NT") +
  geom_seqlogo() + geom_msaBar()
dev.off()
