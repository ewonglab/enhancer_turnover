# Script used to select the final examples of human candidate enhancers and their aligning regions in the mouse genome.
# The selected candidate enhancers pass the Domain Adaptive prediction threshold (>= 0.9) and contain motif instances (CEBPA or HNF4A) identified by FIMO


library(Biostrings)
library("seqinr")
library(data.table)

# Reading results from FIMO to select sequences based on motif score (human candidate enhancers and mouse non-functional)
Cebpa_mouse <- read.table(file = "./fimo_out/mm10_non_functional_Cebpa/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')
Cebpa_human <- read.table(file = "./fimo_out/hg38_enh_ov_CEBPA_Cebpa/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')
Hnf4a_mouse <- read.table(file = "./fimo_out/mm10_non_functional_Hnf4a/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')
Hnf4a_human <- read.table(file = "./fimo_out/hg38_enh_ov_Hnf4a_Hnf4a/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')


# Order motif matches by score
Cebpa_human <- Cebpa_human[with(Cebpa_human, order(-score)), ]

Hnf4a_human <- Hnf4a_human[with(Hnf4a_human, order(-score)), ]


# Read prediction results from Domain Adaptive model
da_predicted <- read.table(file = "./results/result_data1_fixed/add_predictions_annot"
                           , header = T, stringsAsFactors = F, sep ='\t')
# Subset to keep those predicted to have a TF binding site
da_predicted_class1 <- da_predicted[da_predicted$predicted_class == 1, ]


# Reading coordinates of candidate enhancers and non-functional regions
human_hg38 <- fread("Hsap_Enhancers_not_Mmus_alignMmus_hg38.bed", header=F)# candidate human enhancers
human_mm10 <- fread("Hsap_Enhancers_not_Mmus_alignMmus_mm10.bed", header=F)# not functional in mouse

mouse_hg38 <- fread("Mmus_Enhancers_not_Hsap_alignHsap_hg38.bed", header=F)# not functional in mouse
mouse_mm10 <- fread("Mmus_Enhancers_not_Hsap_alignHsap_mm10.bed", header=F)# candidate enhancer mouse

# Get central 500bp
# The Domain Adaptive model was trained and tested on regions of 500bp

peak_centre <- function(x){
  x$width <- with(x, V3-V2)
  x$peak_center <- with(x, V2+(width/2))
  x$start <- with(x, round(peak_center-250))
  x$end <- with(x, round(peak_center+250))
  x$width <- with(x, end - start)
  print(unique(x$width))
  return(x)
}


human_hg38_500 <- peak_centre(human_hg38)
# [1] 500
mouse_hg38_500 <- peak_centre(mouse_hg38)
# [1] 500


human_mm10_500 <- peak_centre(human_mm10)
# [1] 500
mouse_mm10_500 <- peak_centre(mouse_mm10)
# [1] 500


#  V1        V2        V3 contain the coordinates for the full sequence
# start       end contain the coordinates for the trimmed sequence

human_hg38_500$id <- with(human_hg38_500, paste(V1, paste(start, end, sep ='-'), sep = ":"))
mouse_hg38_500$id <- with(mouse_hg38_500, paste(V1, paste(start, end, sep ='-'), sep = ":"))
human_mm10_500$id <- with(human_mm10_500, paste(V1, paste(start, end, sep ='-'), sep = ":"))
mouse_mm10_500$id <- with(mouse_mm10_500, paste(V1, paste(start, end, sep ='-'), sep = ":"))
da_predicted_class1$id <- with(da_predicted_class1, paste(chr, paste(start, end, sep ='-'), sep = ":"))

# Add prediction results to every set of regions
human_hg38_pred1 <- merge(human_hg38_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id")# candidate human enh
mouse_hg38_pred1 <- merge(mouse_hg38_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id")  # non functional in mouse
human_mm10_pred1 <- merge(human_mm10_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id") # non functional in human
mouse_mm10_pred1 <- merge(mouse_mm10_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id") # candidate mouse enh

# Add id for the full sequence coordinates
human_hg38_pred1$id_fullSeq <- with(human_hg38_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))
mouse_hg38_pred1$id_fullSeq <- with(mouse_hg38_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))
human_mm10_pred1$id_fullSeq <- with(human_mm10_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))
mouse_mm10_pred1$id_fullSeq <- with(mouse_mm10_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))


## Top 1 human CEBPA
head(Cebpa_human[Cebpa_human$sequence_name %in% human_hg38_pred1$id_fullSeq,],1)
# motif_id motif_alt_id        sequence_name start stop strand   score
# 1 M04045_2.00        CEBPA chr1:2493678-2495554   779  790      - 13.9598
# p.value q.value matched_sequence
# 1 7.94e-07       1     GCTTGCGCAACC

# Candidate enhancer passes Domain Adaptive Threshold
"chr1:2493678-2495554" %in% human_hg38_pred1$id_fullSeq 
# TRUE


# 2 hits of the motif in the same enhancer (same position, both strands)
Cebpa_human[Cebpa_human$sequence_name == "chr1:2493678-2495554",]
# motif_id motif_alt_id        sequence_name start stop strand   score
# 1 M04045_2.00        CEBPA chr1:2493678-2495554   779  790      - 13.9598
# 4 M04045_2.00        CEBPA chr1:2493678-2495554   779  790      + 13.5862
# p.value q.value matched_sequence
# 1 7.94e-07       1     GCTTGCGCAACC
# 4 2.05e-06       1     GGTTGCGCAAGC



head(Hnf4a_human[Hnf4a_human$sequence_name %in% human_hg38_pred1$id_fullSeq,],1)
# motif_id motif_alt_id            sequence_name start stop strand   score
# 4 M03357_2.00        HNF4A chr2:120515719-120518417  1155 1170      - 19.6121
# p.value q.value matched_sequence
# 4 1.06e-07   0.208 ATTGGACTTTGAACTG

# candidate enhancer in top 4 passes Domain Adaptive probability threshold (>= 0.9)
which(Hnf4a_human$sequence_name == "chr2:120515719-120518417")# top 4


# only one motif instance
Hnf4a_human[Hnf4a_human$sequence_name == "chr2:120515719-120518417",]
# motif_id motif_alt_id            sequence_name start stop strand   score
# 4 M03357_2.00        HNF4A chr2:120515719-120518417  1155 1170      - 19.6121
# p.value q.value matched_sequence
# 4 1.06e-07   0.208 ATTGGACTTTGAACTG


# Matching again human and mouse regions to select the non-functional aligning regions in mouse genome
# reading enhancers and aligned non-functional elements
human_mm10 <- read.table(file = "Hsap_Enhancers_not_Mmus_alignMmus_mm10.bed"
                         , header = F, stringsAsFactors = F, sep ='\t')# non-functional
human_hg38 <- read.table(file = "Hsap_Enhancers_not_Mmus_alignMmus_hg38.bed"
                         , header = F, stringsAsFactors = F, sep ='\t')# human candidate enhancers

colnames(human_mm10)[1:3] <- c("chr_mm10", "start_mm10", "end_mm10")
colnames(human_hg38)[1:3] <- c("chr_hg38", "start_hg38", "end_hg38")

human_enh <- merge(human_hg38, human_mm10, by="V4")

# add mouse and human region IDs
human_enh$hg38_id <- with(human_enh, paste(chr_hg38, paste(start_hg38, end_hg38, sep = '-'), sep = ':'))
human_enh$mm10_id <- with(human_enh, paste(chr_mm10, paste(start_mm10, end_mm10, sep = '-'), sep = ':'))


# pair of sequences #1
human_enh[human_enh$hg38_id == "chr1:2493678-2495554", c("hg38_id", "mm10_id")]
#                   hg38_id                  mm10_id
# 3708 chr1:2493678-2495554 chr4:154992548-154994291


# pair of sequences #2
human_enh[human_enh$hg38_id == "chr2:120515719-120518417", c("hg38_id", "mm10_id")]
#                       hg38_id                  mm10_id
# 4794 chr2:120515719-120518417 chr1:119246593-119250914

