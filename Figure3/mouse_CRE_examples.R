# Script used to select final examples of mouse candidate enhancers that contain TF binding motifs (FIMO) and are predicted to have a CEBPA or HNF4A binding site by the the Domain Adaptive model

library(data.table)
library("GenomeInfoDb")
library("GenomicFeatures")
library("AnnotationFilter")
library(GenomicRanges)
library(RColorBrewer)

# Reading results from FIMO
Cebpa_human <- read.table(file = "./fimo_out/hg38_non_functional_Cebpa/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')#177  10
Cebpa_mouse <- read.table(file = "./fimo_out/mm10_enh_ov_CEBPA_Cebpa/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')#196  10
Hnf4a_human <- read.table(file = "./fimo_out/hg38_non_functional_Hnf4a/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')#1324   10
Hnf4a_mouse <- read.table(file = "./fimo_out/mm10_enh_ov_Hnf4a_Hnf4a/fimo.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')# 1032   10

# order motif matches by score
Cebpa_mouse <- Cebpa_mouse[with(Cebpa_mouse, order(-score)), ]
Hnf4a_mouse <- Hnf4a_mouse[with(Hnf4a_mouse, order(-score)), ]

# top motif instances
Cebpa_mouse[1, ]
# motif_id motif_alt_id           sequence_name start stop strand   score
# 1 M00993_2.00        Cebpa chr19:37702553-37703638   255  263      + 13.9515
# p.value q.value matched_sequence
# 1 6.62e-06       1        CTTGCGCAA

Cebpa_mouse[2, ] # motif location
# motif_id motif_alt_id           sequence_name start stop strand   score
# 2 M00993_2.00        Cebpa chr14:30860809-30862925  1037 1045      + 13.9515
# p.value q.value matched_sequence
# 2 6.62e-06       1        GTTGCGCAA

Hnf4a_mouse[1:2, ]
# motif_id motif_alt_id             sequence_name start stop strand   score
# 1 M00182_2.00        Hnf4a chr11:106056870-106058963   216  225      - 8.44048
# 2 M00182_2.00        Hnf4a    chr9:63200369-63202508   386  395      + 8.44048
# p.value q.value matched_sequence
# 1 1.12e-06   0.454       AGGGGCCACC
# 2 1.12e-06   0.454       AGGGTCCACC

# Making a data frame with the top sequences above  
sample_mouse_enh <- data.frame(chr = c("chr11", "chr19", "chr14", "chr9")
                               , start = c(106056870, 37702553, 30860809, 63200369)
                               , end = c(106058963, 37703638, 30862925, 63202508))


# Reading prediction results from the Domain Adaptive model
da_predicted <- read.table(file = "./results/result_data1_fixed/add_predictions_annot"
                           , header = T, stringsAsFactors = F, sep ='\t')

# seubset to sequences that were predicted to have a TF binding site
da_predicted_class1 <- da_predicted[da_predicted$predicted_class == 1, ]



### the original sequence coordinates are in the following files:
### read sequences and pair human-mouse
human_hg38 <- fread("Hsap_Enhancers_not_Mmus_alignMmus_hg38.bed", header=F)#actual human enhancers
human_mm10 <- fread("Hsap_Enhancers_not_Mmus_alignMmus_mm10.bed", header=F)#not functional in mouse

mouse_hg38 <- fread("Mmus_Enhancers_not_Hsap_alignHsap_hg38.bed", header=F)#not functional in mouse
mouse_mm10 <- fread("Mmus_Enhancers_not_Hsap_alignHsap_mm10.bed", header=F)#actual enhancer mouse

#get central 500bp
# Domain Adaptive model was trained on 500bp sequences

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

# Add DA model predicted probabilities to each group of sequences
# score > 0.9 indicated the sequence is predicted to have a TF binding site
human_hg38_pred1 <- merge(human_hg38_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id") # candidate human enh
mouse_hg38_pred1 <- merge(mouse_hg38_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id") # non-functional in mouse
human_mm10_pred1 <- merge(human_mm10_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id") # non-functional in human
mouse_mm10_pred1 <- merge(mouse_mm10_500, da_predicted_class1[,c("id", "test_set", "prob", "predicted_class")]
                          , by="id") # candidate mouse enh


# add id for the full sequence coordinates
human_hg38_pred1$id_fullSeq <- with(human_hg38_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))
mouse_hg38_pred1$id_fullSeq <- with(mouse_hg38_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))
human_mm10_pred1$id_fullSeq <- with(human_mm10_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))
mouse_mm10_pred1$id_fullSeq <- with(mouse_mm10_pred1, paste(V1, paste(V2, V3, sep ='-'), sep = ":"))

# id to sample sequences
sample_mouse_enh$id <- with(sample_mouse_enh, paste(chr, paste(start, end, sep = "-"), sep =":"))

# 2 of the selected candidate mouse enhancers pass the Domain Adaptive model threshold (probability >= 0.9)
sample_mouse_enh$id[sample_mouse_enh$id %in% mouse_mm10_pred1$id_fullSeq]
# [1] "chr14:30860809-30862925" "chr9:63200369-63202508"
