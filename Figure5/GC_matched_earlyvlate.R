library(GenomicRanges)
library(dplyr)

add.motifs.to.comp_JASPAR <- function(motifs.tab, motifs.li, comp.df, enh_type){
  early <- motifs.li[["early"]]
  late <- motifs.li[["late"]]
  #getting rows of motif cound matrix for early and late enhancers
  early_motif <-subset(motifs.tab, row.names(motifs.tab) %in% early$V4)
  late_motif <-subset(motifs.tab, row.names(motifs.tab) %in% late$V4)
  #sum motif frequency across enhancers (early and late)
  early_motif_sum <- apply(early_motif, 2, sum)
  late_motif_sum <- apply(late_motif, 2, sum)
  #sum motif frequency across ALL enhancers
  motif_sum <- apply(motifs.tab, 2, sum)
  #number of enhancers with non-zero motif counts
  early_motif_non_zeros <- colSums(early_motif != 0)
  late_motif_non_zeros <- colSums(late_motif != 0)
  #adding motif count in early enhancers, late enhancers and all enhancers
  ncol.comp <- ncol(comp.df)
  comp.df$early_motif_sum <- early_motif_sum[match(comp.df$pwm, names(early_motif_sum))]
  comp.df$late_motif_sum <- late_motif_sum[match(comp.df$pwm, names(late_motif_sum))]
  comp.df$motif_sum <- motif_sum[match(comp.df$pwm, names(motif_sum))]
  #calculating proportion of motifs in early and late enhancers compared to all enhancers
  comp.df$late_motif_prop <- comp.df$late_motif_sum/comp.df$motif_sum
  comp.df$early_motif_prop <- comp.df$early_motif_sum/comp.df$motif_sum
  #add number of enhancers that contain every motif
  comp.df$early_motif_non_zeros <- early_motif_non_zeros[match(comp.df$pwm
                                                               , names(early_motif_non_zeros))]
  comp.df$late_motif_non_zeros <- late_motif_non_zeros[match(comp.df$pwm
                                                             , names(late_motif_non_zeros))]
  #Total number of early and late enhancers
  comp.df$all_peaks_early <- nrow(early_motif)
  comp.df$all_peaks_late <- nrow(late_motif)
  # normalize by number of peaks and rescale
  comp.df$late_motif_norm <- comp.df$late_motif_sum/ nrow(late) * 20
  comp.df$early_motif_norm <- comp.df$early_motif_sum/ nrow(early) *20
  comp.df$comp.dfare <- log2(comp.df$early_motif_norm/comp.df$late_motif_norm)#proportion early/late
  #number of early and late enhancers WHITHOUT motif
  comp.df$early_zero <- comp.df$all_peaks_early - comp.df$early_motif_non_zeros
  comp.df$late_zero <- comp.df$all_peaks_late - comp.df$late_motif_non_zeros
  comp.df <- comp.df[order(-comp.df$late_motif_prop),]
  colnames(comp.df)[(ncol.comp+1):ncol(comp.df)] <- paste(colnames(comp.df)[(ncol.comp+1):ncol(comp.df)]
                                                          , enh_type, sep = "_")
  return(comp.df)}

setwd("./replication_timing/human")
#all jaspar: redundant + non-redundant
recent <- read.table(file = "liver_recent_NULL_earlyvlate_jaspar.txt"
                     , header = T, stringsAsFactors = F)
non.recent <- read.table(file = "liver_non.recent_NULL_earlyvlate_jaspar.txt"
                         , header = T, stringsAsFactors = F)

bed_path <- "./replication_timing/human/"

recent.early.bed <- read.table(file = paste0(bed_path, "recent_early_quintile_H9RT_ext_NULL.bed")
                               , header =F, stringsAsFactors = F, sep = '\t')
recent.late.bed <- read.table(file = paste0(bed_path, "recent_late_quintile_H9RT_ext_NULL.bed")
                              , header =F, stringsAsFactors = F, sep = '\t')
non.recent.early.bed <- read.table(file = paste0(bed_path, "conserved_early_quintile_H9RT_ext_NULL.bed")
                                   , header =F, stringsAsFactors = F, sep = '\t')
non.recent.late.bed <- read.table(file = paste0(bed_path, "conserved_late_quintile_H9RT_ext_NULL.bed")
                                  , header =F, stringsAsFactors = F, sep = '\t')

recent.early.bed$V4 <- with(recent.early.bed, paste(V1, paste(V2 + 1, V3, sep = '-'), sep = ':'))
recent.late.bed$V4 <- with(recent.late.bed, paste(V1, paste(V2 + 1, V3, sep = '-'), sep = ':'))
non.recent.early.bed$V4 <- with(non.recent.early.bed, paste(V1, paste(V2 + 1, V3, sep = '-'), sep = ':'))
non.recent.late.bed$V4 <- with(non.recent.late.bed, paste(V1, paste(V2 + 1, V3, sep = '-'), sep = ':'))
   
## jaspar annotation
annot <- read.table(file = "jaspar_motifs_TFs_annot_all.txt"
                    , header =T, stringsAsFactors = F, sep = '\t')
annot <- annot[grep("Homo", annot$Species),]

jaspar.counts <- dplyr::bind_rows(recent, non.recent)
jaspar.counts[is.na(jaspar.counts)] <- 0
jaspar.counts <- jaspar.counts[,c(colnames(jaspar.counts)[!colnames(jaspar.counts) %in%
                                                            c("is.early", "is.recent")]
                                  , "is.early", "is.recent")]


jaspar.counts <- jaspar.counts[,c(colnames(jaspar.counts)[colnames(jaspar.counts) %in% annot$ID]
                                  ,"is.early", "is.recent")]

# NUCLEOTIDE COMPOSITION OF JASPAR MOTIFS
comp <- read.delim('JASPAR_core_nucl_comp_and_ic.txt'
                   , header =T)

#add motif numbers to nucleotide composition and save
#5538 early enhancers and 5538 late enhancers
comp.full <- add.motifs.to.comp_JASPAR(motifs.tab = jaspar.counts
                                       , motifs.li = list(early = rbind(recent.early.bed, non.recent.early.bed)
                                                          , late = rbind(recent.late.bed, non.recent.late.bed))
                                       , comp.df = comp
                                       , enh_type = "liver")
#keeping only jaspar human
comp.full <- comp.full[comp.full$pwm %in% annot$ID,]
write.table(comp.full, 'motif_comp_humanLiver_human_NULL_jaspar.txt'
            , quote = F, sep = '\t')


data.4.heatmap_jas <- function(comp.df, enh_type){
  enh_cols <- paste(c('early_motif_non_zeros', 'early_zero',  'late_motif_non_zeros', 'late_zero')
                    , enh_type, sep = "_")
  df<- comp.df[,enh_cols]
  row.names(df) <- comp.df$pwm
  df <- df[complete.cases(df), ]
  print(nrow(df))
  comp.df$earlyvlate <- log2((comp.df[,enh_cols[1]] / comp.df[,enh_cols[2]]) /
                               (comp.df[,enh_cols[3]] / comp.df[,enh_cols[4]]) )
  comp.df <- merge(comp.df, df, by.x='pwm', by.y=0)#
  #plot earlyvlate with % AT %GC
  df <- comp.df[,c('earlyvlate','C.proportion','G.proportion'
                   ,'GC.proportion','A.proportion','T.proportion','AT.proportion'
                   , "Family", "Class", "Species")]
  row.names(df) <- comp.df$pwm
  df <- df[is.finite(rowSums(df[,c(1:7)])),]#
  print(nrow(df))
  df <- df[order(df$earlyvlate), ]
  return(df)
}

comp <- read.table(file = "motif_comp_humanLiver_human_NULL_jaspar.txt"
                   , header = T, stringsAsFactors = F, sep = '\t')

all.heat <-  data.4.heatmap_jas(comp, "liver")
write.table(x = all.heat, file = "human_liver_all_enh_NULL_rt_humanJaspar_heat"
            , quote = F, sep = '\t')

