library(GenomicRanges)
library(dplyr)

add.motifs.to.comp_JASPAR <- function(motifs.tab, motifs.li, comp.df, enh_type){
  motifs.tab <- motifs.tab[,colnames(motifs.tab)[!colnames(motifs.tab) %in%
                                                   c("type", "rel_rt")]]
  early <- motifs.li[["early"]]
  late <- motifs.li[["late"]]
  #getting rows of motif cound matrix for early and late enhancers
  early_motif <- subset(motifs.tab, row.names(motifs.tab) %in% early)
  print(dim(early_motif))
  late_motif <- subset(motifs.tab, row.names(motifs.tab) %in% late)
  print(dim(late_motif))
  #sum motif frequency across enhancers (early and late)
  early_motif_sum <- apply(early_motif, 2, sum)
  late_motif_sum <- apply(late_motif, 2, sum)
  #sum motif frequency across ALL enhancers
  motif_sum <- apply(motifs.tab, 2, sum)
  #number of enhancers with non-zero motif counts
  early_motif_non_zeros <- colSums(early_motif != 0)##
  late_motif_non_zeros <- colSums(late_motif != 0)#
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
  comp.df$late_motif_norm <- comp.df$late_motif_sum/ length(late) * 20
  comp.df$early_motif_norm <- comp.df$early_motif_sum/ length(early) *20
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
homer <- read.table(file = "liver_homeodomain_higherQuintile_humanJAS.txt"
                    , header =T, stringsAsFactors = F) 

set.seed(1)
late_sample <- sample(x = rownames(homer[homer$rel_rt == "late",])
                      , size = 5513, replace = F)

homer <- homer[rownames(homer) %in% late_sample | homer$rel_rt == "early",]

counts_sd <- apply(homer[,colnames(homer)[!colnames(homer) %in% c("type", "rel_rt")]]
                  , 1, sd)
table(counts_sd > 0)
# TRUE
# 11026

## read jaspar annotation and keep only human motifs
annot <- read.table(file = "jaspar_motifs_TFs_annot_all.txt"
                    , header =T, stringsAsFactors = F, sep = '\t')
annot <- annot[grep("Homo", annot$Species),]

# READING JASPAR PWM NUCLEOTIDE COMPOSITION
comp <- read.delim('JASPAR_core_nucl_comp_and_ic.txt', header =T)

#add motif numbers to nucleotide composition and save
comp.full <- add.motifs.to.comp_JASPAR(motifs.tab = homer
                                       , motifs.li = list(early = rownames(homer[homer$rel_rt == "early",])
                                                          , late = rownames(homer[homer$rel_rt == "late",]))
                                       , comp.df = comp
                                       , enh_type = "all")

#keeping only jaspar human
comp.full <- comp.full[comp.full$pwm %in% annot$ID,]
write.table(comp.full, 'motif_comp_liver_homeodomain_higherQuintile_humanJAS.txt'
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
  comp.df <- merge(comp.df, df, by.x='pwm', by.y=0)
  #plot earlyvlate with % AT %GC
  df <- comp.df[,c('earlyvlate','C.proportion','G.proportion'
                   ,'GC.proportion','A.proportion','T.proportion','AT.proportion'
                   , "Family", "Class", "Species")]
  row.names(df) <- comp.df$pwm
  df <- df[is.finite(rowSums(df[,c(1:7)])),]
  print(nrow(df))
  df <- df[order(df$earlyvlate), ]
  return(df)
}

comp <- read.table(file = "motif_comp_liver_homeodomain_higherQuintile_humanJAS.txt"
                   , header = T, stringsAsFactors = F, sep = '\t')# 810  30
### ALL
all.heat <-  data.4.heatmap_jas(comp, "all")

write.table(x = all.heat
            , file = "human_liver_homeodomain_highestQuintile_humanJAS_heat"
            , quote = F, sep = '\t')
