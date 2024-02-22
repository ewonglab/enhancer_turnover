library(GenomicRanges)
library(dplyr)

### add counts to table of motif composition
add.motifs.to.comp_JASPAR <- function(motifs.tab, motifs.li, comp.df, enh_type){
  motifs.tab <- motifs.tab[,colnames(motifs.tab)[!colnames(motifs.tab) %in%
                                                   c("age", "rel_rt")]]
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



#all jaspar: redundant + non-redundant
homer <- read.table(file = "allTissues_homeodomain_topQuintile_humanJAS.txt"
                    , header =T, stringsAsFactors = F) # 68729   139

# enhancers
all_enh <- read.table(file = "mouse_enh_meanRT_by_germ_line_cons2sp.txt"
                      , header = T, stringsAsFactors = F, sep ='\t')#  134195     11
all_enh$id <- with(all_enh, paste(chr, paste(start, end, sep ='-'), sep =':'))
all_enh$id <- sub("chr", "", all_enh$id)

table(rownames(homer) %in% (all_enh$id))
# TRUE
# 68729

# add conservation
homer <- merge(homer, all_enh[,c("id", "age")], by.x = 0, by.y = "id", all.x = T)

table(homer$age)
# Conserved    Recent
# 39925     28804

table(homer$rel_rt)
# early  late
# 61071  7658

early_freq <- table(homer[homer$rel_rt == "early", "age"])
# Conserved    Recent
#     37680     23391
late_freq <- table(homer[homer$rel_rt == "late", "age"])
# Conserved    Recent
#     2245      5413

n_cons <- min(c(early_freq["Conserved"], late_freq["Conserved"]))
n_recent <- min(c(early_freq["Recent"], late_freq["Recent"]))

n_cons
# [1] 2245
n_recent
# [1] 5413 # subsample from early enhancers

set.seed(1)
early.recent_sample <- sample(x = homer[homer$rel_rt == "early" &
                                          homer$age == "Recent", "Row.names"]
                              , size = n_recent, replace = F)
set.seed(1)
early.cons_sample <- sample(x = homer[homer$rel_rt == "early" &
                                        homer$age == "Conserved", "Row.names"]
                            , size = n_cons, replace = F)

homer_sub <- homer[homer$rel_rt == "late" |
                     homer$Row.names %in% c(early.recent_sample, early.cons_sample),]# 15316   141

as.data.frame(table(homer_sub[,c("age", "rel_rt")]))
#         age rel_rt Freq
# 1 Conserved  early 2245
# 2    Recent  early 5413
# 3 Conserved   late 2245
# 4    Recent   late 5413

as.data.frame(table(homer_sub[,c("rel_rt")]))
# Var1 Freq
# 1 early 7658
# 2  late 7658

homer_sub_sd <- apply(homer_sub[,colnames(homer_sub)[!colnames(homer_sub) %in%
                                                       c("Row.names", "age", "rel_rt")]]
                      , 1, sd)
table(homer_sub_sd > 0)
# TRUE
# 15316


## read jaspar annotation and keep only human motifs
annot <- read.table(file = "/g/data/zk16/cc3704/GBM/jaspar_motifs_TFs_annot_all.txt"
                    , header =T, stringsAsFactors = F, sep = '\t')
annot <- annot[grep("Homo", annot$Species),]#810   7

# READING JASPAR PWM NUCLEOTIDE COMPOSITION
comp <- read.delim('JASPAR_core_nucl_comp_and_ic.txt'
                   , header =T)# 1964   16

rownames(homer_sub) <- homer_sub$Row.names
homer_sub$Row.names <- NULL
#add motif numbers to nucleotide composition and save
comp.full <- add.motifs.to.comp_JASPAR(motifs.tab = homer_sub
                                       , motifs.li = list(early = rownames(homer_sub[homer_sub$rel_rt == "early",])
                                                          , late = rownames(homer_sub[homer_sub$rel_rt == "late",]))
                                       , comp.df = comp
                                       , enh_type = "all")
# [1] 7658  138
# [1] 7658  138

#keeping only jaspar human
comp.full <- comp.full[comp.full$pwm %in% annot$ID,]# 810  30
write.table(comp.full, 'motif_comp_allTissues_homeodomain_topQuintile_humanJAS.txt'
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

comp <- read.table(file = "motif_comp_allTissues_homeodomain_topQuintile_humanJAS.txt"
                   , header = T, stringsAsFactors = F, sep = '\t')# 810  30
### ALL
all.heat <-  data.4.heatmap_jas(comp, "all")


write.table(x = all.heat
            , file = "mouse_allTissues_homeodomain_topQuintile_humanJAS_heat"
            , quote = F, sep = '\t')
