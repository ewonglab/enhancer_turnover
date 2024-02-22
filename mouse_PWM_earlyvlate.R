library(GenomicRanges)
library(dplyr)

### add counts to table of motif composition
add.motifs.to.comp_JASPAR <- function(motifs.tab, motifs.li, comp.df, enh_type){
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



setwd("./replication_timing/mouse")
#all jaspar: redundant + non-redundant
early <- read.table(file = "mouse_recent_and_cons2sp_EARLY_humanJAS2020_homer.txt"
                    , header =T, stringsAsFactors = F) 
late <- read.table(file = "mouse_recent_and_cons2sp_LATE_humanJAS2020_homer.txt"
                   , header =T, stringsAsFactors = F) 

# enhancers
all_enh <- read.table(file = "mouse_enh_meanRT_by_germ_line_cons2sp.txt"
                      , header = T, stringsAsFactors = F, sep ='\t')
all_enh$start <- all_enh$start + 1 # to match HOMER id format
all_enh$id <- with(all_enh, paste(chr, paste(start, end, sep ='-'), sep =':'))
all_enh$id <- sub("chr", "", all_enh$id)

# add conservation
early <- merge(early, all_enh[,c("id", "age")], by.x = 0, by.y = "id", all.x = T)
late <- merge(late, all_enh[,c("id", "age")], by.x = 0, by.y = "id", all.x = T)


early_sd <- apply(early[,colnames(early)[!colnames(early) %in% c("Row.names", "age")]]
                  , 1, sd)
late_sd <- apply(late[,colnames(late)[!colnames(late) %in% c("Row.names", "age")]]
                 , 1, sd)
table(early_sd > 0)
# FALSE  TRUE
# 6 85174
table(late_sd > 0)
# TRUE
# 10898

early <- early[-which(early_sd == 0), ]# 85174   798

set.seed(1)
recent_sample <- sample(x = early[early$age == "Recent", "Row.names"]
                        , size = 8051, replace = F)
set.seed(1)
cons_sample <- sample(x = early[early$age == "Conserved", "Row.names"]
                      , size = 2847, replace = F)


# motif counts at early and late RT enhancers
counts <-  dplyr::bind_rows(early[early$Row.names %in% c(recent_sample, cons_sample), ]
                            , late)
counts[is.na(counts)] <- 0
rownames(counts) <- counts$Row.names
counts$Row.names <- NULL

## read jaspar annotation and keep only human PWMs
annot <- read.table(file = "/g/data/zk16/cc3704/GBM/jaspar_motifs_TFs_annot_all.txt"
                    , header =T, stringsAsFactors = F, sep = '\t')
annot <- annot[grep("Homo", annot$Species),]

#calculate sd by enhancer and remove enhancers with sd = 0
counts_sd <- apply(counts[,colnames(counts)[colnames(counts)!="age"]]
                   , 1, sd)
range(counts_sd)#
# [1]  0.03544406 13.28639822

write.table(x = counts, file = "mouse_enh_earlyvlate_humanJASPAR_cons2sp.txt"
            , quote = F, sep = '\t')


# READING JASPAR PWM NUCLEOTIDE COMPOSITION
comp <- read.delim('JASPAR_core_nucl_comp_and_ic.txt'
                   , header =T)

# add motif numbers to nucleotide composition and save
counts$age <- NULL
comp.full <- add.motifs.to.comp_JASPAR(motifs.tab = counts
                                       , motifs.li = list(early = early$Row.names
                                                          , late = late$Row.names)
                                       , comp.df = comp
                                       , enh_type = "all")

# keeping only jaspar human
comp.full <- comp.full[comp.full$pwm %in% annot$ID,]# 810  30
write.table(comp.full, 'motif_comp_mouse_enh_ov_germRT_allTissues_humanJAS_cons2sp.txt'
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

comp <- read.table(file = "motif_comp_mouse_enh_ov_germRT_allTissues_humanJAS_cons2sp.txt"
                   , header = T, stringsAsFactors = F, sep = '\t')
### ALL
all.heat <-  data.4.heatmap_jas(comp, "all")


write.table(x = all.heat
            , file = "mouse_all_tissues_germ.rt_humanJAS_cons2sp_heat"
            , quote = F, sep = '\t')

library(pheatmap)                  
library(RColorBrewer)

color.palette <- function(steps, n.steps.between=NULL, ...){
 
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
 
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
 
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
 
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}

setwd("./replication_timing/mouse")
mouse <- read.table(file = "mouse_all_tissues_germ.rt_humanJAS_cons2sp_heat"
                    , header = T, stringsAsFactors = F, sep = '\t')# already ordered by earlyvlate

steps <- rev(brewer.pal(n = 8, name = "PuOr"))
pal <- color.palette(steps, c(20, 20,5,5,5, 20,20), space="rgb")

# all enhancers
pdf('heatmap.earlylate.mouse_allTissues_humanJAS_cons2sp.pdf')#
pheatmap(mouse[,c(2:7)], scale='column'
         , show_rownames = FALSE
         , color = pal(500)
         , show_colnames=TRUE
         ,  main="Mouse all tissues earlyvlate",
         cluster_rows=FALSE, cluster_cols=F, fontsize=5)
dev.off()
