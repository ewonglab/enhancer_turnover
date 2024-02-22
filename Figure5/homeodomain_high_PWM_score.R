library(stringr)

setwd("./replication_timing/human")

jaspar <- read.table(file = "jaspar_motifs_TFs_annot_all.txt"
                     , header =T, stringsAsFactors = F
                     , sep = '\t')
homeodomain <- jaspar[jaspar$Class == "Homeo domain factors",]

# read motif location and scores

cons_early <- read.table(file = "./homer/liver_non.recent_early_quintile_H9RT_ext.humanJAS.bed"
                         , header = F, stringsAsFactors = F, sep ='\t', skip = 1) 
cons_late <- read.table(file = "./homer/liver_non.recent_late_quintile_H9RT_ext.humanJAS.bed"
                        , header = F, stringsAsFactors = F, sep ='\t', skip = 1) 
recent_early <- read.table(file = "./homer/liver_recent_early_quintile_H9RT_ext.humanJAS.bed"
                           , header = F, stringsAsFactors = F, sep ='\t', skip = 1) 
recent_late <- read.table(file = "./homer/liver_recent_late_quintile_H9RT_ext.humanJAS.bed"
                          , header = F, stringsAsFactors = F, sep ='\t', skip = 1) 

# keep only homeodomain factors
cons_early.homeo <- cons_early[cons_early$V4 %in% homeodomain$ID, ] 
cons_late.homeo <- cons_late[cons_late$V4 %in% homeodomain$ID, ] 
recent_early.homeo <- recent_early[recent_early$V4 %in% homeodomain$ID, ] 
recent_late.homeo <- recent_late[recent_late$V4 %in% homeodomain$ID, ] 

# add type
cons_early.homeo$type <- "conserved"
cons_late.homeo$type <- "conserved"
recent_early.homeo$type <- "recent"
recent_late.homeo$type <- "recent"

cons_early.homeo$rel_rt <- "early"
cons_late.homeo$rel_rt <- "late"
recent_early.homeo$rel_rt <- "early"
recent_late.homeo$rel_rt <- "late"

homeodomain_motifs <- rbind(cons_early.homeo, cons_late.homeo
                            , recent_early.homeo, recent_late.homeo) 

summary(homeodomain_motifs$V5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 6.073   7.022   7.444   7.618   8.011  17.030

library("Hmisc")
library(GenomicRanges)

homeodomain_motifs$scorebin <- as.factor(cut2(homeodomain_motifs$V5, g=5))

idkey <- data.frame(bin = names(table(homeodomain_motifs$scorebin))
                    , quintile=c(1,2,3,4,5))
idkey
# bin quintile
# 1 [6.07, 6.92)        1
# 2 [6.92, 7.28)        2
# 3 [7.28, 7.64)        3
# 4 [7.64, 8.18)        4
# 5 [8.18,17.03]        5

homeodomain_motifs <- merge(homeodomain_motifs,  idkey
                            ,by.x='scorebin', by.y='bin', all.x=TRUE)
table(homeodomain_motifs$quintile)
# 1      2      3      4      5
# 322461 322165 321342 321899 321298

top_quintile <- homeodomain_motifs[homeodomain_motifs$quintile == 5, ]
# number of unique motifs
length(unique(top_quintile$V4))
# [1] 139

# read enhancers
bed_path <- "/liver/feb2022/"
 
cons_early_bed <- read.table(file = paste0(bed_path, "liver_non.recent_early_quintile_H9RT_ext.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')
cons_late_bed <- read.table(file = paste0(bed_path, "liver_non.recent_late_quintile_H9RT_ext.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')
recent_early_bed <- read.table(file = paste0(bed_path, "liver_recent_early_quintile_H9RT_ext.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')
recent_late_bed <- read.table(file = paste0(bed_path, "liver_recent_late_quintile_H9RT_ext.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')

cons_early_bed$type <- "conserved"
cons_late_bed$type <- "conserved"
recent_early_bed$type <- "recent"
recent_late_bed$type <- "recent"

cons_early_bed$rel_rt <- "early"
cons_late_bed$rel_rt <- "late"
recent_early_bed$rel_rt <- "early"
recent_late_bed$rel_rt <- "late"

all_enh <- rbind(cons_early_bed, cons_late_bed
                 , recent_early_bed, recent_late_bed)

# gr objects
top_quintile.gr <- with(data = top_quintile, expr = GRanges(V1, IRanges(V2 + 1, V3)))
all_enh.gr <- with(data = all_enh, expr = GRanges(V1, IRanges(V2 + 1, V3)))

x <- as.data.frame(findOverlaps(top_quintile.gr, all_enh.gr))
motifs_enh <- cbind(top_quintile[x$queryHits, "V4", drop=F]
                    , all_enh[x$subjectHits, ])
colnames(motifs_enh)[1] <- "ID"
motif_freq <- as.data.frame(table(motifs_enh[,c("V4", "ID")]))

library(tidyr)
motif_freq <- tidyr::spread(motif_freq, ID, Freq) 

motif_freq <- merge(motif_freq, all_enh[, c("V4", "type", "rel_rt")]
                    , by = 'V4', all.x=T)
rownames(motif_freq) <- motif_freq$V4
motif_freq$V4 <- NULL

write.table(x = motif_freq, file = "liver_homeodomain_higherQuintile_humanJAS.txt"
            , quote = F, sep = '\t')
