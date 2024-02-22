library(stringr)

jaspar <- read.table(file = "jaspar_motifs_TFs_annot_all.txt"
                     , header =T, stringsAsFactors = F
                     , sep = '\t')
homeodomain <- jaspar[jaspar$Class == "Homeo domain factors",]

# read motif location and scores

early <- read.table(file = "mouse_recent_and_cons2sp_EARLY_human.jaspar.bed"
                         , header = F, stringsAsFactors = F, sep ='\t', skip = 1) 
late <- read.table(file = "mouse_recent_and_cons2sp_LATE_human.jaspar.bed"
                        , header = F, stringsAsFactors = F, sep ='\t', skip = 1) 

# keep only homeodomain factors
early.homeo <- early[early$V4 %in% homeodomain$ID, ] 
late.homeo <- late[late$V4 %in% homeodomain$ID, ] 

# add type
early.homeo$rel_rt <- "early"
late.homeo$rel_rt <- "late"

homeodomain_motifs <- rbind(early.homeo, late.homeo) 

summary(homeodomain_motifs$V5)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 6.073   7.027   7.446   7.624   7.986  17.847

# get the motifs on the highest quintile of scores
library("Hmisc")
library(GenomicRanges)

homeodomain_motifs$scorebin <- as.factor(cut2(homeodomain_motifs$V5, g=5))

idkey <- data.frame(bin = names(table(homeodomain_motifs$scorebin))
                    , quintile=c(1,2,3,4,5))
idkey
# bin quintile
# 1 [6.07, 6.93)        1
# 2 [6.93, 7.29)        2
# 3 [7.29, 7.65)        3
# 4 [7.65, 8.21)        4
# 5 [8.21,17.85]        5

homeodomain_motifs <- merge(homeodomain_motifs,  idkey
                            ,by.x='scorebin', by.y='bin', all.x=TRUE)
table(homeodomain_motifs$quintile)
# 1      2      3      4      5
# 666523 667973 665689 662870 664892

top_quintile <- homeodomain_motifs[homeodomain_motifs$quintile == 5, ]
# number of unique motifs
length(unique(top_quintile$V4))
# [1] 138

# read enhancers
early_bed <- read.table(file = "mouse_enh_meanRT_by_germ_line_cons2sp_EARLY.bed"
                             , header = F, stringsAsFactors = F, sep ='\t')
late_bed <- read.table(file = "mouse_enh_meanRT_by_germ_line_cons2sp_LATE.bed"
                            , header = F, stringsAsFactors = F, sep ='\t')

early_bed$rel_rt <- "early"
late_bed$rel_rt <- "late"

all_enh <- rbind(early_bed, late_bed)

# gr objects
top_quintile.gr <- with(data = top_quintile, expr = GRanges(V1, IRanges(V2 + 1, V3)))
all_enh.gr <- with(data = all_enh, expr = GRanges(V1, IRanges(V2 + 1, V3)))
all_enh$enh_id <- with(all_enh, paste(V1, paste(V2, V3, sep = '-'), sep = ':'))


x <- as.data.frame(findOverlaps(top_quintile.gr, all_enh.gr))
motifs_enh <- cbind(top_quintile[x$queryHits, "V4", drop=F]
                    , all_enh[x$subjectHits, ])
colnames(motifs_enh)[1] <- "ID"
# add enhancer id
motif_freq <- as.data.frame(table(motifs_enh[,c("enh_id", "ID")]))


library(tidyr)
motif_freq <- tidyr::spread(motif_freq, ID, Freq) 
dim(motif_freq[complete.cases(motif_freq),])
# [1] 68729   139

motif_freq <- merge(motif_freq, all_enh[, c("enh_id", "rel_rt")]
                    , by = 'enh_id', all.x=T)
rownames(motif_freq) <- motif_freq$enh_id
motif_freq$enh_id <- NULL

# save
write.table(x = motif_freq, file = "allTissues_homeodomain_topQuintile_humanJAS.txt"
            , quote = F, sep = '\t')
