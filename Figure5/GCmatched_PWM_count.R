library(stringr)

setwd("./replication_timing/human")

cons_early <- read.table(file = "./homer/conserved_early_quintile_H9RT_ext_NULL_humanJAS.tsv"
                         , header = T, stringsAsFactors = F, sep ='\t')
cons_late <- read.table(file = "./homer/conserved_late_quintile_H9RT_ext_NULL_humanJAS.tsv"
                        , header = T, stringsAsFactors = F, sep ='\t')
recent_early <- read.table(file = "./homer/recent_early_quintile_H9RT_ext_NULL_humanJAS.tsv"
                           , header = T, stringsAsFactors = F, sep ='\t')
recent_late <- read.table(file = "./homer/recent_late_quintile_H9RT_ext_NULL_humanJAS.tsv"
                          , header = T, stringsAsFactors = F, sep ='\t')


count_motifs <- function(tsv){
  motif_cols <- grep("Distance.From", colnames(tsv), value=T)
 
  for(motif in motif_cols){
    motif.name <- gsub(pattern = "\\.Distance\\.From\\.Peak\\.sequence\\.strand\\.conservation."
                       , replacement = "", x = motif)
    # print(motif.name)
    counts <- apply(X = (tsv[, motif, drop=F]), MARGIN = 1
                    , FUN = function (x) str_count(string = x, pattern = "\\)"))
    # print(class(counts))
    # adding new column with counts
    tsv[,motif.name] <- counts
  }
  tsv$peakid <- with(tsv, paste(Chr, paste(Start, End, sep = "-"), sep = ":"))
  all_motifs <- gsub(pattern = "\\.Distance\\.From\\.Peak\\.sequence\\.strand\\.conservation."
                     , replacement = "", x = motif_cols)
  motif_counts <- tsv[,c("peakid", all_motifs)]
  motif_counts[is.na(motif_counts)] <- 0
  rownames(motif_counts) <- motif_counts$peakid
  motif_counts$peakid <- NULL
  countSum <- colSums(motif_counts)
  motif_counts <- motif_counts[,names(countSum[countSum > 0])]
 
  return(motif_counts)
}

cons_early.counts <- count_motifs(cons_early)
cons_late.counts <- count_motifs(cons_late)
recent_early.counts <- count_motifs(recent_early)
recent_late.counts <- count_motifs(recent_late)

write.table(x = cons_early.counts, file = "cons_EARLY_NULL_humanJAS2020_homer.txt"
            , quote = F, sep = '\t')
write.table(x = cons_late.counts, file = "cons_LATE_NULL_humanJAS2020_homer.txt"
            , quote = F, sep = '\t')
write.table(x = recent_early.counts, file = "recent_EARLY_NULL_humanJAS2020_homer.txt"
            , quote = F, sep = '\t')
write.table(x = recent_late.counts, file = "recent_LATE_NULL_humanJAS2020_homer.txt"
            , quote = F, sep = '\t')

                    
library(tidyr)                    
early.recent <- read.table(file = "recent_EARLY_NULL_humanJAS2020_homer.txt"
                           , header = T, stringsAsFactors = F)
late.recent <- read.table(file = "recent_LATE_NULL_humanJAS2020_homer.txt"
                          , header = T, stringsAsFactors = F)

early.nonRecent <- read.table(file = "cons_EARLY_NULL_humanJAS2020_homer.txt"
                              , header = T, stringsAsFactors = F)
late.nonRecent <- read.table(file = "cons_LATE_NULL_humanJAS2020_homer.txt"
                             , header = T, stringsAsFactors = F)

early.recent$is.early <- 1
late.recent$is.early <- 0
early.nonRecent$is.early <- 1
late.nonRecent$is.early <- 0

early.recent$is.recent <- 1
late.recent$is.recent <- 1
early.nonRecent$is.recent <- 0
late.nonRecent$is.recent <- 0

recent <- dplyr::bind_rows(early.recent, late.recent)
non_recent <- dplyr::bind_rows(early.nonRecent, late.nonRecent)

recent[is.na(recent)] <- 0 
non_recent[is.na(non_recent)] <- 0 

recent <- recent[,c(colnames(recent)[!colnames(recent) %in% c("is.early", "is.recent")]
                    , "is.early", "is.recent")]
non_recent <- non_recent[,c(colnames(non_recent)[!colnames(non_recent) %in% c("is.early", "is.recent")]
                            , "is.early", "is.recent")]

write.table(x = recent, file = "liver_recent_NULL_earlyvlate_jaspar.txt", quote = F)
write.table(x = non_recent, file = "liver_non.recent_NULL_earlyvlate_jaspar.txt", quote = F)

# table(recent$is.early)
# 0    1
# 1999 1999
# table(non_recent$is.early)
#
# 0    1
# 3539 3539
