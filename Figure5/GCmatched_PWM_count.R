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
