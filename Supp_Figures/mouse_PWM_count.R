library(stringr)

early <- read.table(file = "mouse_recent_and_cons2sp_EARLY_human.jaspar.tsv"
                    , header = T, stringsAsFactors = F, sep ='\t')
late <- read.table(file = "mouse_recent_and_cons2sp_LATE_human.jaspar.tsv"
                   , header = T, stringsAsFactors = F, sep ='\t')

# EARLY
motif_cols <- grep("Distance.From", colnames(early), value=T)

for(motif in motif_cols){
  motif.name <- gsub(pattern = "\\.Distance\\.From\\.Peak\\.sequence\\.strand\\.conservation."
                     , replacement = "", x = motif)
  print(motif.name)
  counts <- apply(X = (early[,motif, drop=F]), MARGIN = 1
                  , FUN = function (x) str_count(string = x, pattern = "\\)"))
  # print(class(counts))
  # adding new column with counts
  early[,motif.name] <- counts
}

# add peak id
early$peakid <- with(early, paste(Chr, paste(Start, End, sep = "-"), sep = ":"))
all_motifs <- gsub(pattern = "\\.Distance\\.From\\.Peak\\.sequence\\.strand\\.conservation."
                   , replacement = "", x = motif_cols)
early_counts <- early[,c("peakid", all_motifs)]
early_counts[is.na(early_counts)] <- 0
rownames(early_counts) <- early_counts$peakid
early_counts$peakid <- NULL
countSum <- colSums(early_counts)
early_counts <- early_counts[,names(countSum[countSum > 0])]
write.table(x = early_counts, file = "mouse_recent_and_cons2sp_EARLY_humanJAS2020_homer.txt"
            , quote = F, sep = '\t')

# LATE
motif_cols <- grep("Distance.From", colnames(late), value=T)#810

for(motif in motif_cols){
  motif.name <- gsub(pattern = "\\.Distance\\.From\\.Peak\\.sequence\\.strand\\.conservation."
                     , replacement = "", x = motif)
  print(motif.name)
  counts <- apply(X = (late[,motif, drop=F]), MARGIN = 1
                  , FUN = function (x) str_count(string = x, pattern = "\\)"))
  # print(class(counts))
  # adding new column with counts
  late[,motif.name] <- counts
}

late$peakid <- with(late, paste(Chr, paste(Start, End, sep = "-"), sep = ":"))
all_motifs <- gsub(pattern = "\\.Distance\\.From\\.Peak\\.sequence\\.strand\\.conservation."
                   , replacement = "", x = motif_cols)
late_counts <- late[,c("peakid", all_motifs)]
late_counts[is.na(late_counts)] <- 0
rownames(late_counts) <- late_counts$peakid
late_counts$peakid <- NULL
countSum <- colSums(late_counts)
late_counts <- late_counts[,names(countSum[countSum > 0])]
write.table(x = late_counts, file = "mouse_recent_and_cons2sp_LATE_humanJAS2020_homer.txt"
            , quote = F, sep = '\t')
