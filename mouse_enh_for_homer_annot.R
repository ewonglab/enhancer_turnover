setwd("./replication_timing/mouse")
all_enh <- read.table(file = "./mouse_enh_meanRT_by_germ_line_cons2sp.txt"
                      , header = T, stringsAsFactors = F, sep ='\t')#  134195     11
all_enh$chr <- sub("chr", "", all_enh$chr)
# separate between early RT and late RT
rt_cols <- c("PGC.female.1.bedGraph", "PGC.male.1.bedGraph"
             , "Sperm.1.bedGraph", "Sperm.2.bedGraph")
all_enh$mean_rt <- apply(all_enh[,rt_cols], 1, mean)
#add relative replication time
all_enh$rel_rt <- ""
all_enh$rel_rt <- ifelse(all_enh$mean_rt > 0.5
                         , "early", all_enh$rel_rt)

all_enh$rel_rt <- ifelse(all_enh$mean_rt < -0.5
                         , "late", all_enh$rel_rt)

#keeping only early and late enhancers
all_enh <- all_enh[all_enh$rel_rt != "",]

write.table(x = all_enh[all_enh$rel_rt == "early",c("chr", "start", "end")]
            , file = "mouse_enh_meanRT_by_germ_line_cons2sp_EARLY.bed"
            , col.names = F, row.names = F, sep = '\t', quote = F)
write.table(x = all_enh[all_enh$rel_rt == "late",c("chr", "start", "end")]
            , file = "mouse_enh_meanRT_by_germ_line_cons2sp_LATE.bed"
            , col.names = F, row.names = F, sep = '\t', quote = F)
