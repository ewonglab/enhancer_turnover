library(data.table)
library(GenomicRanges)
library(ggplot2)

setwd("./replication_timing/human")

#read enhancers
enhancers <- fread("/rt/liver/Hsap_Enhancers_conservationInOthers")
recent_e <- fread("/rt/liver/Hsap_H3K27ac_humanspEnhancers")

enhancers$peakid <- with(enhancers, paste(Chrom, paste(Start, End, sep = '_'), sep = ":"))
recent_e$peakid <- with(recent_e, paste(Chrom, paste(Start, End, sep = '_'), sep = ":"))

enhancers <- as.data.frame(enhancers) 
recent_e <- as.data.frame(recent_e) 

#new definition of conserved enhancers and promoters
species <- c("Fcat", "Csab", "Mdom", "Ddel", "Mbid", "Cpor"
             , "Bbor", "Rnor", "Btau", "Cjac", "Cfam", "Mmus"
             , "Sscr", "Ocun", "Tbel", "Shar", "Mmul", "Mfur")

#separate non-recent
cons_e <- enhancers[!enhancers$peakid %in% recent_e$peakid,] 

#add the number of species where the enhancer is conserved
cons_e$n_species <- apply(cons_e[,species],1, function(x) length(which(x!="-" & !is.na(x))))

# keep only enhancers and promoters conserved in any other two species
cons_enh_2species <- cons_e[cons_e$n_species >= 2, ]
cons_enh_2species$n_species <- NULL
cons_enh_2species$type <- "conserved"
recent_e$type <- "recent"

# conserved in 2 species + recent
enh <- rbind(cons_enh_2species, recent_e)

# H9 DNA replication time                          
h9_rt <- read.delim('RT_H9_ESC_Ext29405702_hg19.bedgraph'
                    , skip=11, header=F)
h9_rt.gr <- with(h9_rt, GRanges(V1, IRanges(V2+1, V3)))#
colnames(h9_rt)[4] <- "h9_rt"

## human ENHANCERS THAT ALIGN TO AT LEAST ONE SPECIES
enh$Chrom <- paste0("chr", enh$Chrom)
enh_gr <- with(enh, GRanges(Chrom, IRanges(Start+1, End)))

x <- as.data.frame(findOverlaps(enh_gr, h9_rt.gr))
enh_gr_with_h9.rt <- cbind(enh[x$queryHits, c("Chrom", "Start", "End", "type")]
                           , h9_rt[x$subjectHits, "h9_rt", drop =F])
enh_gr_with_h9.rt <- aggregate(h9_rt ~., enh_gr_with_h9.rt, mean)


### define relative RT
enh_gr_with_h9.rt$rel_rt <- ""
enh_gr_with_h9.rt[enh_gr_with_h9.rt$h9_rt > 0.5,"rel_rt"] <- "early"
enh_gr_with_h9.rt[enh_gr_with_h9.rt$h9_rt < -0.5,"rel_rt"] <- "late"
enh_gr_with_h9.rt <- enh_gr_with_h9.rt[enh_gr_with_h9.rt$rel_rt != "",]

# make new gr object
enh_gr_with_h9.rt$Chrom <- sub("chr", "", enh_gr_with_h9.rt$Chrom)
enh_gr_with_h9.rt_gr <- with(enh_gr_with_h9.rt, GRanges(Chrom, IRanges(Start+1, End)))

# reading MPRA data

starr_seq <- read.csv(file = "all_tiles_scores_Klein.csv", header = T, stringsAsFactors = F)
starr_seq <- starr_seq[,1:2] # 6859    2

# separate ids into tile coordinates
starr_seq$chr <- sub(":.*", "", starr_seq$Enhancer.Coordinates..hg19.)
starr_seq$start <- sub("-.*", "", sub(".*:", "", starr_seq$Enhancer.Coordinates..hg19.))
starr_seq$end <- sub(".*-", "", sub(".*:", "", starr_seq$Enhancer.Coordinates..hg19.))

# remove "Negative" controls
starr_seq <- starr_seq[starr_seq$Enhancer.Coordinates..hg19. != "Negative",] # 6735    5

starr_seq$start <- as.numeric(starr_seq$start)
starr_seq$end <- as.numeric(starr_seq$end)

# make gr object
starr_seq_gr <- with(starr_seq, GRanges(chr , IRanges(start + 1, end))) # 6735

x <- as.data.frame(findOverlaps(enh_gr_with_h9.rt_gr, starr_seq_gr))

enh_activity <- cbind(enh_gr_with_h9.rt[x$queryHits, ]
                      , starr_seq[x$subjectHits, "log2.Starr.seq.score.", drop =F])
enh_activity <- aggregate(log2.Starr.seq.score. ~., enh_activity, mean) #  524   7


dim(unique(enh_activity[,c("Chrom", "Start", "End")]))
# [1] 524   3
table(enh_activity$rel_rt)
# early  late
# 455 69

# Number of conserved and recent enhancers
table(enh_activity$type)
# conserved    recent
# 270       254


# separate between recent and conserved enhancers

pdf("villar_STARRseq_score_byRTclass_byType.pdf")
ggplot(enh_activity, aes(x = rel_rt, y = log2.Starr.seq.score.)) +
  geom_boxplot() + theme_classic() + facet_wrap(.~type)
dev.off()

# conserved enhancers
x <- wilcox.test(enh_activity[enh_activity$rel_rt == "early" &
                                enh_activity$type == "conserved"
                              , "log2.Starr.seq.score."]
                 , enh_activity[enh_activity$rel_rt == "late" &
                                  enh_activity$type == "conserved"
                                , "log2.Starr.seq.score."])

# Wilcoxon rank sum test with continuity correction
#
# data:  enh_activity[enh_activity$rel_rt == "early" & enh_activity$type == "conserved", "log2.Starr.seq.score."] and enh_activity[enh_activity$rel_rt == "late" & enh_activity$type == "conserved", "log2.Starr.seq.score."]
# W = 2562, p-value = 0.1873
# alternative hypothesis: true location shift is not equal to 0

# recent enhancers
x <- wilcox.test(enh_activity[enh_activity$rel_rt == "early" &
                                enh_activity$type == "recent"
                              , "log2.Starr.seq.score."]
                 , enh_activity[enh_activity$rel_rt == "late" &
                                  enh_activity$type == "recent"
                                , "log2.Starr.seq.score."])

# Wilcoxon rank sum test with continuity correction
#
# data:  enh_activity[enh_activity$rel_rt == "early" & enh_activity$type == "recent", "log2.Starr.seq.score."] and enh_activity[enh_activity$rel_rt == "late" & enh_activity$type == "recent", "log2.Starr.seq.score."]
# W = 6091, p-value = 0.07593
# alternative hypothesis: true location shift is not equal to 0
