library(GenomicRanges)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)

setwd("./replication_timing/mouse")
# ENHANCERS
mouse.enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                        , header = F, stringsAsFactors = F
                        , sep = "\t")
recent_enh <- mouse.enh[mouse.enh$V6 == "Recent",]
active_atLeast2sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast2sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
cons_enh <- rbind(active_atLeast2sp, poised_atLeast2sp)

# SPLIT COORDINATES OF CONSERVED ENHANCERS
cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end)

all.recent <- recent_enh 
all.cons <- cons_enh 

mouse_only <- read.table(file = "mouse_repeatmasker_mouse_only"
                           , header = T, stringsAsFactors = F, sep ='\t')
ancestral <- read.table(file = "mouse_repeatmasker_ancestral"
                          , header = T, stringsAsFactors = F, sep ='\t')

mouse.only_rep_classes <- unique(mouse_only$repeat_class.family)
anc_rep_classes <- unique(ancestral$repeat_class.family)
anc_rep_classes <- anc_rep_classes[grep("LINE", anc_rep_classes, invert = T)]
anc_rep_classes <- anc_rep_classes[grep("SINE", anc_rep_classes, invert = T)]
anc_rep_classes <- anc_rep_classes[grep("LTR", anc_rep_classes, invert = T)]


ancestral <- ancestral[!ancestral$repeat_class.family %in% anc_rep_classes,]
unique(ancestral$repeat_class.family)
# [1] "SINE/B2"       "SINE/B4"       "LTR/ERVL"      "SINE/ID"      
# [5] "LTR/ERVL-MaLR" "LTR/ERVK"      "SINE/Alu"      "LINE/L1"      
# [9] "LTR/Gypsy"     "SINE/MIR"      "LINE/L2"       "LTR/ERV1"    
# [13] "LINE/CR1"      "LTR/Gypsy?"    "LTR?"          "LTR"          
# [17] "SINE/tRNA-RTE" "LINE/Jockey"   "LINE/RTE-X"    "LTR/ERVL?"    
# [21] "LINE/L1-Tx1"


overlap_enh_repeats <- function(enhancers.tab, repeats.tab){
  # colnames(enhancers.tab) <- c("chr", "start", "end")
  enhancers.gr <- with(enhancers.tab, GRanges(chr, IRanges(start + 1, end)))
  repeats.gr <- with(repeats.tab, GRanges(query_sequence
                                          , IRanges(query_start + 1, query_end)))
  enh_ov_repeats <- as.data.frame(findOverlaps(enhancers.gr
                                               , repeats.gr))
  print(dim(enh_ov_repeats))
  out <- cbind(enhancers.tab[enh_ov_repeats$queryHits,]
               , repeats.tab[enh_ov_repeats$subjectHits,])
  out$enhancer.id <- with(out, paste(chr, paste(start, end, sep = "_"), sep = ":"))
  return(out)}

colnames(all.recent)[1:3] <- c("chr", "start", "end")

recent_ov_repeats <- overlap_enh_repeats(all.recent, rbind(ancestral, mouse_only))
cons_ov_repeats <- overlap_enh_repeats(all.cons, rbind(ancestral, mouse_only)) 


df <- data.frame(recent_ov_rep=length(unique(recent_ov_repeats$enhancer.id))
                 , recent_NOT_ov_rep = nrow(unique(all.recent))-
                   length(unique(recent_ov_repeats$enhancer.id)))
df <- melt(df)
# variable value
# 1     recent_ov_rep 39645
# 2 recent_NOT_ov_rep 41259

df <- data.frame(cons_ov_rep=length(unique(cons_ov_repeats$enhancer.id))
                 , cons_NOT_ov_rep = nrow(unique(all.cons))-
                   length(unique(cons_ov_repeats$enhancer.id)))
df <- melt(df)
#      variable value
# 1     cons_ov_rep 45414
# 2 cons_NOT_ov_rep 42323

### proportion of enhancers overlapping old and recent repeats families
cons_recent.rep <- overlap_enh_repeats(all.cons, mouse_only)
cons_ancestral.rep <- overlap_enh_repeats(all.cons, ancestral)

cons_ov_recent_rep <- (unique(cons_recent.rep$enhancer.id))
cons_ov_anc_rep <- (unique(cons_ancestral.rep$enhancer.id))

#3 classes:
#overlaps ONLY recent repeats
#overlaps ONLY conserved repeats
#overlaps both types of repeats
cons_ov_any_rep <- intersect(cons_ov_recent_rep, cons_ov_anc_rep)
cons_ov_recent_rep <- setdiff(cons_ov_recent_rep, cons_ov_any_rep)
cons_ov_anc_rep <- setdiff(cons_ov_anc_rep, cons_ov_any_rep)

df <- data.frame(cons_ov_recent= length(cons_ov_recent_rep)
                 , cons_ov_anc = length(cons_ov_anc_rep)
                 , cons_ov_both = length(cons_ov_any_rep))
df <- melt(df)
#         variable value
# 1 cons_ov_recent  2007
# 2    cons_ov_anc 36458
# 3   cons_ov_both  6949

#removing enhancers overlapping both types of repeats
df <- df[-3,]
#         variable value
# 1 cons_ov_recent  2007
# 2    cons_ov_anc 36458

round((2007*100)/(2007 + 36458), 2)# 5.22
round((36458*100)/(2007 + 36458), 2)# 94.78

pdf("mouse_cons2sp_enh_ov_retro_byRetroClass.pdf")
ggplot(df, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()
dev.off()

#there are enhancers that overlap with both times of repeats
recent_recent.rep <- overlap_enh_repeats(all.recent, mouse_only)
recent_ancestral.rep <- overlap_enh_repeats(all.recent, ancestral)

recent_ov_recent_rep <- (unique(recent_recent.rep$enhancer.id))
recent_ov_anc_rep <- (unique(recent_ancestral.rep$enhancer.id))

#3 classes:
#overlaps ONLY recent repeats
#overlaps ONLY conserved repeats
#overlaps both types of repeats

recent_ov_any_rep <- intersect(recent_ov_recent_rep, recent_ov_anc_rep)
recent_ov_recent_rep <- setdiff(recent_ov_recent_rep, recent_ov_any_rep)
recent_ov_anc_rep <- setdiff(recent_ov_anc_rep, recent_ov_any_rep)

df <- data.frame(recent_ov_recent= length(recent_ov_recent_rep)
                 , recent_ov_anc = length(recent_ov_anc_rep)
                 , recent_ov_both = length(recent_ov_any_rep))
df <- melt(df)
#           variable value
# 1 recent_ov_recent  4897
# 2    recent_ov_anc 29724
# 3   recent_ov_both  5024

#removing enhancers overlapping both types of repeats
df <- df[-3,]
#           variable value
# 1 recent_ov_recent  4897
# 2    recent_ov_anc 29724

round((4897*100)/(4897 + 29724), 2)# 14.14
round((29724*100)/(4897 + 29724), 2)# 85.86

pdf("mouse_recent_enh_ov_retro_byRetroClass.pdf") 
ggplot(df, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()
dev.off()

m <- matrix(data = c(4897, 2007, 29724, 36458), nrow = 2, ncol = 2, byrow = T)
colnames(m) <- c("recent.enh", "cons.enh")
rownames(m) <- c("recent.repeats", "anc.repeats")

x <- fisher.test(m)
x
# Fisher's Exact Test for Count Data
#
# data:  m
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.834223 3.160980
# sample estimates:
# odds ratio
#   2.992695
#
x$p.value
# [1] 0
x$estimate
# odds ratio
# 2.992695
