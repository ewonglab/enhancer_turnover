# Forest plot of conditional probabilities: P(recent | late RT) / P (recent | early RT)

library(data.table)
library(GenomicRanges)
require(dplyr)
library(ggplot2)
library("Hmisc")

mouse.enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                        , header = F, stringsAsFactors = F
                        , sep = "\t")
active_atLeast2sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast2sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t") 
cons_enh <- rbind(active_atLeast2sp, poised_atLeast2sp)
cons_enh$V1 <- paste0("chr", cons_enh$V1)
unique(sub(":.*", "", cons_enh$V1))
# [1] "chr1"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17"
# [10] "chr18" "chr19" "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"
# [19] "chr9"  "chrX"

mouse.enh <- mouse.enh[mouse.enh$V6 == "Recent" | mouse.enh$V4 %in% cons_enh$V1, ] 

rt <- read.table(file = "mouse_sperm_pgc_meanRT.bedGraph"
                 , header = F, stringsAsFactors = F, sep = '\t')
rt <- rt[complete.cases(rt),]

mouse_enh.gr <- with(mouse.enh, GRanges(V1, IRanges(V2, V3)))
rt.gr <- with(rt, GRanges(V1, IRanges(V2, V3)))

#findOverlaps
x <- as.data.frame(findOverlaps(mouse_enh.gr, rt.gr))
colnames(rt)[4] <- "rt"
mouse_enh <- cbind(mouse.enh[x$queryHits,]
                   , rt[x$subjectHits, "rt", drop =F])
length(unique(mouse_enh$V4))

mouse_enh <- aggregate(rt ~., mouse_enh, mean)
#defining early and late enhancers
mouse_enh$rel_rt <- ""
mouse_enh[mouse_enh$rt > 0.5,"rel_rt"] <- "early"
mouse_enh[mouse_enh$rt < -0.5,"rel_rt"] <- "late"

mouse_enh <- mouse_enh[!mouse_enh$rel_rt == "",]



condit.prob_balanced <- function(x, n){
  early_sample <- sample(x =  x[ x$rel_rt == "early","V4"]
                         , size = n, replace = F)
  late_sample <- sample(x =  x[ x$rel_rt == "late","V4"]
                        , size = n, replace = F)
  x <-  x[ x$V4 %in% c(early_sample, late_sample),]
  df <- as.data.frame(table( x[, c("V6", "rel_rt")]))
  colnames(df) <- c("type", "rel_rt", "Freq")
  total_enh <- sum(df$Freq)
  p_early <- sum(df[df$rel_rt == "early", "Freq"])/total_enh
  p_late <- sum(df[df$rel_rt == "late", "Freq"])/total_enh
  p_recent_and_early <- sum(df[df$type == "Recent" &
                                 df$rel_rt == "early", "Freq"])/total_enh
  p_recent_and_late <- sum(df[df$type == "Recent" &
                                df$rel_rt == "late", "Freq"])/total_enh
  p_recent_given_early <- p_recent_and_early/p_early
  p_recent_given_late <- p_recent_and_late/p_late
  cond.p <- p_recent_given_late/p_recent_given_early
  return(cond.p)
}


# only tissue specific enhancers
as.data.frame(table(mouse_enh[mouse_enh$V7 %in% c("brain", "liver", "muscle", "testis")
                              ,c("rel_rt", "V5", "V7")]))

# rel_rt     V5     V7  Freq
# 1   early Active  brain  8248
# 2    late Active  brain   913
# 3   early Poised  brain  9515
# 4    late Poised  brain  2130
# 5   early Active  liver  6140
# 6    late Active  liver   739
# 7   early Poised  liver 10364
# 8    late Poised  liver  2382
# 9   early Active muscle  7557
# 10   late Active muscle   470
# 11  early Poised muscle  9769
# 12   late Poised muscle  1030
# 13  early Active testis  4070
# 14   late Active testis   436
# 15  early Poised testis 12313
# 16   late Poised testis  2112


library(boot)
#### TISSUE SPECIFIC ENHANCERS   ## ACTIVE ENHANCERS

liver.specific.act <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="liver" & mouse_enh$V5=="Active",], n = 400)
  liver.specific.act <- c(liver.specific.act, tmp)}

brain.specific.act <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="brain" & mouse_enh$V5=="Active",], n = 400)
  brain.specific.act <- c(brain.specific.act, tmp)}

muscle.specific.act <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="muscle" & mouse_enh$V5=="Active",], n = 400)
  muscle.specific.act <- c(muscle.specific.act, tmp)}

testis.specific.act <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="testis" & mouse_enh$V5=="Active",], n = 400)
  testis.specific.act <- c(testis.specific.act, tmp)}

#### TISSUE SPECIFIC ENHANCERS   ## POISED ENHANCERS
liver.specific.pois <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="liver" & mouse_enh$V5=="Poised",], n = 400)
  liver.specific.pois <- c(liver.specific.pois, tmp)}

brain.specific.pois <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="brain" & mouse_enh$V5=="Poised",], n = 400)
  brain.specific.pois <- c(brain.specific.pois, tmp)}

muscle.specific.pois <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="muscle" & mouse_enh$V5=="Poised",], n = 400)
  muscle.specific.pois <- c(muscle.specific.pois, tmp)}

testis.specific.pois <- c()
set.seed(42)
for(i in 1:1000){
  tmp <- condit.prob_balanced(x = mouse_enh[mouse_enh$V7=="testis" & mouse_enh$V5=="Poised",], n = 400)
  testis.specific.pois <- c(testis.specific.pois, tmp)}

####

library(reshape2)
library(ggplot2)


#log version
all.p <- data.frame(sample=c("brain.specific.act", "liver.specific.act"
                             , "muscle.specific.act", "testis.specific.act"
                             , "brain.specific.pois", "liver.specific.pois"
                             , "muscle.specific.pois", "testis.specific.pois")
                    , mean.p=c(mean(log(brain.specific.act)), mean(log(liver.specific.act))
                               , mean(log(muscle.specific.act)), mean(log(testis.specific.act))
                               , mean(log(brain.specific.pois)), mean(log(liver.specific.pois))
                               , mean(log(muscle.specific.pois)), mean(log(testis.specific.pois)))
                    , se.p=c(sd(log(brain.specific.act)), sd(log(liver.specific.act))
                             , sd(log(muscle.specific.act)), sd(log(testis.specific.act))
                             , sd(log(brain.specific.pois)), sd(log(liver.specific.pois))
                             , sd(log(muscle.specific.pois)), sd(log(testis.specific.pois)))
                    , mark = c(rep("Active", 4), rep("Poised", 4))
                    , stringsAsFactors = F)

all.p$lower <- all.p$mean - all.p$se.p
all.p$upper <- all.p$mean + all.p$se.p

all.p$sample <- gsub("act$", "", all.p$sample)
all.p$sample <- gsub("pois$", "", all.p$sample)
all.p$sample <- factor(all.p$sample
                       , levels = rev(c("brain.specific.", "liver.specific."
                                        , "muscle.specific.", "testis.specific.")))
all.p$mark <- factor(all.p$mark, levels = c("Poised", "Active"))

pdf("conditional_recent_earlyvlate_forest_onlyMouse_cons2sp_log.pdf")
ggplot(data=all.p, aes(x=sample, y=mean.p, ymin=lower, ymax=upper)) +
  geom_pointrange() +
  # geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  # xlab("Label") + ylab("Mean (95% CI)") +
  theme_classic()  + facet_wrap(mark~.) +
  scale_y_continuous(limits=c(-0.5,0.8), breaks=seq(-0.5,0.8, by = 0.1))
dev.off()
