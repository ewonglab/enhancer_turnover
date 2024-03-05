library(data.table)
library(reshape2)
library(GenomicRanges)
library("Hmisc")
library("GenomicFeatures")
library("MatchIt")

### LNCaP replication time

ln_rt <- read.delim('GSE98730_LNCaP_WA.bed', header=F)

names(ln_rt) <- c('chr_rt','start_rt','end_rt','ln_rt')

ln_rt$rel_rt <- ""
ln_rt_gr <- with(ln_rt, GRanges( chr_rt , IRanges( start_rt, end_rt)))

pr <- read.delim('GSE57498_PrEC_ChromHMM.bed'
                 , header=F, skip=1)

pr_enh <- subset(pr, V4=='Enhancer')
pr_enh_gr <- with(pr_enh, GRanges( V1 , IRanges( V2, V3 )))
LNCaP_H3K27ac <- read.delim('GSM1902615_ChIP-seq_LNCaP_H3K27ac_region_distal_1kb.bed'
                            , header = F)

LNCaP_H3K4me1 <- read.delim('GSM1902616_ChIP-seq_LNCaP_H3K4me1_region_distal_1kb.bed'
                            , header = F)

LNCaP_H3K27ac$mark <- 'H3K27ac'
LNCaP_H3K4me1$mark <- 'H3K4me1'
LNCaP_enh <- rbind(LNCaP_H3K27ac, LNCaP_H3K4me1)
LNCaP_enh_gr <- with(LNCaP_enh, GRanges( V1 , IRanges( V2, V3 )))

x <- as.data.frame(findOverlaps( LNCaP_enh_gr , pr_enh_gr))
gains <- LNCaP_enh_gr[ !LNCaP_enh_gr %over% pr_enh_gr, ] 
losses <- pr_enh_gr[ !pr_enh_gr %over% LNCaP_enh_gr, ] 
conserved <- pr_enh_gr[ pr_enh_gr %over% LNCaP_enh_gr, ] 
gains <- unique(gains)

## reading recombination rate from Baca et al
baca <- read.csv(file = "Baca_etal_prostate_breakpoints_tableS3C.csv"
                 , header =T, stringsAsFactors = F)

baca <- baca[,1:23]
baca <- baca[!is.na(baca$Breakpoint.1.chromosome) |
               !is.na(baca$Breakpoint.2.chromosome),]

baca <- baca[with(baca, order(-Quality)),]
baca <- head(baca, 5596)

tmp1 <- baca[,c("Breakpoint.1.chromosome", "Breakpoint.1.position")]
tmp2 <- baca[,c("Breakpoint.2.chromosome", "Breakpoint.2.position")]

colnames(tmp1) <- colnames(tmp2) <- c("Breakpoint.chromosome", "Breakpoint.position")

baca.breakp <- unique(rbind(tmp1, tmp2))
baca.breakp$Breakpoint.chromosome <- paste0("chr", baca.breakp$Breakpoint.chromosome)
baca.breakp$Breakpoint.chromosome <- gsub("chr23", "chrX", baca.breakp$Breakpoint.chromosome)
baca.breakp$Breakpoint.chromosome <- gsub("chr24", "chrY", baca.breakp$Breakpoint.chromosome)

baca.breakp_gr <- with(baca.breakp
                       , GRanges(Breakpoint.chromosome
                                 , IRanges(Breakpoint.position, Breakpoint.position)))#11179


#overlap baca with replication time
ln_rt$bin_id <- paste0("ID", 1:nrow(ln_rt)) 
ln_rt$width <- with(ln_rt, end_rt-start_rt)
ln_rt_gr <- with(ln_rt, GRanges( chr_rt , IRanges( start_rt, end_rt )))


x <- as.data.frame(findOverlaps(gains, ln_rt_gr))
gains_lncap <- cbind(as.data.frame(gains)[(x$queryHits),], ln_rt[x$subjectHits,])
x <- as.data.frame(findOverlaps(losses, ln_rt_gr))
losses_lncap <- cbind(as.data.frame(losses)[(x$queryHits),], ln_rt[x$subjectHits,])
x <- as.data.frame(findOverlaps(conserved, ln_rt_gr))
conserved_lncap <- cbind(as.data.frame(conserved)[(x$queryHits),], ln_rt[x$subjectHits,])


gains_lncap <- gains_lncap[,c("seqnames", "start", "end", "width", "strand", "ln_rt")]
losses_lncap <- losses_lncap[,c("seqnames", "start", "end", "width", "strand", "ln_rt")]
conserved_lncap <- conserved_lncap[,c("seqnames", "start", "end", "width", "strand", "ln_rt")]

# mean RT by enhancer
gains_lncap <- aggregate(ln_rt~., gains_lncap, mean)
losses_lncap <- aggregate(ln_rt~., losses_lncap, mean)
conserved_lncap <- aggregate(ln_rt~., conserved_lncap, mean)

#add id to every region
gains_lncap$enh_id <- with(gains_lncap, paste(seqnames, paste(start, end, sep = '_'), sep = ":"))
losses_lncap$enh_id <- with(losses_lncap, paste(seqnames, paste(start, end, sep = '_'), sep = ":"))
conserved_lncap$enh_id <- with(conserved_lncap, paste(seqnames, paste(start, end, sep = '_'), sep = ":"))


#overlap with breakpoints

add.breakpoints <- function(enh){
  enh_gr <- with(enh, GRanges( seqnames , IRanges( start, end )))
  x <- as.data.frame(findOverlaps(enh_gr, baca.breakp_gr))
  enh_with_baca.bp <- cbind(enh[x$queryHits,], baca.breakp[x$subjectHits,])#
  enh_n.bp_per_enh <- as.data.frame(table(enh_with_baca.bp$enh_id))#30  2
  
  #add rt and enhancer with
  enh_n.bp_per_enh <- merge(enh_n.bp_per_enh
                            , unique(enh[,c("enh_id", "ln_rt", "width")])
                            , by.x="Var1", by.y = "enh_id", all.x=T)#
  
  #Number of breakpoints normalized by enhancer width
  enh_n.bp_per_enh$n_by_width <- with(enh_n.bp_per_enh, Freq/width)
  colnames(enh_n.bp_per_enh) <- gsub("Var1", "enh_id", colnames(enh_n.bp_per_enh))
  return(enh_n.bp_per_enh)}


gains_lncap_baca <- add.breakpoints(gains_lncap)
losses_lncap_baca <- add.breakpoints(losses_lncap)
conserved_lncap_baca <- add.breakpoints(conserved_lncap)

#### USE MATCHINT HERE - SELECT GAINS-UNCHANGED AND LOSSES UNCHANGED MATCHED BY # BREAKPOINTS

#select a set of conserved enhancers at random and use matchit to select gains and losses
match_by_recomb.0 <- function(unch.df, changed.df, n){
  unch.df$study <- "unchanged"
  changed.df$study <- "changed"
  set.seed(123)
  unch_sample <- sample(x = unch.df$enh_id, size = n, replace = F)
  df_sub <- rbind(unch.df[unch.df$enh_id %in% unch_sample,]
                  ,  changed.df)
  df_sub$type_binary <- ifelse(df_sub$study=="changed",0,1)
  m.out1 <- matchit(type_binary ~ n_by_width, data = df_sub,
                    method = "nearest", distance = "glm")
  
  #imbalance after matching:
  print("imbalance after matching:")
  print(summary(m.out1, un = FALSE))
  
  changed_selected <- df_sub[(m.out1$match.matrix)[,1],]
  #add the selected conserved enhancers
  changed_selected$type_binary <- NULL
  
  print("unchanged recombination")
  print(summary(unch.df[unch.df$enh_id %in% unch_sample,"n_by_width"]))
  
  print("changed recombination")
  print(summary(changed_selected$n_by_width))
  
  matched_enh <- rbind(unch.df[unch.df$enh_id %in% unch_sample,]
                       , changed_selected)#
  
  return(matched_enh)}


gains_unch_m.0 <- match_by_recomb.0(conserved_lncap_baca, gains_lncap_baca, n =140)

# [1] "imbalance after matching:"
# 
# Call:
#   matchit(formula = type_binary ~ n_by_width, data = df_sub, method = "nearest", 
#           distance = "glm")
# 
# Summary of Balance for Matched Data:
#   Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
# distance          0.1744        0.1708          0.0831     2.0736    0.0042
# n_by_width        0.0010        0.0009          0.0677     1.6969    0.0042
# eCDF Max Std. Pair Dist.
# distance     0.0429          0.0948
# n_by_width   0.0429          0.0797
# 

# Sample Sizes:
#   Control Treated
# All           686     140
# Matched       140     140
# Unmatched     546       0
# Discarded       0       0

# 
# [1] "unchanged recombination"
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 9.614e-05 3.332e-04 6.246e-04 9.817e-04 1.248e-03 4.975e-03 
# [1] "changed recombination"
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 5.882e-05 3.375e-04 6.246e-04 9.102e-04 1.248e-03 3.630e-03 


losses_unch_m.0 <- match_by_recomb.0(conserved_lncap_baca, losses_lncap_baca, n =140)

# [1] "imbalance after matching:"
# 
# Call:
#   matchit(formula = type_binary ~ n_by_width, data = df_sub, method = "nearest", 
#           distance = "glm")
# 
# Summary of Balance for Matched Data:
#   Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
# distance           0.448        0.4457          0.0385     1.0734    0.0434
# n_by_width         0.001        0.0010         -0.0348     1.0667    0.0434
# eCDF Max Std. Pair Dist.
# distance     0.1429          0.0524
# n_by_width   0.1429          0.0488
# 
# Sample Sizes:
#   Control Treated
# All           180     140
# Matched       140     140
# Unmatched      40       0
# Discarded       0       0
# 

# [1] "unchanged recombination"
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 9.614e-05 3.332e-04 6.246e-04 9.817e-04 1.248e-03 4.975e-03 
# [1] "changed recombination"
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0001087 0.0004165 0.0006246 0.0010184 0.0012484 0.0049751 


#   MATCHED ENHANCERS

x <- wilcox.test(gains_unch_m.0[gains_unch_m.0$study == "unchanged", "ln_rt"]
                 , gains_unch_m.0[gains_unch_m.0$study == "changed", "ln_rt"]
                 , alternative = "greater")#

# Wilcoxon rank sum test with continuity correction
# 
# data:  gains_unch_m.0[gains_unch_m.0$study == "unchanged", "ln_rt"] and gains_unch_m.0[gains_unch_m.0$study == "changed", "ln_rt"]
# W = 12068, p-value = 0.0004084
# alternative hypothesis: true location shift is greater than 0
x$p.value # 0.0004084439

x <- wilcox.test(losses_unch_m.0[losses_unch_m.0$study == "unchanged", "ln_rt"]
                 , losses_unch_m.0[losses_unch_m.0$study == "changed", "ln_rt"]
                 , alternative = "greater")

# Wilcoxon rank sum test with continuity correction
# 

# data:  losses_unch_m.0[losses_unch_m.0$study == "unchanged", "ln_rt"] and losses_unch_m.0[losses_unch_m.0$study == "changed", "ln_rt"]
# W = 11031, p-value = 0.03466
# alternative hypothesis: true location shift is greater than 0
x$p.value #0.03466107


#plots
losses_unch_m.0[losses_unch_m.0$study=="changed", "study"]<- "losses"
gains_unch_m.0[gains_unch_m.0$study=="changed", "study"]<- "gains"

df <- rbind(losses_unch_m.0, gains_unch_m.0)

#unchanged enhancers are duplicated

df <- unique(df)
df$study <- factor(df$study, levels=c("gains", "unchanged", "losses"))
table( df$study)

# gains unchanged    losses 
# 140       140       140 


#number of breakpoints / enhancer width of selected enhancer
pdf("all_prostate_enh_matched_by_baca_bk_by_enh_width_140.pdf")
ggplot(df, aes(x=study, y=log10(n_by_width))) + 
  geom_boxplot() + theme_classic()
dev.off()


#plotting replication time of selected enhancers
pdf("all_prostate_enh_matched_by_baca_bk_LNCaP_rt_140.pdf")
ggplot(df, aes(x=study, y=log10(ln_rt))) + 
  geom_boxplot() + theme_classic()
dev.off()

# difference in breakpoints/enhancer width
x <- wilcox.test(gains_unch_m.0[gains_unch_m.0$study == "unchanged", "n_by_width"]
                 , gains_unch_m.0[gains_unch_m.0$study == "gains", "n_by_width"])

x$p.value
# [1] 0.9493489


x <- wilcox.test(losses_unch_m.0[losses_unch_m.0$study == "unchanged", "n_by_width"]
                 , losses_unch_m.0[losses_unch_m.0$study == "losses", "n_by_width"])

x$p.value
#[1] 0.2392059
