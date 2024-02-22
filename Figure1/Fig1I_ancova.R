library(GenomicRanges)

cons.recent.ov.cluster <- function(clust_bed, query_clust, enh_tab){
  clust <- subset(clust_bed, kmeans.cluster==query_clust)
  clust_gr <- with(clust, GRanges( chr , IRanges( start+1, end )))
  enh_gr <- with(enh_tab, GRanges( chr , IRanges( start+1, end )))
  x <-  as.data.frame(findOverlaps(enh_gr, clust_gr))
  result <- enh_tab[unique(x$queryHits),]
  return(table(result$enh_type))}

setwd("./replication_timing/mouse")
mouse_RT <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header =T
                       , stringsAsFactors = F, sep = '\t')#8966   27

somatic <- c("X46C", "D3", "TT2", "iPSC", "iPSC_1D4"
             , "iPSC_2D4", "EPL", "EBM3", "EpiSC5", "EpiSC7"
             , "EBM6", "X46CNPC", "TT2NPC", "EBM9", "Mesoderm"
             , "Endoderm", "piPSC_1A2", "piPSC_1B3", "piPSC_V3"
             , "MEF_female", "MEF_male", "Myoblast")
germline <- c("sperm.1", "sperm.2", "PGC.male.1", "PGC.female.1")


mouse_RT$mean_germlineRT <- apply(mouse_RT[,germline], 1, mean)
mouse_RT$mean_somaticRT <- apply(mouse_RT[,somatic], 1, mean)

mean_mouseRT <- aggregate(.~kmeans.cluster
                          , mouse_RT[,c("mean_germlineRT", "mean_somaticRT","kmeans.cluster")]
                          , mean)

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

# separate tissues
library(tidyr)
library(dplyr)

recent_enh <- recent_enh %>%
  mutate(V7 = strsplit(as.character(V7), ",")) %>%
  unnest(V7)

cons_enh <- cons_enh %>%
  mutate(tissue = strsplit(as.character(tissue), ",")) %>%
  unnest(tissue)

recent_enh <- as.data.frame(recent_enh)
cons_enh <- as.data.frame(cons_enh)

cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end)

# bins and cluster membership
all <- mouse_RT[,"kmeans.cluster",drop =F]
all$peakid <- rownames(all)
all$chr <- gsub(":.*", "", all$peakid)
all$start <- gsub("_.*", "", gsub(".*:", "", all$peakid))
all$end <- gsub(".*_", "", gsub(".*:", "", all$peakid))
all$start <- as.integer(all$start)
all$end <- as.integer(all$end)
all$width <- with(all, end-start)

# poised  - all tissues
# keep poised
poised_recent.long <- recent_enh[recent_enh$V5 == "Poised",]
poised_cons <- cons_enh[cons_enh$type == "poised",]
poised_recent.long$enh_type <- "recent"
poised_cons$enh_type <- "conserved"

poised_recent.long <- unique(poised_recent.long)
poised_cons <- unique(poised_cons)

tmp <- data.frame(chr = c(poised_recent.long$V1, poised_cons$chr)
                  , start = c(poised_recent.long$V2, poised_cons$start)
                  , end = c(poised_recent.long$V3, poised_cons$end)
                  , tissue = c(tolower(poised_recent.long$V7), poised_cons$tissue)
                  , enh_type = c(poised_recent.long$enh_type, poised_cons$enh_type)
                  , stringsAsFactors = F)
tmp$chr <- sub("chrG", "G", tmp$chr)

#separate into tissue pairs
tmp$tissue_pair <- ""
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("liver", "testis"), "liver_testis", tmp$tissue_pair)
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("brain", "muscle"), "brain_muscle", tmp$tissue_pair)

tmp$tissue <- NULL
tmp <- unique(tmp)
table(tmp$tissue_pair)
# brain_muscle liver_testis
# 54518        60171

mean_mouseRT$kmeans.cluster
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue_pair=="liver_testis",]))
mean_mouseRT$liver_testis.conserved.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$liver_testis.recent.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue_pair=="brain_muscle",]))
mean_mouseRT$brain_muscle.conserved.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$brain_muscle.recent.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

### active all tissues
active_recent.long <- recent_enh[recent_enh$V5 == "Active",]
active_cons <- cons_enh[cons_enh$type == "active",]

active_recent.long$enh_type <- "recent"
active_cons$enh_type <- "conserved"
active_recent.long <- unique(active_recent.long)
active_cons <- unique(active_cons)
tmp <- data.frame(chr = c(active_recent.long$V1, active_cons$chr)
                  , start = c(active_recent.long$V2, active_cons$start)
                  , end = c(active_recent.long$V3, active_cons$end)
                  , tissue = c(tolower(active_recent.long$V7), active_cons$tissue)
                  , enh_type = c(active_recent.long$enh_type, active_cons$enh_type)
                  , stringsAsFactors = F)# 67772     5

#separate into tissue pairs
tmp$tissue_pair <- ""
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("liver", "testis"), "liver_testis", tmp$tissue_pair)
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("brain", "muscle"), "brain_muscle", tmp$tissue_pair)

tmp$tissue <- NULL
tmp <- unique(tmp)
table(tmp$tissue_pair)
# brain_muscle liver_testis
# 44645        30965

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue_pair=="liver_testis",]))
mean_mouseRT$liver_testis.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$liver_testis.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue_pair=="brain_muscle",]))
mean_mouseRT$brain_muscle.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$brain_muscle.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))


### RECENT/CONSERVED  log FC TISSUES

mean_mouseRT$liver_testis.poised_rec.cons_logFC <- with(mean_mouseRT
                                                        , log(liver_testis.recent.poised/liver_testis.conserved.poised))
mean_mouseRT$liver_testis.active_rec.cons_logFC <- with(mean_mouseRT
                                                        , log(liver_testis.recent.active/liver_testis.conserved.active))
mean_mouseRT$brain_muscle.poised_rec.cons_logFC <- with(mean_mouseRT
                                                        , log(brain_muscle.recent.poised/brain_muscle.conserved.poised))
mean_mouseRT$brain_muscle.active_rec.cons_logFC <- with(mean_mouseRT
                                                        , log(brain_muscle.recent.active/brain_muscle.conserved.active))

write.table(x = mean_mouseRT, file = "mouse_mean_rt_logFC_by_tissuePairs_mm10_consAlLeast2sp"
            , quote = F, sep = "\t")


library(rstatix)
library(reshape2)
mean.RT.2.clust <- read.table(file = "mouse_mean_rt_logFC_by_tissuePairs_mm10_consAlLeast2sp"
                              , stringsAsFactors = F, sep = "\t", header =T)#18 15

tmp <- mean.RT.2.clust[,c(1:3, 12:15)]

tmp <- melt(tmp, id = c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT"))
tmp$variable <- factor(tmp$variable, levels=c("liver_testis.poised_rec.cons_logFC"
                                              , "liver_testis.active_rec.cons_logFC"
                                              , "brain_muscle.poised_rec.cons_logFC"
                                              , "brain_muscle.active_rec.cons_logFC"))

tmp$enh_type <- ifelse(grepl("active", tmp$variable), "Active", "Poised")
colnames(tmp) <- c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT", "group"
                   , "logFC_recent_cons", "enh_type")

#melt replication time
tmp2 <- melt(tmp, id=c("kmeans.cluster", "group", "logFC_recent_cons", "enh_type"))
tmp2$tissue_pair <- ifelse(grepl("liver_testis", tmp2$group), "liver_testis", "brain_muscle")

# germline poised
tmp2[tmp2$variable=="mean_germlineRT" & tmp2$enh_type=="Poised",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)
# ANOVA Table (type II tests)
#
# Effect DFn DFd       F        p p<.05   ges
# 1       tissue_pair   1  32 150.922 1.18e-13     * 0.825
# 2             value   1  32 502.929 3.87e-21     * 0.940
# 3 tissue_pair:value   1  32  45.811 1.20e-07     * 0.589

#germline active
tmp2[tmp2$variable=="mean_germlineRT" & tmp2$enh_type=="Active",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd       F        p p<.05   ges
# 1       tissue_pair   1  32  72.726 9.58e-10     * 0.694
# 2             value   1  32 280.419 2.16e-17     * 0.898
# 3 tissue_pair:value   1  32  23.521 3.07e-05     * 0.424

#germline all
tmp2[tmp2$variable=="mean_germlineRT",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd       F        p p<.05   ges
# 1       tissue_pair   1  68  62.895 2.93e-11     * 0.481
# 2             value   1  68 224.699 3.07e-23     * 0.768
# 3 tissue_pair:value   1  68  19.671 3.45e-05     * 0.224

#somatic poised
tmp2[tmp2$variable=="mean_somaticRT" & tmp2$enh_type=="Poised",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd      F        p p<.05   ges
# 1       tissue_pair   1  32 24.630 2.21e-05     * 0.435
# 2             value   1  32 57.685 1.19e-08     * 0.643
# 3 tissue_pair:value   1  32  5.090 3.10e-02     * 0.137
 
#somatic active
tmp2[tmp2$variable=="mean_somaticRT" & tmp2$enh_type=="Active",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd      F        p p<.05   ges
# 1       tissue_pair   1  32 15.844 3.70e-04     * 0.331
# 2             value   1  32 37.738 7.18e-07     * 0.541
# 3 tissue_pair:value   1  32  3.450 7.20e-02       0.097

#somatic all
tmp2[tmp2$variable=="mean_somaticRT",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)
# ANOVA Table (type II tests)
#
# Effect DFn DFd      F        p p<.05   ges
# 1       tissue_pair   1  68 28.404 1.21e-06     * 0.295
# 2             value   1  68 67.050 9.97e-12     * 0.496
# 3 tissue_pair:value   1  68  6.017 1.70e-02     * 0.081


# Using linear regression

mean.RT.2.clust <- read.table(file = "mouse_mean_rt_logFC_by_tissuePairs_mm10_consAlLeast2sp"
                              , stringsAsFactors = F, sep = "\t", header =T)#18 15

tmp <- mean.RT.2.clust[,c(1:3, 12:15)]

tmp <- melt(tmp, id = c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT"))
tmp$variable <- factor(tmp$variable, levels=c("liver_testis.poised_rec.cons_logFC"
                                              , "liver_testis.active_rec.cons_logFC"
                                              , "brain_muscle.poised_rec.cons_logFC"
                                              , "brain_muscle.active_rec.cons_logFC"))

tmp$enh_type <- ifelse(grepl("active", tmp$variable), "Active", "Poised")
colnames(tmp) <- c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT", "group"
                   , "logFC_recent_cons", "enh_type")

tmp$tissue_pair <- ifelse(grepl("liver_testis", tmp$group)
                          , "liver_testis", "brain_muscle")
tmp$tissue_pair_bin <- ifelse(tmp$tissue_pair=="liver_testis", 1, 0)
tmp$tissue_pair_bin <- as.factor(tmp$tissue_pair_bin)

# germline poised
x <- lm(logFC_recent_cons ~ mean_germlineRT + tissue_pair_bin + tissue_pair_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Poised",])
coef(summary(x))
#                                    Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)                       0.3524864 0.02392230  14.73464 8.205775e-16
# mean_germlineRT                  -0.7090734 0.03161827 -22.42607 3.874590e-21
# tissue_pair_bin1                 -0.3194734 0.02392230 -13.35463 1.244498e-14
# mean_germlineRT:tissue_pair_bin1  0.2140035 0.03161827   6.76835 1.197464e-07

#germline active
x <- lm(logFC_recent_cons ~ mean_germlineRT + tissue_pair_bin + tissue_pair_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Active",])
coef(summary(x))
#                                    Estimate Std. Error    t value     Pr(>|t|)
# (Intercept)                      -0.1259383 0.03031153  -4.154798 2.262693e-04
# mean_germlineRT                  -0.6708828 0.04006295 -16.745718 2.160681e-17
# tissue_pair_bin1                 -0.2818807 0.03031153  -9.299456 1.299605e-10
# mean_germlineRT:tissue_pair_bin1  0.1942978 0.04006295   4.849813 3.073348e-05

#somatic poised
x <- lm(logFC_recent_cons ~ mean_somaticRT + tissue_pair_bin + tissue_pair_bin:mean_somaticRT
        , data = tmp[tmp$enh_type=="Poised",])
coef(summary(x))
#                                   Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)                      0.2883561 0.05834178  4.942531 2.348268e-05
# mean_somaticRT                  -0.5513938 0.07259889 -7.595072 1.185335e-08
# tissue_pair_bin1                -0.2999345 0.05834178 -5.140990 1.318708e-05
# mean_somaticRT:tissue_pair_bin1  0.1637956 0.07259889  2.256172 3.102605e-02

#somatic active
x <- lm(logFC_recent_cons ~ mean_somaticRT + tissue_pair_bin + tissue_pair_bin:mean_somaticRT
        , data = tmp[tmp$enh_type=="Active",])
coef(summary(x))
#                                   Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)                     -0.1889033 0.06398119 -2.952482 5.862310e-03
# mean_somaticRT                  -0.4890916 0.07961641 -6.143101 7.177568e-07
# tissue_pair_bin1                -0.2640832 0.06398119 -4.127514 2.444460e-04
# mean_somaticRT:tissue_pair_bin1  0.1478902 0.07961641  1.857534 7.245859e-02
