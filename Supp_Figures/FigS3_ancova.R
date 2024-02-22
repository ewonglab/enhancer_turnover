library(GenomicRanges)

cons.recent.ov.cluster <- function(clust_bed, query_clust, enh_tab){
  clust <- subset(clust_bed, kmeans.cluster==query_clust)
  clust_gr <- with(clust, GRanges( chr , IRanges( start+1, end )))
  enh_gr <- with(enh_tab, GRanges( chr , IRanges( start+1, end )))
  x <-  as.data.frame(findOverlaps(enh_gr, clust_gr))
  result <- enh_tab[unique(x$queryHits),]
  return(table(result$enh_type))}

setwd("/g/data/zk16/cc3704/replication_timing/mouse")
mouse_RT <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header =T
                       , stringsAsFactors = F, sep = '\t')

somatic <- c("X46C", "D3", "TT2", "iPSC", "iPSC_1D4"
             , "iPSC_2D4", "EPL", "EBM3", "EpiSC5", "EpiSC7"
             , "EBM6", "X46CNPC", "TT2NPC", "EBM9", "Mesoderm"
             , "Endoderm", "piPSC_1A2", "piPSC_1B3", "piPSC_V3"
             , "MEF_female", "MEF_male", "Myoblast")
germline <- c("sperm.1", "sperm.2", "PGC.male.1", "PGC.female.1")

#mean RT
mouse_RT$mean_germlineRT <- apply(mouse_RT[,germline], 1, mean)
mouse_RT$mean_somaticRT <- apply(mouse_RT[,somatic], 1, mean)

mean_mouseRT <- aggregate(.~kmeans.cluster
                          , mouse_RT[,c("mean_germlineRT", "mean_somaticRT","kmeans.cluster")]
                          , mean)


# ENHANCERS
mouse.enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                        , header = F, stringsAsFactors = F
                        , sep = "\t")#
recent_enh <- mouse.enh[mouse.enh$V6 == "Recent",]
# conserved in min 4 other species
active_atLeast4sp <- read.table(file = "./roller/active_conserved_atLeast4sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast4sp <- read.table(file = "./roller/poised_conserved_atLeast4sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
cons_enh <- rbind(active_atLeast4sp, poised_atLeast4sp)

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

# coordinated of conserved enh
cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end)

#200kb bins and cluster membership
all <- mouse_RT[,"kmeans.cluster",drop =F]
all$peakid <- rownames(all)
all$chr <- gsub(":.*", "", all$peakid)
all$start <- gsub("_.*", "", gsub(".*:", "", all$peakid))
all$end <- gsub(".*_", "", gsub(".*:", "", all$peakid))
all$start <- as.integer(all$start)
all$end <- as.integer(all$end)
all$width <- with(all, end-start)

#poised  - all tissues
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
unique(tmp$chr)
# [1] "chr1"          "chr10"         "chr11"         "chr12"        
# [5] "chr13"         "chr14"         "chr15"         "chr16"        
# [9] "chr17"         "chr18"         "chr19"         "chr2"        
# [13] "chr3"          "chr4"          "chr5"          "chr6"        
# [17] "chr7"          "chr8"          "chr9"          "GL456216.1"  
# [21] "GL456233.1"    "GL456378.1"    "chrJH584304.1" "chrX"        
# [25] "chrY"        

#separate into tissue pairs
tmp$tissue_pair <- ""
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("liver", "testis"), "liver_testis", tmp$tissue_pair)
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("brain", "muscle"), "brain_muscle", tmp$tissue_pair)

tmp$tissue <- NULL
tmp <- unique(tmp)
table(tmp$tissue_pair)
# brain_muscle liver_testis
#       40427        47644

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
                  , stringsAsFactors = F)

#separate into tissue pairs
tmp$tissue_pair <- ""
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("liver", "testis"), "liver_testis", tmp$tissue_pair)
tmp$tissue_pair <- ifelse(tmp$tissue %in% c("brain", "muscle"), "brain_muscle", tmp$tissue_pair)

tmp$tissue <- NULL
tmp <- unique(tmp)
table(tmp$tissue_pair)
# brain_muscle liver_testis
# 34148        24328

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

write.table(x = mean_mouseRT, file = "mouse_mean_rt_logFC_by_tissuePairs_mm10_consAlLeast4sp"
            , quote = F, sep = "\t")

library(rstatix)
library(reshape2)
mean.RT.2.clust <- read.table(file = "mouse_mean_rt_logFC_by_tissuePairs_mm10_consAlLeast4sp"
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

tmp2 <- melt(tmp, id=c("kmeans.cluster", "group", "logFC_recent_cons", "enh_type"))
tmp2$tissue_pair <- ifelse(grepl("liver_testis", tmp2$group), "liver_testis", "brain_muscle")

# germline poised
tmp2[tmp2$variable=="mean_germlineRT" & tmp2$enh_type=="Poised",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd       F        p p<.05   ges
# 1       tissue_pair   1  32 109.940 7.02e-12     * 0.775
# 2             value   1  32 505.417 3.60e-21     * 0.940
# 3 tissue_pair:value   1  32  27.949 8.63e-06     * 0.466

#germline active
tmp2[tmp2$variable=="mean_germlineRT" & tmp2$enh_type=="Active",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd       F        p p<.05   ges
# 1       tissue_pair   1  32  60.211 7.55e-09     * 0.653
# 2             value   1  32 295.258 1.03e-17     * 0.902
# 3 tissue_pair:value   1  32  22.667 3.98e-05     * 0.415

#germline all
tmp2[tmp2$variable=="mean_germlineRT",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd       F        p p<.05   ges
# 1       tissue_pair   1  68  38.900 3.26e-08     * 0.364
# 2             value   1  68 184.570 4.72e-21     * 0.731
# 3 tissue_pair:value   1  68  12.082 8.92e-04     * 0.151

#somatic poised
tmp2[tmp2$variable=="mean_somaticRT" & tmp2$enh_type=="Poised",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd      F        p p<.05   ges
# 1       tissue_pair   1  32 19.273 1.16e-04     * 0.376
# 2             value   1  32 63.421 4.33e-09     * 0.665
# 3 tissue_pair:value   1  32  3.692 6.40e-02       0.103

#somatic active
tmp2[tmp2$variable=="mean_somaticRT" & tmp2$enh_type=="Active",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd      F        p p<.05   ges
# 1       tissue_pair   1  32 12.428 1.00e-03     * 0.280
# 2             value   1  32 37.428 7.72e-07     * 0.539
# 3 tissue_pair:value   1  32  2.797 1.04e-01       0.080

#somatic all
tmp2[tmp2$variable=="mean_somaticRT",] %>%
  anova_test(logFC_recent_cons ~ tissue_pair*value)

# ANOVA Table (type II tests)
#
# Effect DFn DFd      F        p p<.05   ges
# 1       tissue_pair   1  68 19.754 3.33e-05     * 0.225
# 2             value   1  68 62.296 3.43e-11     * 0.478
# 3 tissue_pair:value   1  68  4.099 4.70e-02     * 0.057

# Using linear regression

mean.RT.2.clust <- read.table(file = "mouse_mean_rt_logFC_by_tissuePairs_mm10_consAlLeast4sp"
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

#germline poised
x <- lm(logFC_recent_cons ~ mean_germlineRT + tissue_pair_bin + tissue_pair_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Poised",])
coef(summary(x))
#                                   Estimate Std. Error    t value
# (Intercept)                       0.7676904 0.04037444  19.014267
# mean_germlineRT                  -0.6488208 0.05336317 -12.158587
# tissue_pair_bin1                  0.6454475 0.05709808  11.304189
# mean_germlineRT:tissue_pair_bin1 -0.3989670 0.07546692  -5.286647
# Pr(>|t|)
# (Intercept)                      5.318647e-19
# mean_germlineRT                  1.550081e-13
# tissue_pair_bin1                 1.040380e-12
# mean_germlineRT:tissue_pair_bin1 8.629700e-06

#germline active
x <- lm(logFC_recent_cons ~ mean_germlineRT + tissue_pair_bin + tissue_pair_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Active",])
coef(summary(x))
# Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)                       0.1156185 0.05149347  2.245304 3.178745e-02
# mean_germlineRT                  -0.5978149 0.06805926 -8.783741 4.896856e-10
# tissue_pair_bin1                  0.6210566 0.07282276  8.528330 9.567350e-10
# mean_germlineRT:tissue_pair_bin1 -0.4582479 0.09625033 -4.761001 3.975355e-05

#somatic poised
x <- lm(logFC_recent_cons ~ mean_somaticRT + tissue_pair_bin + tissue_pair_bin:mean_somaticRT
        , data = tmp[tmp$enh_type=="Poised",])
coef(summary(x))
# Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)                      0.7090481 0.09500348  7.463391 1.704923e-08
# mean_somaticRT                  -0.5050900 0.11821968 -4.272469 1.619508e-04
# tissue_pair_bin1                 0.6101370 0.13435520  4.541223 7.498460e-05
# mean_somaticRT:tissue_pair_bin1 -0.3212587 0.16718788 -1.921543 6.360782e-02

#somatic active
x <- lm(logFC_recent_cons ~ mean_somaticRT + tissue_pair_bin + tissue_pair_bin:mean_somaticRT
        , data = tmp[tmp$enh_type=="Active",])
coef(summary(x))
#                                   Estimate Std. Error    t value    Pr(>|t|)
# (Intercept)                      0.05957954  0.1116682  0.5335409 0.597346004
# mean_somaticRT                  -0.43679727  0.1389568 -3.1434040 0.003589524
# tissue_pair_bin1                 0.57766769  0.1579226  3.6579154 0.000905840
# mean_somaticRT:tissue_pair_bin1 -0.32865537  0.1965145 -1.6724226 0.104192174
