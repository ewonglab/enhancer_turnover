library(GenomicRanges)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)

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

# mean RT
mouse_RT$mean_germlineRT <- apply(mouse_RT[,germline], 1, mean)
mouse_RT$mean_somaticRT <- apply(mouse_RT[,somatic], 1, mean)

mean_mouseRT <- aggregate(.~kmeans.cluster
                          , mouse_RT[,c("mean_germlineRT", "mean_somaticRT","kmeans.cluster")]
                          , mean)

# kmeans.cluster mean_germlineRT mean_somaticRT
# 1               1       0.2557718     -0.3377223
# 2               2      -1.5276388     -0.6201144
# 3               3       0.5344305      0.5811477
# 4               4       1.1588596      1.5392447
# 5               5      -0.2522994     -0.7056573
# 6               6       1.0143551      1.3136903
# 7               7      -1.3777489     -0.8931174
# 8               8      -0.3946769     -0.7273448
# 9               9       0.5564480      0.5738714
# 10             10       0.7491716      0.6339254
# 11             11       0.8332980      1.0137935
# 12             12      -0.6998528     -1.0728049
# 13             13       0.3041689      0.1175160
# 14             14       0.2708385     -0.1796346
# 15             15       0.1021640     -0.8683291
# 16             16       0.8287007      1.0084502
# 17             17       0.4266030      0.2032035
# 18             18      -0.1720332     -0.3165335

# ENHANCERS
mouse.enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                        , header = F, stringsAsFactors = F
                        , sep = "\t")
# keep only recent enhancers
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

# SPLIT COORDINATES OF CONSERVED ENHANCERS
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

### adding other tissues and types of enhancers
#poised  - all tissues
poised_recent.long <- recent_enh[recent_enh$V5 == "Poised",]
poised_cons <- cons_enh[cons_enh$type == "poised",]

poised_recent.long$enh_type <- "recent"
poised_cons$enh_type <- "conserved"
poised_recent.long <- unique(poised_recent.long)
poised_cons <- unique(poised_cons)
tmp <- data.frame(chr = c(poised_recent.long$V1, poised_cons$chr)
                  , start = c(poised_recent.long$V2, poised_cons$start)
                  , end = c(poised_recent.long$V3, poised_cons$end)
                  #ignore tissue
                  # , tissue = c(tolower(poised_recent.long$tissue), poised_cons$type)
                  , enh_type = c(poised_recent.long$enh_type, poised_cons$enh_type)
                  , stringsAsFactors = F)
tmp <- unique(tmp)
#remove tissue and keep unique

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp))
mean_mouseRT$all.conserved.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$all.recent.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

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
                  #ignore tissue
                  # , tissue = c(tolower(active_recent.long$tissue), active_cons$type)
                  , enh_type = c(active_recent.long$enh_type, active_cons$enh_type)
                  , stringsAsFactors = F)#
tmp <- unique(tmp)

#remove tissue and keep unique
cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp))
mean_mouseRT$all.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$all.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

### SCATTERPLOTS
# tmp <- mean_mouseRT[,c(1:3, 8:9)]
tmp <- melt(mean_mouseRT, id = c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT"))
# unique(tmp$variable)
# [1] all.conserved.poised all.recent.poised    all.conserved.active
# [4] all.recent.active

tmp$enh_type <- ifelse(grepl("active", tmp$variable), "Active", "Poised")
tmp$enh_type <- factor(tmp$enh_type, levels = c("Active", "Poised"))
tmp$age <- ifelse(grepl("conserved", tmp$variable), "conserved", "recent")
tmp$age <- factor(tmp$age, levels = c("conserved", "recent"))


pdf("mouse_rt_recent_cons_nEnh_allTissues_by_type_mm10_germ_v_somatic_cons2sp.pdf")
ggplot(tmp, aes(x=mean_germlineRT, y=log10(value), group=age, color=age)) +
  geom_point() + theme_classic() +
  geom_smooth(method="lm", se=T, fullrange=TRUE
              , aes(fill=age, color=age))+
  facet_wrap(~enh_type, nrow = 1, ncol = 2)
dev.off()

#testing the difference in slope
tmp$age_bin <- ifelse(tmp$age == "recent", 1, 0)
                                                
tmp$age_bin <- as.factor(tmp$age_bin)
tmp$log10_nEnh <- log10(tmp$value)
                                                
#all
x <- lm(log10_nEnh ~ mean_germlineRT + age_bin + age_bin:mean_germlineRT, data = tmp)
coef(summary(x))
#                            Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)               3.0605183 0.03710187 82.489587 6.753350e-70
# mean_germlineRT           0.5259409 0.04903780 10.725215 2.889851e-16
# age_bin1                  0.0643074 0.05246997  1.225604 2.245757e-01
# mean_germlineRT:age_bin1 -0.2967984 0.06934992 -4.279723 6.001833e-05

#poised
x <- lm(log10_nEnh ~ mean_germlineRT + age_bin + age_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Poised",])
coef(summary(x))
#                            Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)               3.1465802 0.03752144 83.860863 4.404080e-39
# mean_germlineRT           0.4672868 0.04959234  9.422561 9.517578e-11
# age_bin1                  0.1786614 0.05306333  3.366947 1.991259e-03
# mean_germlineRT:age_bin1 -0.3117735 0.07013415 -4.445388 9.876408e-05

#active
x <- lm(log10_nEnh ~ mean_germlineRT + age_bin + age_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Active",])
coef(summary(x))
#                             Estimate Std. Error    t value     Pr(>|t|)
# (Intercept)               2.97445636 0.03740745 79.5150827 2.398499e-38
# mean_germlineRT           0.58459500 0.04944168 11.8239319 3.233420e-13
# age_bin1                 -0.05004661 0.05290212 -0.9460228 3.512262e-01
# mean_germlineRT:age_bin1 -0.28182329 0.06992109 -4.0305907 3.213621e-04

