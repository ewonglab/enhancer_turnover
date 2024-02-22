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

#mean RT
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
active_atLeast4sp <- read.table(file = "./roller/active_conserved_atLeast4sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")#
poised_atLeast4sp <- read.table(file = "./roller/poised_conserved_atLeast4sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t") #
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

# SPLIT COORDINATES OF CONSERVED ENHANCERS
cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end)

#bins and cluster membership
all <- mouse_RT[,"kmeans.cluster",drop =F]
all$peakid <- rownames(all)
all$chr <- gsub(":.*", "", all$peakid)
all$start <- gsub("_.*", "", gsub(".*:", "", all$peakid))
all$end <- gsub(".*_", "", gsub(".*:", "", all$peakid))
all$start <- as.integer(all$start)
all$end <- as.integer(all$end)
all$width <- with(all, end-start)

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
                  , enh_type = c(poised_recent.long$enh_type, poised_cons$enh_type)
                  , stringsAsFactors = F)
tmp <- unique(tmp)

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
                  , enh_type = c(active_recent.long$enh_type, active_cons$enh_type)
                  , stringsAsFactors = F)
tmp <- unique(tmp)

#remove tissue and keep unique
cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp))
mean_mouseRT$all.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$all.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

##plot number of recent and conserved enhancers -poised vs active


### SCATTERPLOTS
tmp <- melt(mean_mouseRT, id = c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT"))

tmp$enh_type <- ifelse(grepl("active", tmp$variable), "Active", "Poised")
tmp$enh_type <- factor(tmp$enh_type, levels = c("Active", "Poised"))
tmp$age <- ifelse(grepl("conserved", tmp$variable), "conserved", "recent")
tmp$age <- factor(tmp$age, levels = c("conserved", "recent"))

pdf("mouse_rt_recent_cons_nEnh_allTissues_by_type_mm10_germ_v_somatic_cons4sp.pdf")
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
#                           Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)               2.7710273 0.03586053 77.272339 5.481941e-68
# mean_germlineRT           0.5848212 0.04739711 12.338754 5.105620e-19
# age_bin1                  0.3537984 0.05071445  6.976284 1.565065e-09
# mean_germlineRT:age_bin1 -0.3556787 0.06702963 -5.306291 1.320941e-06

#poised
x <- lm(log10_nEnh ~ mean_germlineRT + age_bin + age_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Poised",])
coef(summary(x))
#                           Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)               2.8140369 0.03736809 75.305887 1.354867e-37
# mean_germlineRT           0.5283690 0.04938965 10.697971 4.241675e-12
# age_bin1                  0.5112047 0.05284645  9.673396 5.076284e-11
# mean_germlineRT:age_bin1 -0.3728558 0.06984751 -5.338140 7.428075e-06

#active
x <- lm(log10_nEnh ~ mean_germlineRT + age_bin + age_bin:mean_germlineRT
        , data = tmp[tmp$enh_type=="Active",])
coef(summary(x))
#                           Estimate Std. Error   t value     Pr(>|t|)
# (Intercept)               2.7280176 0.03669509 74.342846 2.040713e-37
# mean_germlineRT           0.6412734 0.04850015 13.222092 1.632894e-14
# age_bin1                  0.1963921 0.05189470  3.784436 6.391404e-04
# mean_germlineRT:age_bin1 -0.3385017 0.06858957 -4.935178 2.398955e-05
