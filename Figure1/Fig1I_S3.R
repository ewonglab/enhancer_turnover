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
                        , sep = "\t")
recent_enh <- mouse.enh[mouse.enh$V6 == "Recent",]
active_atLeast4sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast4sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
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

# SPLIT COORDINATES OF CONSERVED ENHANCERS
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
                  , stringsAsFactors = F)# 94810     5
tmp$chr <- sub("chrG", "G", tmp$chr)
unique(tmp$chr)
# [1] "chr1"          "chr10"         "chr11"         "chr12"        
# [5] "chr13"         "chr14"         "chr15"         "chr16"        
# [9] "chr17"         "chr18"         "chr19"         "chr2"        
# [13] "chr3"          "chr4"          "chr5"          "chr6"        
# [17] "chr7"          "chr8"          "chr9"          "GL456216.1"  
# [21] "GL456233.1"    "GL456378.1"    "chrJH584304.1" "chrX"        
# [25] "chrY"

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="liver",]))
mean_mouseRT$liver.conserved.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$liver.recent.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="brain",]))
mean_mouseRT$brain.conserved.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$brain.recent.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="muscle",]))
mean_mouseRT$muscle.conserved.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$muscle.recent.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="testis",]))
mean_mouseRT$testis.conserved.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$testis.recent.poised <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

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
                  , stringsAsFactors = F)#67772     5

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="liver",]))
mean_mouseRT$liver.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$liver.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="brain",]))
mean_mouseRT$brain.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$brain.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="muscle",]))
mean_mouseRT$muscle.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$muscle.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

cons_rec <- lapply(1:18, function(x) cons.recent.ov.cluster(all, x, tmp[tmp$tissue=="testis",]))
mean_mouseRT$testis.conserved.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['conserved']))
mean_mouseRT$testis.recent.active <- unlist(lapply(1:18, function(x) cons_rec[[x]]['recent']))

### RECENT/CONSERVED  log FC TISSUES
mean_mouseRT$liver.poised_rec.cons_logFC <- with(mean_mouseRT
                                                 , log(liver.recent.poised/liver.conserved.poised))
mean_mouseRT$liver.active_rec.cons_logFC <- with(mean_mouseRT
                                                 , log(liver.recent.active/liver.conserved.active))

mean_mouseRT$brain.poised_rec.cons_logFC <- with(mean_mouseRT
                                                 , log(brain.recent.poised/brain.conserved.poised))
mean_mouseRT$brain.active_rec.cons_logFC <- with(mean_mouseRT
                                                 , log(brain.recent.active/brain.conserved.active))

mean_mouseRT$muscle.poised_rec.cons_logFC <- with(mean_mouseRT
                                                  , log(muscle.recent.poised/muscle.conserved.poised))
mean_mouseRT$muscle.active_rec.cons_logFC <- with(mean_mouseRT
                                                  , log(muscle.recent.active/muscle.conserved.active))

mean_mouseRT$testis.poised_rec.cons_logFC <- with(mean_mouseRT
                                                  , log(testis.recent.poised/testis.conserved.poised))
mean_mouseRT$testis.active_rec.cons_logFC <- with(mean_mouseRT
                                                  , log(testis.recent.active/testis.conserved.active))

### R2 values vs germline RT

rsq <- function (x, y) cor(x, y) ^ 2

## GERMLINE RT
round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$liver.poised_rec.cons_logFC), 2)
# [1] 0.93
round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$liver.active_rec.cons_logFC), 2)
# [1] 0.91

round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$brain.poised_rec.cons_logFC), 2)
# [1] 0.92
round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$brain.active_rec.cons_logFC), 2)
# [1] 0.86

round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$muscle.poised_rec.cons_logFC), 2)
# [1] 0.91
round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$muscle.active_rec.cons_logFC), 2)
# [1] 0.76
 
round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$testis.poised_rec.cons_logFC), 2)
# [1] 0.95
round(rsq(mean_mouseRT$mean_germlineRT, mean_mouseRT$testis.active_rec.cons_logFC), 2)
# [1] 0.88

# SOMATIC RT
round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$liver.poised_rec.cons_logFC), 2)
# [1] 0.62
round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$liver.active_rec.cons_logFC), 2)
# [1] 0.57

round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$brain.poised_rec.cons_logFC), 2)
# [1] 0.57
round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$brain.active_rec.cons_logFC), 2)
# [1] 0.44

round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$muscle.poised_rec.cons_logFC), 2)
# [1] 0.77
round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$muscle.active_rec.cons_logFC), 2)
# [1] 0.57
 
round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$testis.poised_rec.cons_logFC), 2)
# [1] 0.71
round(rsq(mean_mouseRT$mean_somaticRT, mean_mouseRT$testis.active_rec.cons_logFC), 2)
# [1] 0.58

write.table(x = mean_mouseRT, file = "mouse_mean_rt_logFC_by_tissue_mm10_consAtLeast2sp"
            , quote = F, sep = "\t")

### SCATTERPLOTS

library(reshape2)
tmp <- mean_mouseRT[,c(1:3, 20:27)]

tmp <- melt(tmp, id = c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT"))
tmp$variable <- factor(tmp$variable, levels = c("liver.poised_rec.cons_logFC", "liver.active_rec.cons_logFC"
                                                , "brain.poised_rec.cons_logFC", "brain.active_rec.cons_logFC"
                                                , "muscle.poised_rec.cons_logFC", "muscle.active_rec.cons_logFC"
                                                , "testis.poised_rec.cons_logFC", "testis.active_rec.cons_logFC"))

tmp$enh_type <- ifelse(grepl("active", tmp$variable), "Active", "Poised")
tmp$tissue <- ""
tmp[grep("liver", tmp$variable), "tissue"] <- "Liver"
tmp[grep("brain", tmp$variable), "tissue"] <- "Brain"
tmp[grep("muscle", tmp$variable), "tissue"] <- "Muscle"
tmp[grep("testis", tmp$variable), "tissue"] <- "Testis"
tmp$enh_type <- factor(tmp$enh_type, levels = c("Active", "Poised"))

colnames(tmp) <- c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT", "sample"
                   , "logFC_recent_cons", "enh_type", "tissue")

tmp2 <- melt(tmp, id=c("kmeans.cluster", "sample", "logFC_recent_cons", "enh_type", "tissue"))

#variable = somatic or germline
library(ggplot2)
pdf("mouse_rt_recent_cons_logFC_by_type_mm10_germ_v_somatic_cons2sp.pdf")
ggplot(tmp2, aes(x=value, y=logFC_recent_cons, group=1, color=tissue)) +
  geom_point() + theme_classic() +
  geom_smooth(method="lm", se=T, fullrange=TRUE
              , aes(fill=tissue, color=tissue, group=tissue))+
  facet_wrap(~enh_type + variable, nrow = 2, ncol = 2)
dev.off()

