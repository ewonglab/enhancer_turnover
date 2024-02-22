library(pheatmap)
library("mclust")
library(data.table)
library("GenomeInfoDb")
library("GenomicFeatures")
library("AnnotationFilter")
library(GenomicRanges)
library(RColorBrewer)

# read bins with somatic and germ line RT (k = 18 clusters defined based on germ line and somatic RT)
setwd("./replication_timing/mouse")
mouse_rt <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header = T
                       , stringsAsFactors = F, sep = '\t')
mouse_rt <- mouse_rt[complete.cases(mouse_rt),]#

# add bin coordinates
mouse_rt$chr <- sub(":.*", "", rownames(mouse_rt))
mouse_rt$start <- sub("_.*", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$end <- sub(".*_", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$start <- as.integer(mouse_rt$start)
mouse_rt$end <- as.integer(mouse_rt$end)
mouse_rt$peakid <- rownames(mouse_rt)

# read enhancers 
enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                  , header = F, stringsAsFactors = F, sep ='\t')
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

enh <- enh[enh$V6 == "Recent" | enh$V4 %in% cons_enh$V1, ] 

# add a column indicating tissue specificity
enh$is_specific <- ifelse(grepl(",", enh$V7), 0, 1)

enh_cons <- unique(enh[enh$V6 == "Conserved", c("V1", "V2", "V3", "is_specific")])
enh_recent <- unique(enh[enh$V6 == "Recent", c("V1", "V2", "V3", "is_specific")])

## make gr objects of RT bins and enhancers
# make gr objects
enh_cons$enh_id <- with(enh_cons, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))
enh_recent$enh_id <- with(enh_recent, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))

enh_cons.gr <- with(enh_cons, GRanges(V1 , IRanges( V2+1, V3 )))
enh_recent.gr <- with(enh_recent, GRanges(V1 , IRanges( V2+1, V3 )))

mouse_rt.gr <- with(mouse_rt, GRanges(chr , IRanges( start+1, end )))

# overlap RT bins and enhancers
x <- as.data.frame(findOverlaps(mouse_rt.gr, enh_cons.gr))
RT_cons <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop = F]
                 , enh_cons[x$subjectHits, c("enh_id", "is_specific")])

x <- as.data.frame(findOverlaps(mouse_rt.gr, enh_recent.gr))
RT_recent <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop = F]
                   , enh_recent[x$subjectHits, c("enh_id", "is_specific")])

# count unique enhancers per cluster
RT_cons <- unique(RT_cons) # 71614     3
RT_recent <- unique(RT_recent) # 63692     3

# frequency of tissue specific and non-tissue specific enhancers per RT bin
# add a column to specify tissue specificity

cons.freq <- as.data.frame(table(RT_cons[,c("kmeans.cluster", "is_specific")]))
recent.freq <- as.data.frame(table(RT_recent[,c("kmeans.cluster", "is_specific")]))

cons.freq$is_specific <- ifelse(cons.freq$is_specific == 1, "specific", "non_specific")
recent.freq$is_specific <- ifelse(recent.freq$is_specific == 1, "specific", "non_specific")


# library(tidyr)
cons.freq <- tidyr::spread(cons.freq, is_specific, Freq)
recent.freq <- tidyr::spread(recent.freq, is_specific, Freq)

colnames(cons.freq)[2:3] <- paste0(colnames(cons.freq)[2:3], "_cons")
colnames(recent.freq)[2:3] <- paste0(colnames(recent.freq)[2:3], "_recent")

all_freq <- merge(cons.freq, recent.freq, by = "kmeans.cluster", all=T)#


# calculate enhancer turnover rate for 1) tissue-specific and 2) NON-tissue-specific enhancers
# use direction recent / conserved
all_freq$enhrate_specific <- with(all_freq, log((specific_recent+1)/(specific_cons+1)))
all_freq$enhrate_nonspecific <- with(all_freq, log((non_specific_recent+1)/(non_specific_cons+1)))

## add mean replication time back - ALL SAMPLES
rt_columns <- c("X46C", "D3", "TT2", "iPSC", "iPSC_1D4", "iPSC_2D4", "EPL", "EBM3"
                , "EpiSC5", "EpiSC7", "EBM6", "X46CNPC", "TT2NPC", "EBM9", "Mesoderm"
                , "Endoderm", "piPSC_1A2", "piPSC_1B3", "piPSC_V3", "MEF_female"
                , "MEF_male", "Myoblast", "PGC.female.1", "PGC.male.1", "sperm.1", "sperm.2")#26

## ONLY GERM LINE RT
rt_columns_germ <- c("PGC.female.1", "PGC.male.1", "sperm.1", "sperm.2")# 4


# mean across the 26 samples
mouse_rt$mean_rt_all <- apply(mouse_rt[,rt_columns], 1, mean)
mouse_rt$mean_rt_germ <- apply(mouse_rt[,rt_columns_germ], 1, mean)

mean_rt <- aggregate(cbind(mean_rt_all, mean_rt_germ) ~ kmeans.cluster, mouse_rt, mean)

all_freq <- merge(all_freq, mean_rt, by = "kmeans.cluster", all.x = T)# 18  9

write.table(x = all_freq, file = "tissue_specificity_meanRT_enhRate_kmeans18"
            , quote = F, sep = '\t')



library(ggplot2)
library(reshape2)
setwd("./replication_timing/mouse")
all_freq <- read.table(file = "tissue_specificity_meanRT_enhRate_kmeans18", header = T
                       , stringsAsFactors = F, sep = '\t')
# use mean germ line RT
df <- melt(all_freq[,c("enhrate_specific", "enhrate_nonspecific", "mean_rt_germ", "kmeans.cluster")]
           , c("mean_rt_germ", "kmeans.cluster"))#  36  4


pdf("meanGermRT_vs_enhRate_by_tissueSpec_by_kmeansClust_k18_cons2sp.pdf")
ggplot(df, aes(x = mean_rt_germ, y = value, color = variable)) +
  geom_point()  + theme_classic() +  geom_smooth(method=lm)
dev.off()

# ancova test
library(rstatix)
library(reshape2)

all_freq <- read.table(file = "tissue_specificity_meanRT_enhRate_kmeans18", header = T
                       , stringsAsFactors = F, sep = '\t')
df <- melt(all_freq[,c("enhrate_specific", "enhrate_nonspecific", "mean_rt_germ", "kmeans.cluster")]
           , c("mean_rt_germ", "kmeans.cluster"))


# ancova test for the difference in slopes

anova_result <- df %>% anova_test(value ~ variable + mean_rt_germ + variable:mean_rt_germ)
anova_result

# ANOVA Table (type II tests)
#
# Effect DFn DFd       F        p p<.05   ges
# 1              variable   1  32 295.114 1.03e-17     * 0.902
# 2          mean_rt_germ   1  32 286.343 1.60e-17     * 0.899
# 3 variable:mean_rt_germ   1  32  26.561 1.27e-05     * 0.454

# correlations
rsq <- function (x, y) cor(x, y) ^ 2
rsq(all_freq$mean_rt_germ, all_freq$enhrate_specific)
# 0.9512257
rsq(all_freq$mean_rt_germ, all_freq$enhrate_nonspecific)
# 0.7802071

x <- cor.test(all_freq$mean_rt_germ, all_freq$enhrate_specific)
x$estimate
# cor
# -0.975308
x$estimate ** 2
# cor
# 0.9512257
x$p.value
# 6.430635e-12

x <- cor.test(all_freq$mean_rt_germ, all_freq$enhrate_nonspecific)
x$estimate
# cor
# -0.8832933
x$estimate ** 2
# cor
# 0.7802071
x$p.value
# 1.192687e-06
