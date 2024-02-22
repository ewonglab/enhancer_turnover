library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

setwd("./replication_timing/mouse")

# read enhancers (roller)
enh <- read.table(file = "all_enh_roller_overlapping_TE", header = T
                  , stringsAsFactors = F, sep ='\t')
# keep recent enhancers and correct conserved enhancers
active_atLeast2sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast2sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
cons_enh <- rbind(active_atLeast2sp, poised_atLeast2sp)
cons_enh$V1 <- paste0("chr", cons_enh$V1)

enh <- enh[enh$V6 == "Recent" | enh$V4 %in% cons_enh$V1, ] 

# percentage of enhancers that overlap TE

table(enh$overlap_anyTE)
# 0     1
# 80385 88256
round((88256/(88256 + 80385)) * 100) # 52

# read replication time data (including clustering membership for 18 kmeans clusters)

mouse_rt <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header = T
                       , stringsAsFactors = F, sep = '\t')
mouse_rt <- mouse_rt[complete.cases(mouse_rt),]

# add bin coordinates
mouse_rt$chr <- sub(":.*", "", rownames(mouse_rt))
mouse_rt$start <- sub("_.*", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$end <- sub(".*_", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$start <- as.integer(mouse_rt$start)
mouse_rt$end <- as.integer(mouse_rt$end)
mouse_rt$peakid <- rownames(mouse_rt)

# separate enhancers into:
# recent that overlap TE
# recent that do not overlap TE
# conserved that overlap TE
# conserved that dont overlap TE

recent_TE <- enh[enh$V6 == "Recent" & enh$overlap_anyTE == 1, ] 
recent_notTE <- enh[enh$V6 == "Recent" & enh$overlap_anyTE == 0, ] 
cons_TE <- enh[enh$V6 == "Conserved" & enh$overlap_anyTE == 1, ] 
cons_notTE <- enh[enh$V6 == "Conserved" & enh$overlap_anyTE == 0, ] 

# make gr objects
recent_TE.gr <- with(recent_TE, GRanges(V1 , IRanges( V2+1, V3 )))
recent_notTE.gr <- with(recent_notTE, GRanges(V1 , IRanges( V2+1, V3 )))
cons_TE.gr <- with(cons_TE, GRanges(V1 , IRanges( V2+1, V3 )))
cons_notTE.gr <- with(cons_notTE, GRanges(V1 , IRanges( V2+1, V3 )))

mouse_rt.gr <- with(mouse_rt, GRanges(chr , IRanges( start+1, end )))

# overlap RT bins and enhancers
x <- as.data.frame(findOverlaps(mouse_rt.gr, recent_TE.gr))
RT_recent_TE <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop =F], recent_TE[x$subjectHits, "V4", drop=F])

x <- as.data.frame(findOverlaps(mouse_rt.gr, recent_notTE.gr))
RT_recent_notTE <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop =F], recent_notTE[x$subjectHits,"V4", drop=F])

x <- as.data.frame(findOverlaps(mouse_rt.gr, cons_TE.gr))
RT_cons_TE <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop =F], cons_TE[x$subjectHits,"V4", drop=F])

x <- as.data.frame(findOverlaps(mouse_rt.gr, cons_notTE.gr))
RT_cons_notTE <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop =F], cons_notTE[x$subjectHits,"V4", drop=F])

# frequency of each enhancer type per kmeans cluster

RT_recent_TE <- unique(RT_recent_TE)
RT_recent_notTE <- unique(RT_recent_notTE)
RT_cons_TE <- unique(RT_cons_TE)
RT_cons_notTE <- unique(RT_cons_notTE)

RT_recent_TE.freq <- as.data.frame(table(RT_recent_TE$kmeans.cluster))
RT_recent_notTE.freq <- as.data.frame(table(RT_recent_notTE$kmeans.cluster))
RT_cons_TE.freq <- as.data.frame(table(RT_cons_TE$kmeans.cluster))
RT_cons_notTE.freq <- as.data.frame(table(RT_cons_notTE$kmeans.cluster))

# add class
RT_recent_TE.freq$class <- "recent_TE"
RT_recent_notTE.freq$class <- "recent_notTE"
RT_cons_TE.freq$class <- "cons_TE"
RT_cons_notTE.freq$class <- "cons_notTE"

all_enh_freq <- rbind(RT_recent_TE.freq, RT_recent_notTE.freq, RT_cons_TE.freq, RT_cons_notTE.freq)
all_enh_freq <- tidyr::spread(all_enh_freq, class, Freq)

# add mean RT
rt_columns <- c("PGC.female.1", "PGC.male.1", "sperm.1", "sperm.2")

mouse_rt$mean_rt <- apply(mouse_rt[,rt_columns], 1, mean)
mean_rt <- aggregate(mean_rt ~ kmeans.cluster, mouse_rt, mean)

all_enh_freq <- merge(all_enh_freq, mean_rt,by.x="Var1", by.y="kmeans.cluster", all.x = T)

# calculate enhancer turnover rates for 1) enhancers that overlap TE and 2) enhancers that do not overlap TE
all_enh_freq$enhrate_TE <- with(all_enh_freq, log((recent_TE+1)/(cons_TE+1)))
all_enh_freq$enhrate_notTE <- with(all_enh_freq, log((recent_notTE+1)/(cons_notTE+1)))

write.table(x = all_enh_freq, file = "allEnh_by_overlapwithTE_meanGermRT_enhRate_kmeans18_cons2sp"
            , quote = F, sep = '\t')


### plot enhancer rate versus GERMLINE mean DNA rep time

library(ggplot2)
library(reshape2)

all_enh_freq <- read.table(file = "allEnh_by_overlapwithTE_meanGermRT_enhRate_kmeans18_cons2sp"
                           , header = T, sep = '\t', stringsAsFactors = F)
rsq <- function (x, y) cor(x, y) ^ 2

rsq(all_enh_freq$mean_rt, all_enh_freq$enhrate_TE)
# [1] 0.9385226
rsq(all_enh_freq$mean_rt, all_enh_freq$enhrate_notTE)
# [1] 0.9464222

x <- cor.test(all_enh_freq$mean_rt, all_enh_freq$enhrate_TE)
x$p.value
# [1] 4.121298e-11
x$estimate
# cor
# -0.9687737
x$estimate ^ 2
# cor
# 0.9385226

x <- cor.test(all_enh_freq$mean_rt, all_enh_freq$enhrate_notTE)
x$p.value
# [1] 1.366399e-11
x$estimate
# cor
# -0.9728423
x$estimate ^ 2
# cor
# 0.9464222


all_enh_freq_melt <- melt(all_enh_freq[,c("mean_rt", "enhrate_TE", "enhrate_notTE", "Var1")], c("mean_rt", "Var1"))

p <- ggplot(all_enh_freq_melt, aes(x = mean_rt, y = value, color = variable, group = variable)) +
  geom_point()  + theme_classic() +  geom_smooth(method=lm, aes(fill = variable))


pdf("meanGermRT_vs_enhRate_by_overlapTE_by_kmeansClust_k18_cons2sp.pdf")
p
dev.off()

library(rstatix)
library(reshape2)

all_enh_freq <- read.table(file = "allEnh_by_overlapwithTE_meanGermRT_enhRate_kmeans18_cons2sp"
                           , header = T, stringsAsFactors = F, sep = '\t')
# Var1 = kmeans cluster
df <- melt(all_enh_freq[,c("Var1", "enhrate_TE", "enhrate_notTE", "mean_rt")]
           , c("Var1", "mean_rt"))# 36  4


# ancova test for the difference in slopes
# I am testing whether the mean RT has different effects on the enhancer rate depending on whether they are tissue
# specific or non-tissue-specific enhancers

anova_result <- df %>% anova_test(value ~ variable + mean_rt + variable:mean_rt)
anova_result

# ANOVA Table (type II tests)
#
#             Effect DFn DFd       F        p p<.05   ges
# 1         variable   1  32   1.482 2.32e-01       0.044
# 2          mean_rt   1  32 504.452 3.70e-21     * 0.940
# 3 variable:mean_rt   1  32   9.538 4.00e-03     * 0.230
