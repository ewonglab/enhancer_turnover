library(pheatmap)
library("mclust")
library(data.table)
library("GenomeInfoDb")
library("GenomicFeatures")
library("AnnotationFilter")
library(GenomicRanges)
library(RColorBrewer)

setwd("./replication_timing/mouse")
mouse_rt <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header = T
                       , stringsAsFactors = F, sep = '\t')
mouse_rt$mean.rt <- NULL
mouse_rt <- mouse_rt[complete.cases(mouse_rt),]

# bin coordinates
mouse_rt$chr <- sub(":.*", "", rownames(mouse_rt))
mouse_rt$start <- sub("_.*", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$end <- sub(".*_", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$start <- as.integer(mouse_rt$start)
mouse_rt$end <- as.integer(mouse_rt$end)
mouse_rt$peakid <- rownames(mouse_rt)

# read conserved enhancers
active_atLeast2sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast2sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
cons_enh <- rbind(active_atLeast2sp, poised_atLeast2sp)
cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end)

# RECENT ENHANCERS
enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                  , header = F, stringsAsFactors = F, sep ='\t')#
enh_recent <- enh[enh$V6 == "Recent", ]#

# read enhancers mapped to hg19
enh_hg19 <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue_hg19.bed"
                       , header = F, stringsAsFactors = F, sep ='\t')

# read conserved enhancers
# for mouse-specific enhancers separate between those that align to hg19 and those that don't align
# to hg19
enh_recent_align <- enh_recent[enh_recent$V4 %in% enh_hg19$V4, ]
enh_recent_dontalign <- enh_recent[!enh_recent$V4 %in% enh_hg19$V4, ]

## make gr objects of RT bins and enhancers
# make gr objects
cons_enh.gr <- with(cons_enh, GRanges(chr , IRanges( start+1, end )))
enh_recent_align.gr <- with(enh_recent_align, GRanges(V1 , IRanges( V2+1, V3 )))
enh_recent_dontalign.gr <- with(enh_recent_dontalign, GRanges(V1 , IRanges( V2+1, V3 )))
mouse_rt.gr <- with(mouse_rt, GRanges(chr , IRanges( start+1, end )))

# overlap RT bins and enhancers
# count unique enhancers overlapping each cluster
x <- as.data.frame(findOverlaps(mouse_rt.gr, cons_enh.gr))
RT_cons <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop = F]
                 , cons_enh[x$subjectHits, "V1", drop = F])# V1 = enhancer id

x <- as.data.frame(findOverlaps(mouse_rt.gr, enh_recent_align.gr))
RT_recent_align <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop = F]
                         , enh_recent_align[x$subjectHits, "V4", drop = F])

x <- as.data.frame(findOverlaps(mouse_rt.gr, enh_recent_dontalign.gr))
RT_recent_dontalign <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop = F]
                             , enh_recent_dontalign[x$subjectHits, "V4", drop=F]) 

# unique enhancers per cluster
RT_cons <- unique(RT_cons)=
RT_recent_align <- unique(RT_recent_align) 
RT_recent_dontalign <- unique(RT_recent_dontalign)

# frequency of enhancers per RT bin
cons.freq <- as.data.frame(table(RT_cons$kmeans.cluster))
recent_align.freq <- as.data.frame(table(RT_recent_align$kmeans.cluster))
recent_dontalign.freq <- as.data.frame(table(RT_recent_dontalign$kmeans.cluster))

# change column names
colnames(cons.freq) <- c("kmeans.cluster", "N_cons")
colnames(recent_align.freq) <- c("kmeans.cluster", "N_recent_align")
colnames(recent_dontalign.freq) <- c("kmeans.cluster", "N_recent_dontalign")

all_freq <- merge(cons.freq, recent_align.freq, by = "kmeans.cluster")
all_freq <- merge(all_freq, recent_dontalign.freq, by = "kmeans.cluster")
all_freq <- all_freq[with(all_freq, order(kmeans.cluster)), ]

#     kmeans.cluster N_cons N_recent_align N_recent_dontalign
# 1               1   3286           2272               1459
# 11              2    195            302                418
# 12              3   3542           1583               1767
# 13              4  12517           2927               4598
# 14              5   2174           1689                993
# 15              6  12670           3100               4229
# 16              7    709           1287               1659
# 17              8   1755           1459               1425
# 18              9   3078           1197               1222
# 2              10   4137           1513               1292
# 3              11   5045           1812               2337
# 4              12   1877           2361               1748
# 5              13   2986           1460               1675
# 6              14   3206           1525               1112
# 7              15   2927           1923               1122
# 8              16   6432           2407               2578
# 9              17   3622           1682               1271
# 10             18   1456           1028               1260

# calculate enhancer turnover rate for 1) enhancers that align and 2) enhancers that don't align
all_freq$enhrate_recent_align <- with(all_freq, log((N_recent_align+1)/(N_cons+1)))
all_freq$enhrate_recent_dontalign <- with(all_freq, log((N_recent_dontalign+1)/(N_cons+1)))


## add mean replication time back
mouse_rt$mean_rt <- apply(mouse_rt[,c("PGC.female.1", "PGC.male.1", "sperm.1", "sperm.2")], 1, mean)
mean_rt <- aggregate(mean_rt~kmeans.cluster, mouse_rt, mean)

all_freq <- merge(all_freq, mean_rt, by="kmeans.cluster")
write.table(x = all_freq, file = "mean_germRT_recentEnh_align_hg19_mimMatch0.6_kmeans18_cons2sp"
            , quote = F, sep = '\t')


library(ggplot2)
library(reshape2)
setwd("./replication_timing/mouse")
all_freq <- read.table(file = "mean_germRT_recentEnh_align_hg19_mimMatch0.6_kmeans18_cons2sp"
                       , header = T, stringsAsFactors = F, sep = '\t')

all_freq_melt <- melt(all_freq[,c("kmeans.cluster", "enhrate_recent_align"
                                  , "enhrate_recent_dontalign", "mean_rt")]
                      , c("kmeans.cluster", "mean_rt"))


pdf("meanGermRT_vs_enhRate_by_kmeansClust18_recent.alignhg19_cons2sp.pdf")
ggplot(all_freq_melt, aes(x = mean_rt, y = value, color=variable, group=variable, fill=variable)) +
  geom_point()  + theme_classic() + geom_smooth(method=lm)
dev.off()

rsq <- function (x, y) cor(x, y) ^ 2

rsq(all_freq$mean_rt, all_freq$enhrate_recent_align)
#  0.9539876
rsq(all_freq$mean_rt, all_freq$enhrate_recent_dontalign)
#[1] 0.8183637

x <- cor.test(all_freq$mean_rt, all_freq$enhrate_recent_align)
round((x$estimate)**2, 2)
# cor
# 0.95
(x$p.value)
# 4.028771e-12


x <- cor.test(all_freq$mean_rt, all_freq$enhrate_recent_dontalign)
round((x$estimate)**2, 2)
# cor
# 0.82
(x$p.value)
# [1] 2.54115e-07
