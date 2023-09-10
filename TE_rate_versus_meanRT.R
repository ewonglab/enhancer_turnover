# Correlation between the rate of TE turnover and mean germline DNA replication time across mouse DNA replication time k-means clusters (k = 18 clusters)
# TE rate was calculated as log (number of ancestral TE / number of mouse-specific TE)

library(data.table)
library(GenomicRanges)
library(ggplot2)

mouse_rt <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header =T
                       , stringsAsFactors = F, sep = '\t')
mouse_rt <- mouse_rt[complete.cases(mouse_rt),]

# add bin coordinates
mouse_rt$chr <- sub(":.*", "", rownames(mouse_rt))
mouse_rt$start <- sub("_.*", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$end <- sub(".*_", "", sub(".*:", "", rownames(mouse_rt)))
mouse_rt$start <- as.integer(mouse_rt$start)
mouse_rt$end <- as.integer(mouse_rt$end)
mouse_rt$peakid <- rownames(mouse_rt)

# read TE
mouseonly_te <- read.table(file = "mouse_repeatmasker_mouse_only", header = T, stringsAsFactors = F, sep ='\t')
mouseanc_te <- read.table(file = "mouse_repeatmasker_ancestral", header = T, stringsAsFactors = F, sep ='\t')

# Calculate mean germ line DNA replication time 
mouse_rt$mean_rt <- apply(mouse_rt[,c("PGC.female.1", "PGC.male.1", "sperm.1", "sperm.2")]
                          , 1, mean)
mouse_rt$bin_id <- rownames(mouse_rt)

# GenomicRanges objects
mouse_rt.gr <- with(mouse_rt, GRanges(chr , IRanges( start+1, end )))

mouseonly_te.gr <- with(mouseonly_te, GRanges(query_sequence
                                              , IRanges( query_start+1, query_end)))
mouseanc_te.gr <- with(mouseanc_te, GRanges(query_sequence
                                            , IRanges( query_start+1, query_end)))

# count number of mouseonly and ancestral repeats

# add unique ID 
mouseonly_te$repeat_id <- 1:nrow(mouseonly_te)
mouseanc_te$repeat_id <- 1:nrow(mouseanc_te)

x <- as.data.frame(findOverlaps(mouse_rt.gr, mouseonly_te.gr))
rt_with_mouseonly <- cbind(mouse_rt[x$queryHits,"kmeans.cluster", drop=F]
                           , mouseonly_te[x$subjectHits,"repeat_id", drop=F])# 414765      2


x <- as.data.frame(findOverlaps(mouse_rt.gr, mouseanc_te.gr))
rt_with_anc <- cbind(mouse_rt[x$queryHits, "kmeans.cluster", drop=F]
                     , mouseanc_te[x$subjectHits, "repeat_id", drop=F])# 1595211       2

# keep unique repeat ids and cluster number so that every repeat is counted
# only once per k-means cluster

rt_with_mouseonly <- unique(rt_with_mouseonly)#414640      2
rt_with_anc <- unique(rt_with_anc)#1594778       2

# frequency of mouseonly and ancestral repeats per cluster
rt_with_mouseonly.freq <- as.data.frame(table(rt_with_mouseonly$kmeans.cluster))
rt_with_anc.freq <- as.data.frame(table(rt_with_anc$kmeans.cluster))

colnames(rt_with_mouseonly.freq) <- c("kmeans.cluster", "N_mouseonly")
colnames(rt_with_anc.freq) <- c("kmeans.cluster", "N_ancestral")

repeat_freq <- merge(rt_with_mouseonly.freq, rt_with_anc.freq)
# TE rate per cluster
repeat_freq$TE_rate <- with(repeat_freq, log(N_ancestral/N_mouseonly))

mouse_rt_kmean <- aggregate(mean_rt ~kmeans.cluster, mouse_rt, mean)

# merge TE rate and mean RT
repeat_freq <- merge(repeat_freq, mouse_rt_kmean, by = "kmeans.cluster")

pdf("mouse_TErate_vs_meanGermRT.pdf")
ggplot(repeat_freq, aes(x=mean_rt, y=TE_rate)) +
  geom_point() + geom_smooth(method=lm) + theme_classic()
dev.off()

# calculate correlation
cor_t <- cor.test(repeat_freq$mean_rt, repeat_freq$TE_rate)
cor_t$estimate
# cor
# 0.9779688
(cor_t$estimate)**2
# cor
# 0.9564231

round((cor_t$estimate)**2, 2)
# 0.96
cor_t$p.value
# [1] 2.604556e-12
