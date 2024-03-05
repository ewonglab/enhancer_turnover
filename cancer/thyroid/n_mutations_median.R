library(RColorBrewer)
library(data.table)
library(tidyr)
library(reshape2)
library(ggplot2)

library(GenomicRanges)
library(GenomicFeatures)
library("Formula")
library("Hmisc")
library(data.table)
library(reshape2)

mut <- read.table(file = 'plot_vcf_RT_olaps.txt', header = T, stringsAsFactors = F)

#taking unique enhancers as in the last version of the boxplots made with 
df <- (mut[,c("region", "rt", "Peaks", "width_norm_log")])
df$Peaks <- factor(df$Peaks, levels = c("Gains", "Unchanged", "Losses"))

# add coordinates to peaks
df$region <- sub("-", ":", df$region)
df$chr <- paste0("chr", gsub(":.*", "", df$region))
df$start <- gsub("-.*", "", gsub(".*:", "", df$region))
df$end <- gsub(".*-", "", gsub(".*:", "", df$region))
df$start <- as.numeric(df$start)
df$end <- as.numeric(df$end)

#read replication time and define quintiles
rt <- fread("../replicon_thyroid_hg19.bed", header =F)
#limits of quintiles
rt$rtbin <- as.factor(cut2(rt$V4, g=5))
idkey <- data.frame( bin=names(table(rt$rtbin)), quintile=c(1,2,3,4,5))

# bin quintile
# 1 [-4.91,-3.18)        1
# 2 [-3.18,-2.86)        2
# 3 [-2.86,-2.52)        3
# 4 [-2.52,-2.18)        4
# 5 [-2.18, 2.00]        5

rt <- merge(rt,  idkey ,by.x='rtbin', by.y='bin', all.x=TRUE)
colnames(rt) <- c('rtbin','chr_rt','start_rt','end_rt','rt','rt_quintile')
rt_gr <- with(rt, GRanges( chr_rt , IRanges( start_rt+1, end_rt )))

enh_gr <- with(df, GRanges( chr , IRanges( start+1, end )))

x <- as.data.frame(findOverlaps(enh_gr, rt_gr))
enh_with_rt <- cbind(df[x$queryHits,], rt[x$subjectHits,])

df <- unique(enh_with_rt[,c("region", "Peaks", "width_norm_log", "rt_quintile")])
df$rt_quintile <- as.factor(df$rt_quintile)

#plot median values
df_medians <- aggregate(width_norm_log~(rt_quintile + Peaks), df, median)
colnames(df_medians) <- c("rt_quintile", "study", "median")
df_medians$rt_quintile <- factor(df_medians$rt_quintile, levels = c(1,2,3,4,5))

se <- function(x) sqrt(var(x) / length(x))
df_se <- aggregate(width_norm_log~(rt_quintile + Peaks), df, se)
colnames(df_se) <- c("rt_quintile", "study", "se")
df_se$rt_quintile <- factor(df_se$rt_quintile, levels = c(1,2,3,4,5))

df_medians <- merge(df_medians, df_se, by=c("rt_quintile", "study"))


pdf("thyroid_enh_median_n_by_enhWidth_thy.replicon_RT.1.pdf")
ggplot(df_medians, aes(x=rt_quintile, y=(median), group=study, color=study)) +
  geom_line() + theme_classic() +
  geom_errorbar(aes(ymin=(median-se), ymax=(median+se)), width=.2,
                position=position_dodge(0.05))
dev.off()
