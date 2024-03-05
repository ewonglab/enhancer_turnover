library(RColorBrewer)
library(data.table)
library(tidyr)
library(reshape2)
library(ggplot2)ewer)
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library("Formula")
library("Hmisc")
library(reshape2)
library(reshape2)

aml_mut <- read.table(file = 'aml_plot_vcf_RT_olaps.txt', header = T, stringsAsFactors = F)

df <- (aml_mut[,c("region", "rt", "Peaks", "width_norm_log")])
df$Peaks <- factor(df$Peaks, levels = c("Gains", "Unchanged", "Losses"))

#add coordinates to peaks
df$region <- sub("-", ":", df$region)
df$chr <- paste0("chr", gsub(":.*", "", df$region))
df$start <- gsub("-.*", "", gsub(".*:", "", df$region))
df$end <- gsub(".*-", "", gsub(".*:", "", df$region))
df$start <- as.numeric(df$start)
df$end <- as.numeric(df$end)

#read replication time and define quintiles
rt <- fread("../replicon_pHSC_hg19.bed", header =F)
rt$rtbin <- as.factor(cut2(rt$V4, g=5))

idkey <- data.frame( bin=names(table(rt$rtbin)), quintile=c(1,2,3,4,5))
# bin quintilelot.R
# 1 [-4.47,-3.173)        1
# 2 [-3.17,-2.873)        2
# 3 [-2.87,-2.544)        3
# 4 [-2.54,-2.174)        4
# 5 [-2.17,-0.416]        5

rt <- merge(rt,  idkey ,by.x='rtbin', by.y='bin', all.x=TRUE)#6127658       6
colnames(rt) <- c('rtbin','chr_rt','start_rt','end_rt','rt','rt_quintile')
rt_gr <- with(rt, GRanges( chr_rt , IRanges( start_rt+1, end_rt )))

enh_gr <- with(df, GRanges( chr , IRanges( start+1, end )))

x <- as.data.frame(findOverlaps(enh_gr, rt_gr))
enh_with_rt <- cbind(df[x$queryHits,], rt[x$subjectHits,])

df <- unique(enh_with_rt[,c("region", "Peaks", "width_norm_log", "rt_quintile")])#1771    4
df$rt_quintile <- as.factor(df$rt_quintile)
table(df$Peaks)

# Gains    Losses Unchanged 
# 592       303       876 

df_medians <- aggregate(width_norm_log~(rt_quintile + Peaks), df, median)
colnames(df_medians) <- c("rt_quintile", "study", "median")
df_medians$rt_quintile <- factor(df_medians$rt_quintile, levels = c(1,2,3,4,5))

se <- function(x) sqrt(var(x) / length(x))

df_se <- aggregate(width_norm_log~(rt_quintile + Peaks), df, se)
colnames(df_se) <- c("rt_quintile", "study", "se")
df_se$rt_quintile <- factor(df_se$rt_quintile, levels = c(1,2,3,4,5))

df_medians <- merge(df_medians, df_se, by=c("rt_quintile", "study"))

pdf("AML_enh_median_n_by_enhWidth_pHSC_RT.1.pdf")
ggplot(df_medians, aes(x=rt_quintile, y=(median), group=study, color=study)) +
  geom_line() + theme_classic() +
  geom_errorbar(aes(ymin=(median-se), ymax=(median+se)), width=.2,
                position=position_dodge(0.05))
dev.off()

df_unique <- unique(df[,c("region", "Peaks")])
table(df_unique$Peaks)

# Gains Unchanged    Losses 
# 592       876       303 
