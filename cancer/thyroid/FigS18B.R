library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library("Formula")
library("Hmisc")
library(reshape2)
library(ggplot2)

gains_with_rt <- read.table(file = "gains_enh_with_thyroid_rt.txt"
                            , header=T, sep = '\t', stringsAsFactors = F)
losses_with_rt <- read.table(file = "losses_enh_with_thyroid_rt.txt"
                             , header=T, sep = '\t', stringsAsFactors = F)
unch_with_rt <- read.table(file = "unch_enh_with_thyroid_rt.txt"
                           , header=T, sep = '\t', stringsAsFactors = F)

df <- rbind(gains_with_rt, losses_with_rt, unch_with_rt)
df$study <- factor(df$study, levels = c("Gains", "Conserved", "Losses"))

df <- unique(df[,c("V1", "V2", "V3", "rt", "study")])
df <- aggregate(rt~., df, mean) 

pdf("thyroid_union_1kb_RT.pdf")
ggplot(df, aes(x=study, y=rt)) +
  geom_boxplot() + theme_classic()
dev.off()

x <- wilcox.test(df[df$study=="Gains", "rt"], df[df$study=="Conserved", "rt"])
x$p.value
# [1] 3.764045e-70
x <- wilcox.test(df[df$study=="Losses", "rt"], df[df$study=="Conserved", "rt"])
x$p.value
# [1] 2.127589e-86

table(df$study)
# Gains Conserved    Losses 
# 41046     28426     34228 
