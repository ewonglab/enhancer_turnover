library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library(data.table)
library(reshape2)

#replication time of union gains, loses and unchanged enhancers

gains_with_rt <- read.table(file = "gains_enh_with_pHSC_rt.txt"
                            , header=T, sep = '\t', stringsAsFactors = F)
losses_with_rt <- read.table(file = "losses_enh_with_pHSC_rt.txt"
                             , header=T, sep = '\t', stringsAsFactors = F)
unch_with_rt <- read.table(file = "unch_enh_with_pHSC_rt.txt"
                           , header=T, sep = '\t', stringsAsFactors = F)

gains_with_rt$rt_quintile <- NULL
losses_with_rt$rt_quintile <- NULL
unch_with_rt$rt_quintile <- NULL 

gains_with_rt <- aggregate(rt~., gains_with_rt, mean)
losses_with_rt <- aggregate(rt~., losses_with_rt, mean)
unch_with_rt <- aggregate(rt~., unch_with_rt, mean)

df <- rbind(gains_with_rt, losses_with_rt, unch_with_rt)
df$study <- factor(df$study, levels = c("Gains", "Conserved", "Losses"))

pdf("AML_union_1kb_pHSC_RT.pdf")
ggplot(df, aes(x=study, y=rt)) +
  geom_boxplot() + theme_classic()
dev.off()

x <- wilcox.test(df[df$study=="Gains", "rt"], df[df$study=="Conserved", "rt"])
x$p.value #2.747729e-109
x <- wilcox.test(df[df$study=="Losses", "rt"], df[df$study=="Conserved", "rt"])
x$p.value #4.205047e-205

