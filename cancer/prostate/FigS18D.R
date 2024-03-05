library(data.table)
library(GenomicRanges)
library("Formula")
library("Hmisc")
library(reshape2)

gains1 <- read.table(file = "LNCaP_enh_gains.bed", header = F, stringsAsFactors = F
                     , sep = '\t')
losses1 <- read.table(file = "PReC_enh_losses.bed", header = F, stringsAsFactors = F
                      , sep = '\t')
cons1 <- read.table(file = "PReC_enh_conserved.bed", header = F, stringsAsFactors = F
                    , sep = '\t')

pr_rt <- read.delim('GSE98730_PrEC_WA.bed', header=F)
pr_rt$rtbin <- as.factor(cut2(pr_rt$V4, g=5))
idkey <- data.frame( bin=names(table(pr_rt$rtbin)), quintile=c(1,2,3,4,5))

pr_rt <- merge(pr_rt,  idkey ,by.x='rtbin', by.y='bin', all.x=TRUE)#2784248       6
names(pr_rt) <- c('rtbin','chr_rt','start_rt','end_rt','pr_rt','rt_quintile')
pr_rt_gr <- with(pr_rt, GRanges( chr_rt , IRanges( start_rt, end_rt )))

gains1 <- unique(gains1[1:3])
losses1 <- unique(losses1[1:3])
cons1 <- unique(cons1[1:3])

gains1_gr <- with(gains1, GRanges(V1 , IRanges(V2, V3)))
losses1_gr <- with(losses1, GRanges(V1 , IRanges(V2, V3)))
cons1_gr <- with(cons1, GRanges(V1 , IRanges(V2, V3)))

x <- as.data.frame(findOverlaps(gains1_gr, pr_rt_gr))
gains1_with_rt <- cbind(gains1[x$queryHits,], pr_rt[x$subjectHits,"pr_rt", drop =F])

x <- as.data.frame(findOverlaps(losses1_gr, pr_rt_gr))
losses1_with_rt <- cbind(losses1[x$queryHits,], pr_rt[x$subjectHits,"pr_rt", drop =F])

x <- as.data.frame(findOverlaps(cons1_gr, pr_rt_gr))
cons1_with_rt <- cbind(cons1[x$queryHits,], pr_rt[x$subjectHits,"pr_rt", drop =F])

gains1_with_rt <- aggregate(pr_rt~., gains1_with_rt, mean)
losses1_with_rt <- aggregate(pr_rt~., losses1_with_rt, mean)
cons1_with_rt <- aggregate(pr_rt~., cons1_with_rt, mean)

gains1_with_rt$study <- "Gains"
losses1_with_rt$study <- "Losses"
cons1_with_rt$study <- "Conserved"

df <- rbind(gains1_with_rt[,c( "V1", "V2", "V3", "pr_rt", "study")]
            , losses1_with_rt[,c( "V1", "V2", "V3", "pr_rt", "study")]
            , cons1_with_rt[,c( "V1", "V2", "V3", "pr_rt", "study")])
df$study <- factor(df$study, levels = c("Gains", "Conserved", "Losses"))

table(df$study)

# Gains Conserved    Losses 
# 159204     40658     68566 

pdf("PrEC_vs_LNCaP_enh1_PrER_RT.pdf")
ggplot(df, aes(x=study, y=log10(pr_rt))) + 
  geom_boxplot() + theme_classic()
dev.off()

pdf("PrEC_vs_LNCaP_enh1_PrER_RT.0.pdf")
ggplot(df, aes(x=study, y=(pr_rt))) + 
  geom_boxplot() + theme_classic()
dev.off()


x <- wilcox.test(cons1_with_rt$pr_rt, gains1_with_rt$pr_rt, alternative = "greater")
# Wilcoxon rank sum test with continuity correction
# 
# data:  cons1_with_rt$pr_rt and gains1_with_rt$pr_rt
# W = 4656265216, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0
x$p.value
# [1] 0

x <- wilcox.test(cons1_with_rt$pr_rt, losses1_with_rt$pr_rt, alternative = "greater")
# Wilcoxon rank sum test with continuity correction
# 
# data:  cons1_with_rt$pr_rt and losses1_with_rt$pr_rt
# W = 1446836330, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0
x$p.value
# [1] 3.755266e-26

# TWO-SIDED
x <- wilcox.test(cons1_with_rt$pr_rt, gains1_with_rt$pr_rt)
x
# Wilcoxon rank sum test with continuity correction
# 
# data:  cons1_with_rt$pr_rt and gains1_with_rt$pr_rt
# W = 4656265216, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

x$p.value # 0

x <- wilcox.test(cons1_with_rt$pr_rt, losses1_with_rt$pr_rt)
# Wilcoxon rank sum test with continuity correction
# 
# data:  cons1_with_rt$pr_rt and losses1_with_rt$pr_rt
# W = 1446836330, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

x$p.value
# [1] 7.510533e-26

