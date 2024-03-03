library("Hmisc")
library("GenomicRanges")

gains1 <- read.table(file = "LNCaP_enh_gains.bed", header = F
                     , stringsAsFactors = F, sep = '\t')
losses1 <- read.table(file = "PReC_enh_losses.bed", header = F
                      , stringsAsFactors = F, sep = '\t')
cons1 <- read.table(file = "PReC_enh_conserved.bed", header = F
                    , stringsAsFactors = F, sep = '\t')

ln_rt <- read.delim('GSE98730_LNCaP_WA.bed', header=F)
ln_rt$rtbin <- as.factor( cut2(ln_rt$V4, g = 5))
idkey <- data.frame( bin = names(table(ln_rt$rtbin)), quintile = c(1,2,3,4,5))

ln_rt <- merge(ln_rt,  idkey ,by.x = 'rtbin', by.y = 'bin', all.x = TRUE)
names(ln_rt) <- c('rtbin','chr_rt','start_rt','end_rt','ln_rt','rt_quintile')
ln_rt_gr <- with(ln_rt, GRanges( chr_rt , IRanges( start_rt, end_rt )))

gains1 <- unique(gains1[1:3])
losses1 <- unique(losses1[1:3])
cons1 <- unique(cons1[1:3])

gains1_gr <- with(gains1, GRanges(V1 , IRanges(V2, V3)))
losses1_gr <- with(losses1, GRanges(V1 , IRanges(V2, V3)))
cons1_gr <- with(cons1, GRanges(V1 , IRanges(V2, V3)))

x <- as.data.frame(findOverlaps(gains1_gr, ln_rt_gr))
gains1_with_rt <- cbind(gains1[x$queryHits,], ln_rt[x$subjectHits,"ln_rt", drop =F])

x <- as.data.frame(findOverlaps(losses1_gr, ln_rt_gr))
losses1_with_rt <- cbind(losses1[x$queryHits,], ln_rt[x$subjectHits,"ln_rt", drop =F])

x <- as.data.frame(findOverlaps(cons1_gr, ln_rt_gr))
cons1_with_rt <- cbind(cons1[x$queryHits,], ln_rt[x$subjectHits,"ln_rt", drop =F])

gains1_with_rt <- aggregate(ln_rt~., gains1_with_rt, mean)
losses1_with_rt <- aggregate(ln_rt~., losses1_with_rt, mean)
cons1_with_rt <- aggregate(ln_rt~., cons1_with_rt, mean)

gains1_with_rt$study <- "Gains"
losses1_with_rt$study <- "Losses"
cons1_with_rt$study <- "Conserved"

df <- rbind(gains1_with_rt[,c( "V1", "V2", "V3", "ln_rt", "study")]
            , losses1_with_rt[,c( "V1", "V2", "V3", "ln_rt", "study")]
            , cons1_with_rt[,c( "V1", "V2", "V3", "ln_rt", "study")])
df$study <- factor(df$study, levels = c("Gains", "Conserved", "Losses"))

pdf("PrEC_vs_LNCaP_enh1_LNCaP_RT.pdf")
ggplot(df, aes(x=study, y=log10(ln_rt))) +
  geom_boxplot() + theme_classic()
dev.off()

pdf("PrEC_vs_LNCaP_enh1_LNCaP_RT.0.pdf")
ggplot(df, aes(x=study, y=(ln_rt))) +
  geom_boxplot() + theme_classic()
dev.off()

x <- wilcox.test(gains1_with_rt$ln_rt,cons1_with_rt$ln_rt)
# Wilcoxon rank sum test with continuity correction
# 
# data:  gains1_with_rt$ln_rt and cons1_with_rt$ln_rt
# W = 2.324e+09, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
x$p.value#0

x <- wilcox.test(losses1_with_rt$ln_rt,cons1_with_rt$ln_rt)
# Wilcoxon rank sum test with continuity correction
# 
# data:  losses1_with_rt$ln_rt and cons1_with_rt$ln_rt
# W = 1284499236, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
# 
x$p.value
# [1] 1.519255e-104

x <- wilcox.test(gains1_with_rt$ln_rt,losses1_with_rt$ln_rt)
# Wilcoxon rank sum test with continuity correction
# 
# data:  gains1_with_rt$ln_rt and losses1_with_rt$ln_rt
# W = 4405603688, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
x$p.value
# [1] 0
