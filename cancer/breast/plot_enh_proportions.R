library(GenomicRanges)
library(GenomicRanges)
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library(ggplot2)
library(reshape2)

gains_with_rt <- read.table(file = "gains_enh_with_MCF7_rt.txt"
                            , header =T, sep = '\t', stringsAsFactors = F)#221183     13
losses_with_rt <- read.table(file = "losses_enh_with_MCF7_rt.txt"
                             , header =T, sep = '\t', stringsAsFactors = F)#70710    13
cons_with_rt <- read.table(file = "cons_enh_with_MCF7_rt.txt"
                           , header =T, sep = '\t', stringsAsFactors = F)#78625    13

gains_with_rt <- unique(gains_with_rt[,c(1:3,12)])
losses_with_rt <- unique(losses_with_rt[,c(1:3,12)])
cons_with_rt <- unique(cons_with_rt[,c(1:3,12)])

df <- data.frame(gains_enh = table(gains_with_rt$rt_quintile)
                 , losses_enh = table(losses_with_rt$rt_quintile)
                 , cons_enh = table(cons_with_rt$rt_quintile))

df$gains_enh_prop <- with(df, (gains_enh.Freq)/(gains_enh.Freq+cons_enh.Freq))
df$losses_enh_prop <- with(df, (losses_enh.Freq)/(losses_enh.Freq+cons_enh.Freq))
df$losses_enh_prop <- -1*(df$losses_enh_prop)
df$gains_enh.Var1 <- factor(df$gains_enh.Var1, levels = c(1,2,3,4,5))

library(reshape2)
df <- df[,c("gains_enh.Var1", "gains_enh_prop", "losses_enh_prop")]
colnames(df)[1] <- "rt_quintile" 
df <- melt(df)

pdf("breast_epith_vs_mcf7_enh_mcf7_RT_prop_final.pdf")
ggplot(df, aes(fill=variable, y=value, x=rt_quintile)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_classic() + ylim(-1, 1)
dev.off()

#Fisher's exact tests for the difference in proportion of gains and losses between replication time quintiles

breast.c <- data.frame(quintile=1:5

                               , gains=c(3030, 12867, 25224, 35794, 43112)

                               , losses=c(2119, 4397, 7414, 10159, 15075)

                               , unchanged=c(281, 1550, 4775, 10029, 19820))



Fisher.rt.quintiles <- function(x, test_enh){
  for(i in 1:4){
    m <- as.matrix(x[i:(i+1), c(test_enh, "unchanged")])
    out <- fisher.test(m, alternative = "greater")
    print(paste("test:", paste0("Q", i), "vs", paste0("Q", i+1)))
    print(paste("odds.ratio", out$estimate))
    print(paste("p.value", out$p.value))}}


#gains

Fisher.rt.quintiles(breast.c, "gains")

# [1] "test: Q1 vs Q2"
# [1] "odds.ratio 1.29891281520728"
# [1] "p.value 4.65522324644107e-05"

# [1] "test: Q2 vs Q3"
# [1] "odds.ratio 1.57144956965454"
# [1] "p.value 2.23078437700231e-50"

# [1] "test: Q3 vs Q4"
# [1] "odds.ratio 1.48007643729781"
# [1] "p.value 2.07445868424655e-93"

# [1] "test: Q4 vs Q5"
# [1] "odds.ratio 1.64078377797833"
# [1] "p.value 9.92209944612851e-274"


#losses
Fisher.rt.quintiles(breast.c, "losses")

# [1] "test: Q1 vs Q2"
# [1] "odds.ratio 2.65794127658506"
# [1] "p.value 4.46423552266495e-51"

# [1] "test: Q2 vs Q3"
# [1] "odds.ratio 1.82695296296269"
# [1] "p.value 1.0995946085427e-69"

# [1] "test: Q3 vs Q4"
# [1] "odds.ratio 1.53275235005619"
# [1] "p.value 4.73862708875979e-76"

# [1] "test: Q4 vs Q5"
# [1] "odds.ratio 1.33177511856019"
# [1] "p.value 6.14042674340386e-59"
