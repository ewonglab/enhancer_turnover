
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library("Formula")
library("Hmisc")
library(reshape2)
library(ggplot2)

gains <- fread("gains_thyroid_union_1kb.bed", header = F)
loses <- fread("loses_thyroid_union_1kb.bed", header = F)
unch <- fread("unchanged_thyroid_union_1kb.bed", header = F)

thyroid.rt <- fread(file = "../replicon_thyroid_hg19.bed", header =F)
thyroid.rt$rtbin <- as.factor( cut2(thyroid.rt$V4, g=5))
idkey <- data.frame( bin=names(table(thyroid.rt$rtbin)), quintile=c(1,2,3,4,5))

thyroid.rt <- merge(thyroid.rt,  idkey ,by.x='rtbin', by.y='bin', all.x=TRUE)
colnames(thyroid.rt) <- c('rtbin','chr_rt','start_rt','end_rt','rt','rt_quintile')
thyroid.rt_gr <- with(thyroid.rt, GRanges( chr_rt , IRanges( start_rt+1, end_rt )))

#make gr objects for enhancers
gains_gr <- with(gains, GRanges( V1 , IRanges( V2+1, V3 )))
loses_gr <- with(loses, GRanges( V1 , IRanges( V2+1, V3 )))
unch_gr <- with(unch, GRanges( V1 , IRanges( V2+1, V3 )))

# OVERLAPPING ENHANCERS WITH THYROID REPLICATION TIME QUINTILES

x <- as.data.frame(findOverlaps(gains_gr, thyroid.rt_gr))
gains_with_rt <- cbind(gains[x$queryHits,]
                       , thyroid.rt[x$subjectHits,c("rt", "rt_quintile")])
x <- as.data.frame(findOverlaps(loses_gr, thyroid.rt_gr))
losses_with_rt <- cbind(loses[x$queryHits,]
                        , thyroid.rt[x$subjectHits,c("rt", "rt_quintile")])
x <- as.data.frame(findOverlaps(unch_gr, thyroid.rt_gr))
unch_with_rt <- cbind(unch[x$queryHits,]
                      , thyroid.rt[x$subjectHits,c("rt", "rt_quintile")])

gains_with_rt$study <- "Gains"
losses_with_rt$study <- "Losses"
unch_with_rt$study <- "Conserved"

write.table(x = gains_with_rt, file = "gains_enh_with_thyroid_rt.txt"
            , row.names = F, sep = '\t', quote = F)
write.table(x = losses_with_rt, file = "losses_enh_with_thyroid_rt.txt"
            , row.names = F, sep = '\t', quote = F)
write.table(x = unch_with_rt, file = "unch_enh_with_thyroid_rt.txt"
            , row.names = F, sep = '\t', quote = F)


### remove rt and plot number of enhancers by rt quintile
gains_with_rt$rt <- NULL
losses_with_rt$rt <- NULL
unch_with_rt$rt <- NULL

### remove rt and plot number of enhancers by rt quintile
gains_with_rt <- unique(gains_with_rt)#41056     5
losses_with_rt <- unique(losses_with_rt)#34234     5
unch_with_rt <- unique(unch_with_rt)#28431     5


df <- data.frame(gains_enh = table(gains_with_rt$rt_quintile)
                 , losses_enh = table(losses_with_rt$rt_quintile)
                 , cons_enh = table(unch_with_rt$rt_quintile))

df
# gains_enh.Var1 gains_enh.Freq losses_enh.Var1 losses_enh.Freq cons_enh.Var1
# 1              1           2519               1            2044             1
# 2              2           4544               2            3866             2
# 3              3           7992               3            6826             3
# 4              4          11494               4            9738             4
# 5              5          14507               5           11760             5
# cons_enh.Freq
# 1          1091
# 2          2562
# 3          5184
# 4          8134
# 5         11460

df$gains_enh_prop <- with(df, (gains_enh.Freq)/(gains_enh.Freq+cons_enh.Freq))
df$losses_enh_prop <- with(df, (losses_enh.Freq)/(losses_enh.Freq+cons_enh.Freq))
df$losses_enh_prop <- -1*(df$losses_enh_prop)
df$gains_enh.Var1 <- factor(df$gains_enh.Var1, levels = c(1,2,3,4,5))

df <- df[,c("gains_enh.Var1", "gains_enh_prop", "losses_enh_prop")]
colnames(df)[1] <- "rt_quintile"
df <- melt(df)

pdf("thyroid_enh_proportions_along_RT.pdf")
ggplot(df, aes(fill=variable, y=value, x=rt_quintile)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() + ylim(-1, 1)
dev.off()

# fisher test
tmp <- data.frame(gains_enh = table(gains_with_rt$rt_quintile)
                  , losses_enh = table(losses_with_rt$rt_quintile)
                  , cons_enh = table(unch_with_rt$rt_quintile))

#gains vs unchanged
for(i in 1:4){
  print(paste("quintile", i, "vs quintile", (i+1)))
  m <- matrix(data = c(tmp[tmp$gains_enh.Var1 ==i, "gains_enh.Freq"]
                       , tmp[tmp$gains_enh.Var1 ==i, "cons_enh.Freq"]
                       , tmp[tmp$gains_enh.Var1 ==(i+1), "gains_enh.Freq"]
                       , tmp[tmp$gains_enh.Var1 ==(i+1), "cons_enh.Freq"])
              , nrow = 2, ncol = 2, byrow = F)
  x <- fisher.test(m, alternative = "greater")
  print(paste("p=", x$p.value))
  print(paste("estimate=", x$estimate))
}

# [1] "quintile 1 vs quintile 2" # must be the p-values of Fig6E but they are not shown in the figure
# [1] "p= 8.28113904175537e-10"  
# [1] "estimate= 1.30176419208928"
# [1] "quintile 2 vs quintile 3"
# [1] "p= 2.17608141781746e-06"
# [1] "estimate= 1.15044346483633"
# [1] "quintile 3 vs quintile 4"
# [1] "p= 7.80243327053553e-05"
# [1] "estimate= 1.09099338010851"
# [1] "quintile 4 vs quintile 5"
# [1] "p= 4.66251530561649e-09"
# [1] "estimate= 1.11627841531644"

#loses vs unchanged
for(i in 1:4){
  print(paste("quintile", i, "vs quintile", (i+1)))
  m <- matrix(data = c(tmp[tmp$gains_enh.Var1 ==i, "losses_enh.Freq"]
                       , tmp[tmp$gains_enh.Var1 ==i, "cons_enh.Freq"]
                       , tmp[tmp$gains_enh.Var1 ==(i+1), "losses_enh.Freq"]
                       , tmp[tmp$gains_enh.Var1 ==(i+1), "cons_enh.Freq"])
              , nrow = 2, ncol = 2, byrow = F)
  x <- fisher.test(m, alternative = "greater")
  print(paste("p=", x$p.value))
  print(paste("estimate=", x$estimate))
}

# [1] "quintile 1 vs quintile 2"
# [1] "p= 9.23304869407604e-07"
# [1] "estimate= 1.24154522504648"
# [1] "quintile 2 vs quintile 3"
# [1] "p= 7.62743833806008e-06"
# [1] "estimate= 1.14598174760095"
# [1] "quintile 3 vs quintile 4"
# [1] "p= 3.24664265027635e-05"
# [1] "estimate= 1.09985139508034"
# [1] "quintile 4 vs quintile 5"
# [1] "p= 5.75375877516542e-15"
# [1] "estimate= 1.16664821715224"
