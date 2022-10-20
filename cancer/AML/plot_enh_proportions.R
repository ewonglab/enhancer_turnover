library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library(data.table)
library(reshape2)

setwd("/g/data/zk16/cc3704/replication_timing/human/AML")

gains <- fread("gains_AML_union_1kb.bed", header = F)
loses <- fread("loses_AML_union_1kb.bed", header = F)
unch <- fread("unchanged_AML_union_1kb.bed", header = F)

eryth.rt <- fread(file = "../replicon_pHSC_hg19.bed", header =F)
eryth.rt$rtbin <- as.factor( cut2(eryth.rt$V4, g=5))
idkey <- data.frame( bin=names(table(eryth.rt$rtbin)), quintile=c(1,2,3,4,5))

# bin quintile
# 1 [-4.47,-3.173)        1
# 2 [-3.17,-2.873)        2
# 3 [-2.87,-2.544)        3
# 4 [-2.54,-2.174)        4
# 5 [-2.17,-0.416]        5

eryth.rt <- merge(eryth.rt,  idkey ,by.x='rtbin', by.y='bin', all.x=TRUE)#6127658       6
colnames(eryth.rt) <- c('rtbin','chr_rt','start_rt','end_rt','rt','rt_quintile')
eryth.rt_gr <- with(eryth.rt, GRanges( chr_rt , IRanges( start_rt+1, end_rt )))

#make gr objects for enhancers
gains_gr <- with(gains, GRanges( V1 , IRanges( V2+1, V3 )))
loses_gr <- with(loses, GRanges( V1 , IRanges( V2+1, V3 )))
unch_gr <- with(unch, GRanges( V1 , IRanges( V2+1, V3 )))

# no overlap with proximal regions
hg19_txdb <- makeTxDbFromGFF(file="/g/data/zk16/cc3704/human_data/hg19.ensGene.gtf"
                             , format="gtf",organism="Homo sapiens")
hg19_prom <- GenomicFeatures::promoters(x = (hg19_txdb)
                                        , upstream = 1000
                                        , downstream = 1000)

x <- as.data.frame(findOverlaps(gains_gr, hg19_prom))
dim((x))#none
x <- as.data.frame(findOverlaps(loses_gr, hg19_prom))
dim(x)#none
x <- as.data.frame(findOverlaps(unch_gr, hg19_prom))
dim(x)#none

#OVERLAPPING ENHANCERS WITH THYROID REPLICATION TIME QUINTILES

x <- as.data.frame(findOverlaps(gains_gr, eryth.rt_gr))
gains_with_rt <- cbind(gains[x$queryHits,]
                       , eryth.rt[x$subjectHits,c("rt", "rt_quintile")])
x <- as.data.frame(findOverlaps(loses_gr, eryth.rt_gr))
losses_with_rt <- cbind(loses[x$queryHits,]
                        , eryth.rt[x$subjectHits,c("rt", "rt_quintile")])
x <- as.data.frame(findOverlaps(unch_gr, eryth.rt_gr))
unch_with_rt <- cbind(unch[x$queryHits,]
                      , eryth.rt[x$subjectHits,c("rt", "rt_quintile")])

gains_with_rt$study <- "Gains"
losses_with_rt$study <- "Losses"
unch_with_rt$study <- "Conserved"

#save enhancers with RT
### duplicated enhancers contained, aggregate RT
write.table(x = gains_with_rt, file = "gains_enh_with_pHSC_rt.txt"
            , row.names = F, sep = '\t', quote = F)
write.table(x = losses_with_rt, file = "losses_enh_with_pHSC_rt.txt"
            , row.names = F, sep = '\t', quote = F)
write.table(x = unch_with_rt, file = "unch_enh_with_pHSC_rt.txt"
            , row.names = F, sep = '\t', quote = F)

### remove rt and plot number of enhancers by rt quintile
gains_with_rt$rt <- NULL
losses_with_rt$rt <- NULL
unch_with_rt$rt <- NULL

### remove rt and plot number of enhancers by rt quintile
gains_with_rt <- unique(gains_with_rt)#58880     5
losses_with_rt <- unique(losses_with_rt)#48266     5
unch_with_rt <- unique(unch_with_rt)#45216     5

df <- data.frame(gains_enh = table(gains_with_rt$rt_quintile)
                 , losses_enh = table(losses_with_rt$rt_quintile)
                 , cons_enh = table(unch_with_rt$rt_quintile))

df
# gains_enh.Var1 gains_enh.Freq losses_enh.Var1 losses_enh.Freq cons_enh.Var1
# 1              1           4368               1            4137             1
# 2              2           6832               2            6534             2
# 3              3          12721               3           10242             3
# 4              4          16302               4           12853             4
# 5              5          18662               5           14494             5
# cons_enh.Freq
# 1          2327
# 2          4480
# 3          9022
# 4         12805
# 5         16577


df$gains_enh_prop <- with(df, (gains_enh.Freq)/(gains_enh.Freq+cons_enh.Freq))
df$losses_enh_prop <- with(df, (losses_enh.Freq)/(losses_enh.Freq+cons_enh.Freq))
df$losses_enh_prop <- -1*(df$losses_enh_prop)
df$gains_enh.Var1 <- factor(df$gains_enh.Var1, levels = c(1,2,3,4,5))

df <- df[,c("gains_enh.Var1", "gains_enh_prop", "losses_enh_prop")]
colnames(df)[1] <- "rt_quintile"
df <- melt(df)

pdf("AML_enh_proportions_along_pHSC_RT.pdf")
ggplot(df, aes(fill=variable, y=value, x=rt_quintile)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() + ylim(-1, 1)
dev.off()

# Fisher's exact test for the difference in the proportion of gains and losses compared to unchanged enhancers between RT quintiles

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
  # print(m)
  x <- fisher.test(m, alternative = "greater")
  print(paste("p=", x$p.value))
  print(paste("estimate=", x$estimate))
}

# [1] "quintile 1 vs quintile 2"
# [1] "p= 4.53066683497165e-11"
# [1] "estimate= 1.23086413172573"

# [1] "quintile 2 vs quintile 3"
# [1] "p= 0.000471464345130215"
# [1] "estimate= 1.08153841998383"

# [1] "quintile 3 vs quintile 4"
# [1] "p= 9.23397739056848e-09"
# [1] "estimate= 1.107530979863"

# [1] "quintile 4 vs quintile 5"
# [1] "p= 5.74494351390717e-15"
# [1] "estimate= 1.13085648300595"


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
# [1] "p= 4.76453060678939e-10"
# [1] "estimate= 1.21893876467906"

# [1] "quintile 2 vs quintile 3"
# [1] "p= 1.60479234133906e-25"
# [1] "estimate= 1.28473244800492"

# [1] "quintile 3 vs quintile 4"
# [1] "p= 5.93811814914877e-11"
# [1] "estimate= 1.1309803679507"

# [1] "quintile 4 vs quintile 5"
# [1] "p= 1.5926571028726e-16"
# [1] "estimate= 1.14799589309067"
