library(GenomicRanges)
library(GenomicRanges)
library(ggplot2)
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")

setwd("/g/data/zk16/cc3704/replication_timing/human/breast_cancer")

b.H3K27ac <- read.table(file = "./breast_epithelium/ENCFF569TYE_H3K27ac.bed"
                        , header =, stringsAsFactors = F, sep = "\t")
b.H3K4me1 <- read.table(file = "./breast_epithelium/ENCFF440IVB_H3K4me1.bed"
                        , header =, stringsAsFactors = F, sep = "\t")
b.H3K4me3 <- read.table(file = "./breast_epithelium/ENCFF739KZZ_H3K4me3.bed"
                        , header =, stringsAsFactors = F, sep = "\t")

mcf7.H3K27ac <- read.table(file = "./MCF7/ENCFF132IWN_H3K27ac.bed", header =, stringsAsFactors = F
                           , sep = "\t")
mcf7.H3K4me1 <- read.table(file = "./MCF7/ENCFF674BKS_H3K4me1.bed", header =, stringsAsFactors = F
                           , sep = "\t")
mcf7.H3K4me3 <- read.table(file = "./MCF7/ENCFF727UPU_H3K4me3.bed", header =, stringsAsFactors = F
                           , sep = "\t")
mcf7.rt <- read.table(file = "./MCF7/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.bedGraph", header =, stringsAsFactors = F
                      , sep = "\t")

mcf7.rt$rtbin <- as.factor( cut2(mcf7.rt$V4, g=5))
idkey <- data.frame( bin=names(table(mcf7.rt$rtbin)), quintile=c(1,2,3,4,5))
# bin quintile
# 1 [-3.95,24.0)        1
# 2 [23.98,38.4)        2
# 3 [38.37,53.2)        3
# 4 [53.18,66.1)        4
# 5 [66.07,86.4]        5

mcf7.rt <- merge(mcf7.rt,  idkey ,by.x='rtbin', by.y='bin', all.x=TRUE)
colnames(mcf7.rt) <- c('rtbin','chr_rt','start_rt','end_rt','rt','rt_quintile')
mcf7.rt_gr <- with(mcf7.rt, GRanges( chr_rt , IRanges( start_rt, end_rt )))

b.enh <- rbind(b.H3K27ac, b.H3K4me1)
mcf7.enh <- rbind(mcf7.H3K27ac, mcf7.H3K4me1)

b.enh_gr <- with(b.enh, GRanges(V1, IRanges(V2+1, V3)))
mcf7.enh_gr <- with(mcf7.enh, GRanges(V1, IRanges(V2+1, V3)))

#removing proximal regions (+/- 1kb from the TSS)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
hg19_txdb <- makeTxDbFromGFF(file="/g/data/zk16/cc3704/human_data/hg19.ensGene.gtf"
                             , format="gtf",organism="Homo sapiens")
hg19_prom <- GenomicFeatures::promoters(x = hg19_txdb
                                        , upstream = 1000
                                        , downstream = 1000)

x <- as.data.frame(findOverlaps(b.enh_gr, hg19_prom)) 
b.enh <- b.enh[-unique(x$queryHits),]
x <- as.data.frame(findOverlaps(mcf7.enh_gr, hg19_prom))
mcf7.enh <- mcf7.enh[-unique(x$queryHits),]

#updating gr objects
b.enh_gr <- with(b.enh, GRanges(V1, IRanges(V2+1, V3)))
mcf7.enh_gr <- with(mcf7.enh, GRanges(V1, IRanges(V2+1, V3)))

dim(unique(mcf7.enh[,1:3]))

gains.df <- as.data.frame(gains)
losses.df <- as.data.frame(losses)
conserved.df <- as.data.frame(conserved)

#back to 0 based
gains.df$start <- gains.df$start -1
losses.df$start <- losses.df$start -1 
conserved.df$start <- conserved.df$start -1

write.table(x = gains.df, file = "breast_cancer_gains.bed", sep = '\t', quote = F)
write.table(x = losses.df, file = "breast_cancer_loses.bed", sep = '\t', quote = F)
write.table(x = conserved.df, file = "breast_cancer_unchanged.bed", sep = '\t', quote = F)

gains.df <- mcf7.enh[-unique(x$queryHits),]
losses.df <- b.enh[-unique(x$subjectHits),]
cons.df <- b.enh[unique(x$subjectHits),]

#OVERLAPPING ENHANCERS WITH MCF7 REPLICATION TIME

x <- as.data.frame(findOverlaps(gains, mcf7.rt_gr))
gains_with_rt <- cbind(gains.df[x$queryHits,]
                       , mcf7.rt[x$subjectHits,c("rt", "rt_quintile")])
x <- as.data.frame(findOverlaps(losses, mcf7.rt_gr))
losses_with_rt <- cbind(losses.df[x$queryHits,]
                        , mcf7.rt[x$subjectHits,c("rt", "rt_quintile")])
x <- as.data.frame(findOverlaps(conserved, mcf7.rt_gr))
cons_with_rt <- cbind(cons.df[x$queryHits,]
                      , mcf7.rt[x$subjectHits,c("rt", "rt_quintile")])

gains_with_rt$study <- "Gains"
losses_with_rt$study <- "Losses"
cons_with_rt$study <- "Conserved"

write.table(x = gains_with_rt, file = "gains_enh_with_MCF7_rt.txt"
            , row.names = F, sep = '\t', quote = F)
write.table(x = losses_with_rt, file = "losses_enh_with_MCF7_rt.txt"
            , row.names = F, sep = '\t', quote = F)
write.table(x = cons_with_rt, file = "cons_enh_with_MCF7_rt.txt"
            , row.names = F, sep = '\t', quote = F)

# remove RT quintile to aggregate rt to unique enhancers
gains_with_rt$rt_quintile <- NULL
losses_with_rt$rt_quintile <- NULL
cons_with_rt$rt_quintile <- NULL 

# mean rt for unique enhancers
gains_with_rt <- aggregate(rt~., gains_with_rt, mean)
losses_with_rt <- aggregate(rt~., losses_with_rt, mean)
cons_with_rt <- aggregate(rt~., cons_with_rt, mean)

df <- rbind(gains_with_rt, losses_with_rt, cons_with_rt)
df$study <- factor(df$study, levels = c("Gains", "Conserved", "Losses"))

pdf("breast_epith_vs_mcf7_enh_mcf7RT_final.pdf")
ggplot(df, aes(x=study, y=log10(rt))) + 
  geom_boxplot() + theme_classic()
dev.off()

#without log10
pdf("breast_epith_vs_mcf7_enh_mcf7RT_final.0.pdf")
ggplot(df, aes(x=study, y=(rt))) + 
  geom_boxplot() + theme_classic()
dev.off()

wilcox.test(gains_with_rt$rt, cons_with_rt$rt)
wilcox.test(losses_with_rt$rt, cons_with_rt$rt)
wilcox.test(gains_with_rt$rt, losses_with_rt$rt)

