library(data.table)
library(BiocGenerics, lib = "/srv/scratch/z5291980/extra_R_packages")
library(S4Vectors, lib = "/srv/scratch/z5291980/extra_R_packages")

library(IRanges, lib = "/srv/scratch/z5291980/extra_R_packages")
library(GenomeInfoDb, lib = "/srv/scratch/z5291980/extra_R_packages")
library(GenomicRanges, lib = "/srv/scratch/z5291980/extra_R_packages")
require(dplyr)
library(ggplot2)

library(mixtools, lib = "/srv/scratch/z5291980/extra_R_packages/R_4.0.0")
library(plyr); library(dplyr)
library(tidyr)

#H9 hESC eplication time
rt <- read.delim('/srv/scratch/wonglab/data/repliseq/RT_H9_ESC_Ext29405702_hg19.bedgraph'
              ,skip=11,header=F)
rt_gr <- with(rt, GRanges( V1 , IRanges( V2+1, V3 )))
names(rt) <- c('chr','start','end','rt')

#reading coordinates of every TF and overlap with H9 replication time
dfplot <- data.frame( ID=character(), rt=numeric() )
filenames <- list.files("/srv/scratch/wonglab/data/K562_TF/hg19/"
                        , pattern="*\\.bed", full.names=TRUE)#71 in total
set.seed(1)
for(i in filenames) {
  obj <- read.delim(i, header=F)
  obj$ID <- paste0('ID_',seq(1:nrow(obj)))
  gr <- with(obj, GRanges( V1 , IRanges( V2+1, V3 )))
  x <- as.data.frame(findOverlaps( rt_gr, gr))
  d <- cbind(rt[x$queryHits,], obj[x$subjectHits,])
  ave <- as.data.frame(group_by(d, ID) %>% summarise(rt = mean(rt)))
  ave$ID <- strsplit(basename(i),'\\.')[[1]][1]
  dfplot <- rbind(dfplot, ave)}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)}

mt <- data.frame(ID=character(), mu1=numeric(), mu2=numeric()
                  , sd1=numeric(), sd2=numeric(),  l1=numeric()
                  , l2=numeric(),  max1=numeric(), max2=numeric())

setwd("/srv/scratch/z5291980/replication_timing/mutations")
set.seed(1)
for (i in filenames) {
  id <- strsplit(basename(i),'\\.')[[1]][1]
  obj <- read.delim(i, header=F)
  obj$ID <- paste0('ID_',seq(1:nrow(obj)))
  gr <- with(obj, GRanges( V1 , IRanges( V2+1, V3 )))
  x <- as.data.frame(findOverlaps( rt_gr, gr))
  d <- cbind(rt[x$queryHits,], obj[x$subjectHits,])
  #we iterate over the TFs (chip-seq sets)
  sub <- d[,c("ID", "rt")]
  sub <- aggregate(rt~ID, sub, mean)
  set.seed(1)
  mixmdl <- normalmixEM(x=sub$rt, k = 2) #this function returns the best mu, sigma and lambda 
  #parameters for a normal mixed distribution of k components
  pdf(sprintf('%s_mixtools.k562_RT_n2.pdf', id))
  print (data.frame(x = mixmdl$x) %>%
           ggplot() +
           geom_histogram(aes(x, ..density..), binwidth = .2, colour = "black", 
                          fill = "white") +
           stat_function(geom = "line", fun = plot_mix_comps,
                         args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                         colour = "red", lwd = 1.5) +
           stat_function(geom = "line", fun = plot_mix_comps,
                         args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                         colour = "blue", lwd = 1.5) +
           ylab("Density") + theme_minimal())
  dev.off()
  max1 <- max(mixmdl$lambda[1] * dnorm(mixmdl$x, mixmdl$mu[1],  mixmdl$sigma[1]))
  max2 <- max(mixmdl$lambda[2] * dnorm(mixmdl$x, mixmdl$mu[2],  mixmdl$sigma[2]))
  mt <- rbind( mt, c(id , mixmdl$mu[1], mixmdl$mu[2] ,  mixmdl$sigma[1], mixmdl$sigma[2]
                      , mixmdl$lambda[1], mixmdl$lambda[2],  max1, max2))
}

names(mt) <- c('ID','mu1', 'mu2',  'sd1','sd2','l1','l2','mode1', 'mode2')

comp <- read.delim('/srv/scratch/z5291980/jaspar_data/JASPAR_core_nucl_comp_and_ic.txt')
comp <- comp[grep("Homo", comp$Species),]
comp <-comp[!duplicated(comp[ , c("motif.name")]),]#639  16

mt <- mt %>% separate(ID, c("Acc", "ID"))

write.table(mt, 'ESC_RT_with_K562_TF.mixture.txt',sep='\t', row.names=F, quote=F)

mt <- read.delim('ESC_RT_with_K562_TF.mixture.txt')
m <- merge( mt, comp, by.x='ID', by.y='motif.name')
m <- m[order(m$AT.proportion),]
m <- m[!duplicated(m[ , c("ID")]),]#71 25

sub <- subset(m, mu1 < mu2)

round(with(sub, cor(mu2, AT.proportion)), 2)#-0.11
round(with(sub, cor(mu2, GC.proportion)), 2)#0.12
round(with(sub, cor(mu1, AT.proportion)), 2)#-0.2
round(with(sub, cor(mu1, GC.proportion)), 2)#0.24

library(cowplot)
library("ggpubr", lib = '/srv/scratch/z5291980/extra_R_packages/R_4.0.0')

pdf('ESC_RT_K562_earlylate.proportion.2.pdf')  
p1 <- ggplot(sub, aes(mu2, AT.proportion)) + 
  geom_point() +  stat_smooth(method = "lm") + theme_classic()
p2 <- ggplot(sub, aes(mu2, GC.proportion)) + 
  geom_point() +  stat_smooth(method = "lm") + theme_classic()
p3 <- ggplot(sub, aes(mu1, AT.proportion)) + 
  geom_point() +  stat_smooth(method = "lm") + theme_classic()
p4 <- ggplot(sub, aes(mu1, GC.proportion)) + 
  geom_point() +  stat_smooth(method = "lm") + theme_classic()
plot_grid(p1, p2 , p3, p4,
          labels = c('AT early', 'GC early', 'AT late', 'GC late'),
          ncol = 2, nrow = 2)
dev.off() 
