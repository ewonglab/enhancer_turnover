library(GenomicRanges)
require(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)

color.palette <- function(steps, n.steps.between=NULL, ...){
 
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
 
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
 
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
 
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}

setwd("./replication_timing/mouse")

gilbert.sc <- read.table(file = "hiratani_plus_germline_18kmeans.txt"
                         , header = T, sep = "\t", stringsAsFactors = F)
gilbert.sc.aggr <- aggregate(.~kmeans.cluster, data=gilbert.sc, FUN=mean)

sample.clust1 <- apply(gilbert.sc.aggr[,c("X46C", "D3", "TT2", "iPSC", "iPSC_1D4"
                                          , "iPSC_2D4", "EPL")], 1, mean)
sample.clust2 <-  apply(gilbert.sc.aggr[,c("EBM3", "EpiSC5", "EpiSC7")], 1, mean)
sample.clust3 <- apply(gilbert.sc.aggr[,c("EBM6", "X46CNPC", "TT2NPC", "EBM9" )], 1, mean)
sample.clust4 <- apply(gilbert.sc.aggr[,c("Mesoderm", "Endoderm")], 1, mean)
sample.clust5 <- apply(gilbert.sc.aggr[,c("piPSC_1A2", "piPSC_1B3", "piPSC_V3")], 1, mean)
sample.clust6 <- apply(gilbert.sc.aggr[,c("MEF_female", "MEF_male", "Myoblast")], 1, mean)
sample.clust7 <- gilbert.sc.aggr$PGC.female.1
sample.clust8 <- gilbert.sc.aggr$PGC.male.1
sample.clust9 <- apply(gilbert.sc.aggr[,c("sperm.1", "sperm.2")], 1, mean)

samples.clust.annot <- data.frame(cluster1 = sample.clust1
                                  , cluster2 = sample.clust2
                                  , cluster3 = sample.clust3
                                  , cluster4 = sample.clust4
                                  , cluster5 = sample.clust5
                                  , cluster6 = sample.clust6
                                  , cluster7 = sample.clust7
                                  , cluster8 = sample.clust8
                                  , cluster9 = sample.clust9
                                  , stringsAsFactors = F)

enh_rt <- read.table(file = "mouse_enh_meanRT_by_germ_line_cons4sp.txt", header =T
                     , stringsAsFactors = F, sep = '\t')

# separate tissue column into multiple rows
library(tidyr)
library(dplyr)

enh_rt <- enh_rt %>%
  mutate(tissue = strsplit(as.character(tissue), ",")) %>%
  unnest(tissue)

enh_rt$type <- tolower(enh_rt$type)
enh_rt <- as.data.frame(enh_rt)

#separate into active and poised
active_rt <- enh_rt[enh_rt$type == "active", c("id", "tissue", "age", "PGC.female.1.bedGraph"
                                               , "PGC.male.1.bedGraph", "Sperm.1.bedGraph"
                                               , "Sperm.2.bedGraph")]
poised_rt <- enh_rt[enh_rt$type == "poised", c("id", "tissue","age", "PGC.female.1.bedGraph"
                                               , "PGC.male.1.bedGraph", "Sperm.1.bedGraph"
                                               , "Sperm.2.bedGraph")]

active_rt$group <- paste(active_rt$tissue, active_rt$age, sep = '_')
poised_rt$group <- paste(poised_rt$tissue, poised_rt$age, sep = '_')


active_rt$tissue <- NULL
active_rt$type <- NULL
poised_rt$tissue <- NULL
poised_rt$type <- NULL
active_rt$id <- NULL
poised_rt$id <- NULL
active_rt$age <- NULL
poised_rt$age <- NULL

poised_rt.heat <- aggregate(.~group, poised_rt, mean)
active_rt.heat <- aggregate(.~group, active_rt, mean)

rownames(poised_rt.heat) <- poised_rt.heat$group
rownames(active_rt.heat) <- active_rt.heat$group

poised_rt.heat$group <- NULL
active_rt.heat$group <- NULL

annotation_active <- data.frame(Tissue = (gsub("_.*", "", rownames(active_rt.heat)))
                                , Type = (gsub(".*_", "", rownames(active_rt.heat))))
annotation_poised <- data.frame(Tissue = (gsub("_.*", "", rownames(poised_rt.heat)))
                                , Type = (gsub(".*_", "", rownames(poised_rt.heat))))

rownames(annotation_active) <- rownames(active_rt.heat)
rownames(annotation_poised) <- rownames(poised_rt.heat)

steps <- rev(brewer.pal(n = 8, name = "RdYlBu"))
pal <- color.palette(rev(steps), c(20, 20,5,5,5, 20,20), space="rgb")

all_rt_values <- c(poised_rt.heat$PGC.female.1.bedGraph, poised_rt.heat$PGC.male.1.bedGraph
                   , poised_rt.heat$Sperm.1.bedGraph, poised_rt.heat$Sperm.2.bedGraph
                   , active_rt.heat$PGC.female.1.bedGraph, active_rt.heat$PGC.male.1.bedGraph
                   , active_rt.heat$Sperm.1.bedGraph, active_rt.heat$Sperm.2.bedGraph)


library("Hmisc")
rt_bins <- cut2(all_rt_values, g=50, onlycuts=T)


pdf("mouse_poised_rt_by_tissue_heatmap_matchedColour_cons4sp.pdf")
pheatmap(poised_rt.heat, cluster_rows=TRUE, cluster_cols=TRUE
         , color = pal(50), annotation_row = annotation_poised
         , breaks = rt_bins)
dev.off()

pdf("mouse_active_rt_by_tissue_heatmap_matchedColour_cons4sp.pdf")
pheatmap(active_rt.heat, cluster_rows=TRUE, cluster_cols=TRUE
         , color = pal(50), annotation_row = annotation_active
         , breaks = rt_bins)
dev.off()
