# Fisher's exact tests for the differential overlap of recent and conserved mouse enhancers with accessible regions of the genome (DNAse-seq data)

library(GenomicRanges)
library(tidyr)
library(dplyr)
library(reshape2)
library(data.table)

active_cons <- read.table(file = "activeEnhancer_all_tissues_conserved_with_mouse_withTissueDup"
                          , header = T, stringsAsFactors = F)#59398     6
active_recent <- read.table(file = "activeEnhancer_all_tissues_mouse_only"
                            , header = T, stringsAsFactors = F)#31618     7
poised_cons <- read.table(file = "primedEnhancer_all_tissues_conserved_with_mouse_withTissueDup"
                          , header = T, stringsAsFactors = F)#65703     6
poised_recent <- read.table(file = "primedEnhancer_all_tissues_mouse_only"
                            , header = T, stringsAsFactors = F)#62315     7

active_cons <- unique(active_cons)#59398     6
active_recent <- unique(active_recent)#25610     7
poised_cons <- unique(poised_cons)#65703     6
poised_recent <- unique(poised_recent)#55294     7

active_recent.long <- active_recent %>%
  mutate(V6 = strsplit(as.character(active_recent$V6), ",")) %>%
  unnest(V6)

active_recent.long <- as.data.frame(active_recent.long)
active_recent.long$tissue <- sub(pattern = "_H.*", "", active_recent.long$V6)
active_recent.long$tissue <- sub(pattern = "Mouse_", "", active_recent.long$tissue)
active_recent.long$V6 <- sub(pattern = ".*_H", "H", active_recent.long$V6)
active_recent.long$V6 <- sub(pattern = "_.*", "", active_recent.long$V6)#32476     8

poised_recent.long <- poised_recent %>%
  mutate(V6 = strsplit(as.character(poised_recent$V6), ",")) %>%
  unnest(V6)

poised_recent.long <- as.data.frame(poised_recent.long)
poised_recent.long$tissue <- sub(pattern = "_H.*", "", poised_recent.long$V6)
poised_recent.long$tissue <- sub(pattern = "Mouse_", "", poised_recent.long$tissue)
poised_recent.long$V6 <- sub(pattern = ".*_H", "H", poised_recent.long$V6)
poised_recent.long$V6 <- sub(pattern = "_.*", "", poised_recent.long$V6)#62600     8

# keep only coordinates and tissue
active_cons <- unique(active_cons[,c("chr", "start", "end", "type")])
active_recent.long <- unique(active_recent.long[,c("V1", "V2", "V3", "tissue")])
poised_cons <- unique(poised_cons[,c("chr", "start", "end", "type")])
poised_recent.long <- unique(poised_recent.long[,c("V1", "V2", "V3", "tissue")])

## same column names
colnames(active_cons) <- c("chr", "start", "end", "tissue")
colnames(active_recent.long) <- c("chr", "start", "end", "tissue")
colnames(poised_cons) <- c("chr", "start", "end", "tissue")
colnames(poised_recent.long) <- c("chr", "start", "end", "tissue")

active_recent.long$tissue <- tolower(active_recent.long$tissue)
poised_recent.long$tissue <- tolower(poised_recent.long$tissue)

# add type (conservation)
active_cons$type <- "conserved"
active_recent.long$type <- "recent"
poised_cons$type <- "conserved"
poised_recent.long$type <- "recent"

# separate enhancers by tissue
# no available adult testis atac-seq or dnase-seq 
brain <- rbind(active_cons[active_cons$tissue=="brain", ]
               , active_recent.long[active_recent.long$tissue=="brain",]
               , poised_cons[poised_cons$tissue=="brain",]
               , poised_recent.long[poised_recent.long$tissue=="brain",])#59458     5
muscle <- rbind(active_cons[active_cons$tissue=="muscle", ]
                , active_recent.long[active_recent.long$tissue=="muscle",]
                , poised_cons[poised_cons$tissue=="muscle",]
                , poised_recent.long[poised_recent.long$tissue=="muscle",])#58677     5
liver <- rbind(active_cons[active_cons$tissue=="liver", ]
               , active_recent.long[active_recent.long$tissue=="liver",]
               , poised_cons[poised_cons$tissue=="liver",]
               , poised_recent.long[poised_recent.long$tissue=="liver",])#53943     5

###

### reading DNAse-seq
brain_dnase <- fread("./encode/ENCFF185BXA_brain_adult_mm10.bed.gz")# 258549     10
liver_dnase <- fread("./encode/ENCFF542LZZ_liver_adult_mm10.bed.gz")# 180516     10
muscle_dnase <- fread("./encode/ENCFF357NZL_muscle_adult_mm10.bed.gz")# 83795    10

##  test using a minimum percentage of basepairs in enhancers overlapping the DNAse-seq data
### Note: provide the percentage in a range from 0 to 1
overlap_with_dnasev2 <- function(enh, dnase, perc){
  dnase <- as.data.frame(unique(dnase[,1:3]))
  # make gr objects
  enh_gr <- with(enh, GRanges( chr , IRanges( start, end )))
  dnase_gr <- with(dnase, GRanges( V1 , IRanges( V2, V3 )))
  hits <- (findOverlaps(enh_gr, dnase_gr))
  overlaps <- pintersect(enh_gr[queryHits(hits)], dnase_gr[subjectHits(hits)])
  # Note: percentage of overlap will be calculated with respect to the enhancer length
  percentOverlap <- width(overlaps) / width(enh_gr[queryHits(hits)])
  hits <- hits[percentOverlap >= perc]
  hits <- as.data.frame(hits)
  enh$class <- ""
  enh[unique(hits$queryHits), "class"] <- "DNAse"
  enh[-unique(hits$queryHits), "class"] <- "NOT_DNAse"
  enh_freq <- as.data.frame(table(enh[,c("type", "class")]))
  enh_freq <- tidyr::spread(enh_freq, class, Freq)
  rownames(enh_freq) <- enh_freq$type
  enh_freq$type <- NULL
  print(enh_freq)
  fisher_out <- fisher.test(as.matrix(enh_freq), alternative = "greater")
  return(fisher_out)
}



###### ********************************************######
######   AT LEAST 30% OF ENHANCER base pairs OVERLAP #####
###### ********************************************######

brain_fisher_30perc <- overlap_with_dnasev2(brain, brain_dnase, 0.3)
#           DNAse NOT_DNAse
# conserved  3204     33944
# recent     3144     19166

brain_fisher_30perc

# Fisher's Exact Test for Count Data
#
# data:  as.matrix(enh_freq)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.5505211       Inf
# sample estimates:
# odds ratio
#  0.5754456
#  

muscle_fisher_30perc <- overlap_with_dnasev2(muscle, muscle_dnase, 0.3)
#           DNAse NOT_DNAse
# conserved   745     35315
# recent      861     21756

muscle_fisher_30perc
# Fisher's Exact Test for Count Data
#
# data:  as.matrix(enh_freq)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.4896707       Inf
# sample estimates:
# odds ratio
#  0.5330478
#  

liver_fisher_30perc <- overlap_with_dnasev2(liver, liver_dnase, 0.3)
#     DNAse NOT_DNAse
# conserved  2392     26755
# recent     4838     19958

liver_fisher_30perc
# Fisher's Exact Test for Count Data
#
# data:  as.matrix(enh_freq)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.3528561       Inf
# sample estimates:
# odds ratio
#  0.3688201
