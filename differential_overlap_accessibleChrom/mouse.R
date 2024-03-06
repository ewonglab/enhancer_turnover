# Fisher's exact tests for the differential overlap of recent and conserved mouse enhancers with accessible regions of the genome (DNAse-seq data)

library(GenomicRanges)
library(tidyr)
library(dplyr)
library(reshape2)
library(data.table)

mouse.enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                        , header = F, stringsAsFactors = F
                        , sep = "\t")

recent_enh <- mouse.enh[mouse.enh$V6 == "Recent",]
active_atLeast2sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast2sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")

cons_enh <- rbind(active_atLeast2sp, poised_atLeast2sp)

library(tidyr)
library(dplyr)

recent_enh <- recent_enh %>%
  mutate(V7 = strsplit(as.character(V7), ",")) %>%
  unnest(V7)

cons_enh <- cons_enh %>%
  mutate(tissue = strsplit(as.character(tissue), ",")) %>%
  unnest(tissue)

recent_enh <- as.data.frame(recent_enh)
cons_enh <- as.data.frame(cons_enh)

cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end) # unique rows

# keep only coordinates and tissue
cons_enh <- unique(cons_enh[,c("chr", "start", "end", "tissue")]) 
recent_enh <- unique(recent_enh[,c("V1", "V2", "V3", "V7")])

## same column names
colnames(cons_enh) <- c("chr", "start", "end", "tissue")
colnames(recent_enh) <- c("chr", "start", "end", "tissue")

cons_enh$tissue <- tolower(cons_enh$tissue)
recent_enh$tissue <- tolower(recent_enh$tissue)

# add type (conservation)
cons_enh$type <- "conserved"
recent_enh$type <- "recent"

brain <- rbind(cons_enh[cons_enh$tissue=="brain", ]
               , recent_enh[recent_enh$tissue=="brain",])
muscle <- rbind(cons_enh[cons_enh$tissue=="muscle", ]
                , recent_enh[recent_enh$tissue=="muscle",])
liver <- rbind(cons_enh[cons_enh$tissue=="liver", ]
               , recent_enh[recent_enh$tissue=="liver",])

### reading DNAse-seq 
brain_dnase <- fread("./encode/ENCFF185BXA_brain_adult_mm10.bed.gz")
liver_dnase <- fread("./encode/ENCFF542LZZ_DNAse_liver_mm10.bed.gz")
muscle_dnase <- fread("./encode/ENCFF357NZL_muscle_adult_mm10.bed.gz") 

overlap_with_dnasev2 <- function(enh, dnase, perc){
  dnase <- as.data.frame(unique(dnase[,1:3]))
  # make gr objects
  enh_gr <- with(enh, GRanges( chr , IRanges( start + 1, end )))
  dnase_gr <- with(dnase, GRanges( V1 , IRanges( V2 + 1, V3 )))
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


brain_fisher_30perc <- overlap_with_dnasev2(brain, brain_dnase, 0.3)
#           DNAse NOT_DNAse
# conserved  3002     31782
# recent     3125     19185

brain_fisher_30perc
# Fisher's Exact Test for Count Data
# 
# data:  as.matrix(enh_freq)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.5543441       Inf
# sample estimates:
# odds ratio 
#  0.5799145 

brain_fisher_30perc$p.value
# [1] 1



muscle_fisher_30perc <- overlap_with_dnasev2(muscle, muscle_dnase, 0.3)
#           DNAse NOT_DNAse
# conserved   671     32822
# recent      856     21761

muscle_fisher_30perc

# Fisher's Exact Test for Count Data
# 
# data:  as.matrix(enh_freq)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.4762088       Inf
# sample estimates:
# odds ratio 
#  0.5197072 
# 
muscle_fisher_30perc$p.value
# [1] 1


liver_fisher_30perc <- overlap_with_dnasev2(liver, liver_dnase, 0.3)
#           DNAse NOT_DNAse
# conserved  2187     24851
# recent     4806     19990

liver_fisher_30perc

# Fisher's Exact Test for Count Data
# 
# data:  as.matrix(enh_freq)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.3497271       Inf
# sample estimates:
# odds ratio 
#  0.3660509 

liver_fisher_30perc$p.value
# [1] 1

