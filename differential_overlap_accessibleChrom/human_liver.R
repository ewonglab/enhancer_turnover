# Fisher's exact tests for the the differential overlap of recent and conserved enhancers with accessible regions of the genome (ATAC-seq)

library(data.table)
library("GenomeInfoDb")
library(GenomicRanges)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

#read enhancers
enhancers <- fread("Hsap_Enhancers_conservationInOthers")
recent_e <- fread("Hsap_H3K27ac_humanspEnhancers")
#read promoters
promoters <- fread("Hsap_Promoters_conservationInOthers")
recent_p <- fread("Hsap_H3K4me3all_humanspPromoters")

enhancers <- as.data.frame(enhancers)
recent_e <- as.data.frame(recent_e)
promoters <- as.data.frame(promoters)
recent_p <- as.data.frame(recent_p)

# subset enhancers and promoters that are conserved with at least 2 species
species <- c("Fcat", "Csab", "Mdom", "Ddel", "Mbid", "Cpor"
             , "Bbor", "Rnor", "Btau", "Cjac", "Cfam", "Mmus"
             , "Sscr", "Ocun", "Tbel", "Shar", "Mmul", "Mfur")

enhancers$n_species_active <- apply(enhancers[,species],1, function(x) length(which(x!="-" & !is.na(x))))
promoters$n_species_active <- apply(promoters[,species],1, function(x) length(which(x!="-" & !is.na(x))))

cons_p <- promoters[promoters$n_species_active >= 2,]# 10166    24
cons_e <- enhancers[enhancers$n_species_active >= 2,]#13329    24

cons_p$Chrom <- paste0("chr", cons_p$Chrom)
cons_e$Chrom <- paste0("chr", cons_e$Chrom)
recent_e$Chrom <- paste0("chr", recent_e$Chrom)
recent_p$Chrom <- paste0("chr", recent_p$Chrom)

cons_p <- unique(cons_p[,c("Chrom", "Start", "End")])
cons_e <- unique(cons_e[,c("Chrom", "Start", "End")])
recent_e <- unique(recent_e[,c("Chrom", "Start", "End")])
recent_p <- unique(recent_p[,c("Chrom", "Start", "End")])


# bed narrow ATAC-seq peaks
hg19_ATAC_peaks <- read.table(file = "ENCFF631JIS_ATAC_hg19.bed"
                              , header = F, stringsAsFactors = F, sep = '\t')# 279525      3

# add enhancer type (conservation)
cons_p$type <- "conserved"
recent_p$type <- "recent"
cons_e$type <- "conserved"
recent_e$type <- "recent"

liver_enh <- rbind(cons_e, recent_e)
liver_prom <- rbind(cons_p, recent_p)


# minimum percentage of overlap
overlap_with_ATACv2 <- function(enh, dnase, perc){
  dnase <- as.data.frame(unique(dnase[,1:3]))
  # make gr objects
  enh_gr <- with(enh, GRanges( Chrom , IRanges( Start +1, End )))
  dnase_gr <- with(dnase, GRanges( V1 , IRanges( V2 +1, V3 )))
  hits <- (findOverlaps(enh_gr, dnase_gr))
  overlaps <- pintersect(enh_gr[queryHits(hits)], dnase_gr[subjectHits(hits)])
  # Note: percentage of overlap will be calculated with respect to the enhancer length
  percentOverlap <- width(overlaps) / width(enh_gr[queryHits(hits)])
  hits <- hits[percentOverlap >= perc]
  hits <- as.data.frame(hits)
  enh$class <- ""
  enh[unique(hits$queryHits), "class"] <- "ATAC"
  enh[-unique(hits$queryHits), "class"] <- "NOT_ATAC"
  enh_freq <- as.data.frame(table(enh[,c("type", "class")]))
  enh_freq <- tidyr::spread(enh_freq, class, Freq)
  rownames(enh_freq) <- enh_freq$type
  enh_freq$type <- NULL
  print(enh_freq)
  fisher_out <- fisher.test(as.matrix(enh_freq), alternative = "greater")
  return(fisher_out)
}


### AT LEAST 30% OVERLAP
prom_fisher_0.3 <- overlap_with_ATACv2(liver_prom, hg19_ATAC_peaks, 0.3)
#           ATAC NOT_ATAC
# conserved 7047     3119
# recent     460      334

prom_fisher_0.3
# Fisher's Exact Test for Count Data
# 
# data:  as.matrix(enh_freq)
# p-value = 4.997e-11
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.44559     Inf
# sample estimates:
# odds ratio 
#    1.64044 

enh_fisher_0.3 <- overlap_with_ATACv2(liver_enh, hg19_ATAC_peaks, 0.3)
#           ATAC NOT_ATAC
# conserved 5119     8210
# recent    4451     5983

enh_fisher_0.3
# Fisher's Exact Test for Count Data
# 
# data:  as.matrix(enh_freq)
# p-value = 1
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.8019218       Inf
# sample estimates:
# odds ratio 
#  0.8381256 

