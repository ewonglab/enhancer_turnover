library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')
library("Formula", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library("Hmisc", lib = "/scratch/zk16/cc3704/R_packages/R_4.0.0")
library(data.table)
library(reshape2)


#Loading macs3 output
chrOrder<-paste0("chr", c(1:22,"X","Y"))
#xsl files from macs are 1-based
dir = "/g/data/zk16/veronika/projects/repli_time/AML/ATAC/peaks/"
sample <- c("SU484_blast", "SU484_pHSC", "SU501_blast"
            , "SU501_pHSC", "SU654_blast", "SU654_pHSC")
macsPeaks <- list()
macsPeaks.gr <- GRangesList()
for (i in 1:length(sample)){
  print(sample[i])
  macsPeaks[[i]] <- read.table(paste0(dir, sample[i], "_peaks.xls"), header = T)
  # print(dim(macsPeaks[[i]]))
  macsPeaks[[i]]$chr <- paste0("chr", macsPeaks[[i]]$chr)
  macsPeaks[[i]] <- macsPeaks[[i]][macsPeaks[[i]]$chr %in% chrOrder,]
  # print(dim(macsPeaks[[i]]))
  macsPeaks.gr[[i]] <- makeGRangesFromDataFrame(macsPeaks[i], keep.extra.columns = T)
}

setwd("/g/data/zk16/cc3704/replication_timing/human/AML")
SU484_blast.gr <- macsPeaks.gr[[1]]
SU484_pHSC.gr <- macsPeaks.gr[[2]]
SU501_blast.gr <- macsPeaks.gr[[3]]
SU501_pHSC.gr <- macsPeaks.gr[[4]]
SU654_blast.gr <- macsPeaks.gr[[5]]
SU654_pHSC.gr <- macsPeaks.gr[[6]]

#keep peaks only in normal and common between normal and cancer
SU484_normal <- subsetByOverlaps(SU484_pHSC.gr, SU484_blast.gr, invert = TRUE)
SU501_normal <- subsetByOverlaps(SU501_pHSC.gr, SU501_blast.gr, invert = TRUE)
SU654_normal <- subsetByOverlaps(SU654_pHSC.gr, SU654_blast.gr, invert = TRUE)

SU484_olaps <- findOverlaps(SU484_blast.gr, SU484_pHSC.gr)
SU484_common_1 <- SU484_blast.gr[queryHits(SU484_olaps)]
SU484_common_2 <- SU484_pHSC.gr[subjectHits(SU484_olaps)]
SU484_common <- GenomicRanges::union(SU484_common_1, SU484_common_2)

SU501_olaps <- findOverlaps(SU501_blast.gr, SU501_pHSC.gr)
SU501_common_1 <- SU501_blast.gr[queryHits(SU501_olaps)]
SU501_common_2 <- SU501_pHSC.gr[subjectHits(SU501_olaps)]
SU501_common <- GenomicRanges::union(SU501_common_1, SU501_common_2)

SU654_olaps <- findOverlaps(SU654_blast.gr, SU654_pHSC.gr)
SU654_common_1 <- SU654_blast.gr[queryHits(SU654_olaps)]
SU654_common_2 <- SU654_pHSC.gr[subjectHits(SU654_olaps)]
SU654_common <- GenomicRanges::union(SU654_common_1, SU654_common_2)

common_peaks <- GenomicRanges::union(SU484_common, c(SU501_common, SU654_common))
normal_peaks <- GenomicRanges::union(SU484_normal, c(SU501_normal, SU654_normal))

#IMPORTANT NOTE: THE UNION METHOD REDUCES GENOMIC RANGES SO MY MAPPING BACK TO MACS
# OUTPUT WON'T BE ONE-TO-ONE

#OVERLAP SELECTED RANGES WITH MACS OUTPUT ANS CALCULATE AVERAGE SIGNAL FOR EVER REGION
names(macsPeaks) <- sample
macsPeaks.df <- do.call("rbind", macsPeaks)
dim(macsPeaks.df)#466086     10
dim(unique(macsPeaks.df[,c("chr", "start", "end")]))#465372      3
#macs output is already 1-based
macsPeaks.df.gr <- with(macsPeaks.df, GRanges(chr , IRanges( start, end )))
macsPeaks.df$sample <- rownames(macsPeaks.df)

x <- as.data.frame(findOverlaps(common_peaks, macsPeaks.df.gr))
common_signal <- cbind(as.data.frame(common_peaks)[x$queryHits,]
                       , macsPeaks.df[x$subjectHits, c("pileup", "fold_enrichment", "sample")])
#start position back to 0-based
common_signal$start <- common_signal$start-1
common_signal$sample <- gsub("\\..*", "", common_signal$sample)
unique(common_signal$sample)
# [1] "SU484_blast" "SU484_pHSC"  "SU501_blast" "SU501_pHSC"  "SU654_blast"
# [6] "SU654_pHSC"

#keeping only overlaps with normal samples
common_signal <- common_signal[common_signal$sample %in%
                                 c("SU484_pHSC", "SU501_pHSC", "SU654_pHSC" ),]
#remove sample and calculate mean signal
common_signal$sample <- NULL
common_mean_pileup <- aggregate(pileup ~ .
                                , common_signal[,c("seqnames", "start", "end", "width"
                                                   , "strand", "pileup")], mean)#66575     6
#fold enrichment
common_mean_FE <- aggregate(fold_enrichment ~. ,common_signal[c("seqnames", "start", "end", "width"
                                             , "strand", "fold_enrichment")], mean)#66575     6

x <- as.data.frame(findOverlaps(normal_peaks, macsPeaks.df.gr))
normal_signal <- cbind(as.data.frame(normal_peaks)[x$queryHits,]
                       , macsPeaks.df[x$subjectHits, c("pileup", "fold_enrichment", "sample")])
#start position back to 0-based
normal_signal$start <- normal_signal$start-1
normal_signal$sample <- gsub("\\..*", "", normal_signal$sample)
unique(normal_signal$sample)
normal_signal <- normal_signal[normal_signal$sample %in%
                                 c("SU484_pHSC", "SU501_pHSC", "SU654_pHSC" ),]
normal_signal$sample <- NULL
normal_mean_pileup <- aggregate(pileup ~ .
                                , normal_signal[,c("seqnames", "start", "end", "width"
                                                   , "strand", "pileup")], mean)#57491     6
normal_mean_FE <- aggregate(fold_enrichment ~., normal_signal[c("seqnames", "start"
                                                                , "end", "width"
                                                                , "strand", "fold_enrichment")]
                            , mean)#57491     6

#correlation between pileup and loggfold enrichment
tmp <- merge(common_mean_pileup, common_mean_FE, by = c("seqnames", "start", "end", "width", "strand"))
cor.test(tmp$pileup, tmp$fold_enrichment)
# Pearson's product-moment correlation
#
# data:  tmp$pileup and tmp$fold_enrichment
# t = 206.73, df = 66573, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6206363 0.6298889
# sample estimates:
#       cor
# 0.6252846

tmp <- merge(normal_mean_pileup, normal_mean_FE, by = c("seqnames", "start", "end", "width", "strand"))
# Pearson's product-moment correlation
#
# data:  tmp$pileup and tmp$fold_enrichment
# t = 212.15, df = 57489, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6580534 0.6672232
# sample estimates:
#       cor
# 0.6626631

#ABout the fold enrichment in MACS
#fold enrichment for this peak summit against random Poisson distribution with local lambda, -log10(qvalue) at peak summit
#it is basically the signal compared to local background
#https://www.biostars.org/p/318105/

#bind signal of common and normal peaks
#there are some duplicated peaks :P
peaks <- cbind(common_mean_FE, normal_mean_FE)
peaks <- aggregate(fold_enrichment ~., peaks, mean)#124044      6
#I don't need to reduce peaks

#scale signal values
peaks$fold_enrichment <- scale(x = peaks$fold_enrichment, center = T, scale = T)
#save fold enrichment values and run replicon
write.table(x = peaks[,c("seqnames", "start", "end", "fold_enrichment")]
            , file = "AML_pHSC_mean_FE.bed", quote = F, col.names = F, row.names = F
            , sep = '\t')
cd /g/data/zk16/cc3704/replication_timing/human
java -jar /g/data/zk16/cc3704/tools/RepliconWrench.jar BedToIPLS -i ./AML/AML_pHSC_mean_FE.bed -o AML_hg19_IPLS
# 19:58:14,349  INFO BedToIPLS:137 - Processing IPLS for chromosome chr9
# 19:58:14,747  INFO BedToIPLS:137 - Processing IPLS for chromosome chrY
# 19:58:14,856  INFO BedToIPLS:137 - Processing IPLS for chromosome chr7
# 19:58:15,043  INFO BedToIPLS:137 - Processing IPLS for chromosome chr8
# 19:58:15,267  INFO BedToIPLS:137 - Processing IPLS for chromosome chrX
# 19:58:15,395  INFO BedToIPLS:137 - Processing IPLS for chromosome chr5
# 19:58:15,633  INFO BedToIPLS:137 - Processing IPLS for chromosome chr6
# 19:58:15,842  INFO BedToIPLS:137 - Processing IPLS for chromosome chr10
# 19:58:15,984  INFO BedToIPLS:137 - Processing IPLS for chromosome chr3
# 19:58:16,150  INFO BedToIPLS:137 - Processing IPLS for chromosome chr21
# 19:58:16,225  INFO BedToIPLS:137 - Processing IPLS for chromosome chr11
# 19:58:16,350  INFO BedToIPLS:137 - Processing IPLS for chromosome chr22
# 19:58:16,425  INFO BedToIPLS:137 - Processing IPLS for chromosome chr4
# 19:58:16,603  INFO BedToIPLS:137 - Processing IPLS for chromosome chr12
# 19:58:16,759  INFO BedToIPLS:137 - Processing IPLS for chromosome chr1
# 19:58:17,004  INFO BedToIPLS:137 - Processing IPLS for chromosome chr2
# 19:58:17,219  INFO BedToIPLS:137 - Processing IPLS for chromosome chr13
# 19:58:17,332  INFO BedToIPLS:137 - Processing IPLS for chromosome chr20
# 19:58:17,414  INFO BedToIPLS:137 - Processing IPLS for chromosome chr18
# 19:58:17,500  INFO BedToIPLS:137 - Processing IPLS for chromosome chr19
# 19:58:17,588  INFO BedToIPLS:137 - Processing IPLS for chromosome chr14
# 19:58:17,694  INFO BedToIPLS:137 - Processing IPLS for chromosome chr15
# 19:58:17,800  INFO BedToIPLS:137 - Processing IPLS for chromosome chr16
# 19:58:17,891  INFO BedToIPLS:137 - Processing IPLS for chromosome chr17


#Note: liver enhancers from villar et al are located in main chromosomes
module load gcc
module load gsl
module load R/4.0.0
cd /g/data/zk16/cc3704/replication_timing
export PATH=/g/data/zk16/xzhang/bin:$PATH

aml_path="/g/data/zk16/cc3704/replication_timing/human/AML_hg19_IPLS/"
aml_base="AML_pHSC_mean_FE"

cd ${aml_path}
# replicon -threads 8 ./chr1/chr1_${aml_base} chr1
# replicon -threads 8 ./chr15/chr15_${aml_base} chr15
# replicon -threads 8 ./chr4/chr4_${aml_base} chr4
# replicon -threads 8 ./chrX/chrX_${aml_base} chrX
# replicon -threads 8 ./chr10/chr10_${aml_base} chr10
# replicon -threads 8 ./chr16/chr16_${aml_base} chr16
# replicon -threads 8 ./chr2/chr2_${aml_base} chr2
# replicon -threads 8 ./chr9/chr9_${aml_base} chr9
# replicon -threads 8 ./chrY/chrY_${aml_base} chrY
# replicon -threads 8 ./chr11/chr11_${aml_base} chr11
# replicon -threads 8 ./chr17/chr17_${aml_base} chr17
# replicon -threads 8 ./chr20/chr20_${aml_base} chr20
# replicon -threads 8 ./chr5/chr5_${aml_base} chr5
# replicon -threads 8 ./chr12/chr12_${aml_base} chr12
# replicon -threads 8 ./chr21/chr21_${aml_base} chr21
# replicon -threads 8 ./chr6/chr6_${aml_base} chr6
# replicon -threads 8 ./chr13/chr13_${aml_base} chr13
# replicon -threads 8 ./chr18/chr18_${aml_base} chr18
# replicon -threads 8 ./chr22/chr22_${aml_base} chr22
# replicon -threads 8 ./chr7/chr7_${aml_base} chr7
# replicon -threads 8 ./chr14/chr14_${aml_base} chr14
# replicon -threads 8 ./chr19/chr19_${aml_base} chr19
# replicon -threads 8 ./chr3/chr3_${aml_base} chr3
# replicon -threads 8 ./chr8/chr8_${aml_base} chr8

#CGHnTiming to timing
rep_path="/g/data/zk16/cc3704/tools/"
cd /g/data/zk16/cc3704/replication_timing/human/AML_hg19_IPLS
for cgh in *.CGHnTiming.csv
do
new_name="${cgh/CGHnTiming.csv/}"
java -jar ${rep_path}RepliconWrench.jar CghToTiming -i ${cgh} > ${new_name}timing
done

