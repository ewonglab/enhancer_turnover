library("rlang", lib= "/srv/scratch/z5291980/extra_R_packages")
library("Rcpp", lib="/srv/scratch/z5291980/extra_R_packages")
library("lifecycle", lib="/srv/scratch/z5291980/extra_R_packages")
library("pillar", lib="/srv/scratch/z5291980/extra_R_packages")
library("ellipsis", lib="/srv/scratch/z5291980/extra_R_packages")
library(data.table)
library("dplyr", lib="/srv/scratch/z5291980/extra_R_packages")
library(ggplot2)
library(reshape2)

library("matrixStats", lib = "/srv/scratch/z5291980/extra_R_packages")
library("MatrixGenerics", lib= "/srv/scratch/z5291980/extra_R_packages")
library(BiocGenerics, lib = "/srv/scratch/z5291980/extra_R_packages")
library(S4Vectors, lib = "/srv/scratch/z5291980/extra_R_packages")
library(IRanges, lib = "/srv/scratch/z5291980/extra_R_packages")
library(GenomeInfoDb, lib = "/srv/scratch/z5291980/extra_R_packages")
library(GenomicRanges, lib = "/srv/scratch/z5291980/extra_R_packages")

library(XVector, lib = "/srv/scratch/z5291980/extra_R_packages")
library(Biostrings, lib = "/srv/scratch/z5291980/extra_R_packages")
library("BiocManager", lib="/srv/scratch/z5291980/extra_R_packages")
library("rtracklayer", lib = "/srv/scratch/z5291980/extra_R_packages")
library(BSgenome, lib = "/srv/scratch/z5291980/extra_R_packages")
library(BSgenome.Hsapiens.UCSC.hg19, lib = "/srv/scratch/z5291980/extra_R_packages")
library("Hmisc")

library("biomaRt", lib="/srv/scratch/z5291980/extra_R_packages")
library("Biobase", lib= "/srv/scratch/z5291980/extra_R_packages")
library("SummarizedExperiment", lib= "/srv/scratch/z5291980/extra_R_packages")
library("Rsamtools", lib= "/srv/scratch/z5291980/extra_R_packages")
library("VariantAnnotation", lib= "/srv/scratch/z5291980/extra_R_packages")
library("AnnotationDbi", lib= "/srv/scratch/z5291980/extra_R_packages")
library("GenomicFeatures", lib="/srv/scratch/z5291980/extra_R_packages")
library("ggplot2", lib= "/srv/scratch/z5291980/extra_R_packages")
library("NMF", lib= "/srv/scratch/z5291980/extra_R_packages")
library("SomaticSignatures", lib= "/srv/scratch/z5291980/extra_R_packages")

hg19 <- BSgenome.Hsapiens.UCSC.hg19
# setwd("/srv/scratch/wonglab/paola/replication_timing")
setwd("/srv/scratch/z5291980/replication_timing/")
#only correct positions
vcf <- fread('/srv/scratch/wonglab/paola/replication_timing/novel_human_mutations_cor.bed.gz', header=T)#9488044       8

### MODIFICATION: I WILL SEPARATE THE MUTATIONS INTO THE MUTATIONS AT CpG and NON-CpG sites
chrs <- names(Hsapiens)[1:24]
#getting start positions of CGs
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))
# There were 23 warnings (use warnings() to see them)
vcf_gr <- with(vcf, GRanges( chr , IRanges( end, end )))

x <- as.data.frame(findOverlaps(vcf_gr, cpgr))#1227052       2
length(unique(x$queryHits))#unique
# [1] 1227052
length(unique(x$subjectHits))#non unique
# [1] 1223755

vcf_cpg <- vcf[unique(x$queryHits),]#1227052       8
table(vcf_cpg$node1) #now it is correct! I only have Cs and Gs at human CpGs
# C      G
# 610898 616154

vcf_non.cpg <- vcf[-unique(x$queryHits),]#8260992       8
table(vcf_non.cpg$node1)
# A       C       G       T
# 2427707 1702662 1704307 2426316

vcf_cpg_gr <- with(vcf_cpg, GRanges( chr , IRanges( end, end )))
vcf_non.cpg_gr <- with(vcf_non.cpg, GRanges( chr , IRanges( end, end )))



#############################################################################################

# read in bed files e.g. exons/intergenic/enhancers/early RT/late RT
human_txdb <- makeTxDbFromGFF(file="/srv/scratch/z5291980/human_data/hg19.ensGene.gtf"
                              , format="gtf",organism="Homo sapiens")
hg19_exons <- unique(exons(human_txdb))# 584914
strand(hg19_exons) <- "*"
# I will reduce the enhancers because some of them overlap between them
hg19_exons <- GenomicRanges::reduce(hg19_exons)#315668

exbygene <- exonsBy(human_txdb, "gene")#60234
intergenicRegions <- gaps(unlist(range(exbygene)))# 50142
intergenicRegions <- unique(intergenicRegions)#they are unique but have strand information
strand(intergenicRegions) <- "*"
intergenicRegions <- GenomicRanges::reduce(intergenicRegions)#14696

#do intergenic regions overlap exonic regions ??
x <- as.data.frame(findOverlaps(hg19_exons, intergenicRegions))#they do
#I guess I should remove exons from intergenic regions as they are not exactly intergenic
intergenicRegions <- GenomicRanges::setdiff(intergenicRegions, hg19_exons)#287926
x <- as.data.frame(findOverlaps(hg19_exons, intergenicRegions))#Now they don't overlap

#enhancers
recent <- read.table(file = "/srv/scratch/z5291980/replication_timing/data/human/villar/Hsap_H3K27ac_humanspEnhancers"
                     , header =T, stringsAsFactors = F, sep = '\t')#10434    22
non.recent <- read.table(file = "/srv/scratch/z5291980/replication_timing/data/human/villar/Hsap_H3K4me3.H3K27Ac_overlap_H3K27Aconly_nonRecent"
                         , header =F, stringsAsFactors = F, sep = '\t')#18743     4
recent <- recent[,1:3]
non.recent <- non.recent[,1:3]
colnames(non.recent) <- colnames(recent)
enhancers <- unique(rbind(recent, non.recent))#29177     3
enhancers$Chrom <- paste0("chr", enhancers$Chrom)#all unique enhancers :)
enhancers_gr <- with(enhancers, GRanges( Chrom , IRanges( Start+1, End )))
#reduce enhancers too :(
enhancers_gr <- GenomicRanges::reduce(enhancers_gr)#28295

hg19_exons.df <- as.data.frame(hg19_exons)#315668      5
intergenicRegions.df <- as.data.frame(intergenicRegions)#287926      5
enhancers.df  <- as.data.frame(enhancers_gr)#28295     5

#add ID to every region
hg19_exons.df$ID <- paste0('ID_',seq(1:nrow(hg19_exons.df)))
intergenicRegions.df$ID <- paste0('ID_',seq(1:nrow(intergenicRegions.df)))
enhancers.df$ID <- paste0('ID_',seq(1:nrow(enhancers.df)))

# setwd("/srv/scratch/wonglab/paola/replication_timing")
setwd("/srv/scratch/z5291980/replication_timing/")

#OVERLAP VARISNTS WITH EXONS
x <- as.data.frame(findOverlaps( vcf_gr, hg19_exons))
df <- cbind(vcf[x$queryHits,], hg19_exons.df[x$subjectHits,])#291893     13
df_gr <- with(df, GRanges( chr , IRanges( end, end )))
sca_vr <- VRanges(
  seqnames = as.factor(df$chr),
  ranges =  IRanges(start=df$end, width=1),
  ref = df$node18,
  alt = df$node1,
  study = df$ID)#291893
sca_vr
# sort(table(sca_vr$study), decreasing = TRUE)
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19)
head(sca_motifs)
# can normalized freq (maybe these are useful for SLIM?)
# but if you don't you can get number of mutations if you do colSum

#### SEPARATE sca_motifs INTO THE MUTATIONS AT CpG AND NON-CpG
sca_motifs.cpg <- sca_motifs[ sca_motifs %over% vcf_cpg_gr, ] #43973
sca_motifs.non.cpg <- sca_motifs[ sca_motifs %over% vcf_non.cpg_gr, ] #247920

sca_mm_exons.cpg <- motifMatrix(sca_motifs.cpg, group = "study", normalize = FALSE)#96 30382
sca_mm_exons.non.cpg <- motifMatrix(sca_motifs.non.cpg, group = "study", normalize = FALSE)# 96 91704

write.table(x = sca_mm_exons.cpg, file = "hg19_exons_mut_context_ancestral_CpG.txt"
            , quote = F, sep = '\t')
write.table(x = sca_mm_exons.non.cpg, file = "hg19_exons_mut_context_ancestral_nonCpG.txt"
            , quote = F, sep = '\t')


###### OVERLAP VARIANTS WITH INTERGENIC
x <- as.data.frame(findOverlaps( vcf_gr, intergenicRegions))
df <- cbind(vcf[x$queryHits,], intergenicRegions.df[x$subjectHits,])# 8968847      13
df_gr <- with(df, GRanges( chr , IRanges( end, end )))
sca_vr <- VRanges(
  seqnames = as.factor(df$chr),
  ranges =  IRanges(start=df$end, width=1),
  ref = df$node18,
  alt = df$node1,
  study = df$ID)#8968847
sca_vr
# sort(table(sca_vr$study), decreasing = TRUE)
sca_motifs <- mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19)#8968847

sca_motifs.cpg <- sca_motifs[ sca_motifs %over% vcf_cpg_gr, ] #1153176
sca_motifs.non.cpg <- sca_motifs[ sca_motifs %over% vcf_non.cpg_gr, ]#7815671

sca_mm_inter.cpg <- motifMatrix(sca_motifs.cpg, group = "study", normalize = FALSE)#96 134675
sca_mm_inter.non.cpg <- motifMatrix(sca_motifs.non.cpg, group = "study", normalize = FALSE)#96 200203

write.table(x = sca_mm_inter.cpg, file = "hg19_intergenic_mut_context_ancestral_CpG.txt"
            , quote = F, sep = '\t')
write.table(x = sca_mm_inter.non.cpg, file = "hg19_intergenic_mut_context_ancestral_nonCpG.txt"
            , quote = F, sep = '\t')

###### OVERLAP VARIANTS WITH ENHANCERS
x <- as.data.frame(findOverlaps( vcf_gr, enhancers_gr))
df <- cbind(vcf[x$queryHits,], enhancers.df[x$subjectHits,])#
df_gr <- with(df, GRanges( chr , IRanges( end, end )))
sca_vr <- VRanges(
  seqnames = as.factor(df$chr),
  ranges =  IRanges(start=df$end, width=1),
  ref = df$node18,
  alt = df$node1,
  study = df$ID)#329070
sca_vr
# sort(table(sca_vr$study), decreasing = TRUE)
sca_motifs <- mutationContext(sca_vr, BSgenome.Hsapiens.UCSC.hg19)#329070

sca_motifs.cpg <- sca_motifs[ sca_motifs %over% vcf_cpg_gr, ] #47300
sca_motifs.non.cpg <- sca_motifs[ sca_motifs %over% vcf_non.cpg_gr, ]#281770

sca_mm_enh.cpg <- motifMatrix(sca_motifs.cpg, group = "study", normalize = FALSE)#96 18737
sca_mm_enh.non.cpg <- motifMatrix(sca_motifs.non.cpg, group = "study", normalize = FALSE)#96 26180

write.table(x = sca_mm_enh.cpg, file = "hg19_enhancers_mut_context_ancestral_CpG.txt"
            , quote = F, sep = '\t')
write.table(x = sca_mm_enh.non.cpg, file = "hg19_enhancers_mut_context_ancestral_nonCpG.txt"
            , quote = F, sep = '\t')

# sca_mm_enh.cpg
# sca_mm_enh.non.cpg
# sca_mm_exons.cpg
# sca_mm_exons.non.cpg
# sca_mm_inter.cpg
# sca_mm_inter.non.cpg

library(tidyr)
library(data.table)
#AGGREGATE BY MUTATION TYPE
agg.mutation_type <- function(x){
  x <- as.data.frame(x)
  # colnames(x)
  x$mutation_type <- rownames(x)
  x$mutation_type <- gsub(" .*", "", x$mutation_type)
  print(unique(x$mutation_type))
  # [1] "CA" "CG" "CT" "TA" "TC" "TG"
  #summing mutations by mutation time (ignoring mutation context)
  # print(head(x))
  x.agg <- aggregate(.~mutation_type, x, sum)#6 40925
  rownames(x.agg) <- x.agg$mutation_type
  x.agg$mutation_type <- NULL
  # print(head(x.agg))
  x.agg <- t(x.agg)#245544      3
  # print(head(x.agg))
  # x.agg <- tidyr::spread(data = x.agg, key=mutation_type, value=value, fill=0)#40924     7
  return(as.data.frame(x.agg))}

interg.mut.agg_cpg <- agg.mutation_type(sca_mm_inter.cpg)#
interg.mut.agg_non.cpg <- agg.mutation_type(sca_mm_inter.non.cpg)#
exon.mut.agg_cpg <- agg.mutation_type(sca_mm_exons.cpg)#30382     6
exon.mut.agg_non.cpg <- agg.mutation_type(sca_mm_exons.non.cpg)#
enh.mut.agg_cpg <- agg.mutation_type(sca_mm_enh.cpg)# 18737     6
enh.mut.agg_non.cpg <- agg.mutation_type(sca_mm_enh.non.cpg)# 26180     6

write.table(x = interg.mut.agg_cpg, file = "intergenic_human_anc_mut_CpG", sep = '\t', quote=F)
write.table(x = interg.mut.agg_non.cpg, file = "intergenic_human_anc_mut_non.CpG", sep = '\t', quote=F)
write.table(x = exon.mut.agg_cpg, file = "exon_human_anc_mut_CpG", sep = '\t', quote=F)
write.table(x = exon.mut.agg_non.cpg, file = "exon_human_anc_mut_non.CpG", sep = '\t', quote=F)
write.table(x = enh.mut.agg_cpg, file = "liverEnh_human_anc_mut_CpG", sep = '\t', quote=F)
write.table(x = enh.mut.agg_non.cpg, file = "liverEnh_human_anc_mut_non.CpG", sep = '\t', quote=F)

#######################
#     continue        #
#######################
interg.mut.agg.cpg <- read.table(file = "intergenic_human_anc_mut_CpG", sep = '\t'
                                 , header =T, stringsAsFactors = F)#134603      6
interg.mut.agg.non.cpg <- read.table(file = "intergenic_human_anc_mut_non.CpG", sep = '\t'
                                     , header =T, stringsAsFactors = F)#200107      6

exon.mut.agg.cpg <- read.table(file = "exon_human_anc_mut_CpG", sep = '\t'
                               , header =T, stringsAsFactors = F)#30350     6
exon.mut.agg.non.cpg <- read.table(file = "exon_human_anc_mut_non.CpG", sep = '\t'
                                   , header =T, stringsAsFactors = F)#91618     6

enh.mut.agg.cpg <- read.table(file = "liverEnh_human_anc_mut_CpG", sep = '\t'
                              , header =T, stringsAsFactors = F)#18726     6
enh.mut.agg.non.cpg <- read.table(file = "liverEnh_human_anc_mut_non.CpG", sep = '\t'
                                  , header =T, stringsAsFactors = F)#26171     6

## PLOT NUMBER OF VARIANTS ACROSS REPLICATION TIME WITHOUT NORMALIZING
rt <- read.delim('/srv/scratch/wonglab/data/repliseq/RT_H9_ESC_Ext29405702_hg19.bedgraph'
                 , skip=11, header=F)#2161380       4
names(rt) <- c('chr','start','end','rt')
rt$rtbin <- as.factor(cut2(rt$rt, g=5))
idkey <- data.frame( bin=names(table(rt$rtbin)), quintile=c(1,2,3,4,5))
# bin quintile
# 1 [-2.420,-0.886)        1
# 2 [-0.886,-0.507)        2
# 3 [-0.507, 0.287)        3
# 4 [ 0.287, 1.073)        4
# 5 [ 1.073, 2.151]        5
rt <- merge(rt,  idkey ,by.x='rtbin', by.y='bin', all.x=TRUE)#2161380       6

### MERGE REGIONS LOCATIONS AND NUMBER OF MUTATIONS

#updated to take reduced elements
# read in bed files e.g. exons/intergenic/enhancers/early RT/late RT
human_txdb <- makeTxDbFromGFF(file="/srv/scratch/z5291980/human_data/hg19.ensGene.gtf"
                              , format="gtf",organism="Homo sapiens")
hg19_exons <- unique(exons(human_txdb))# 584914
strand(hg19_exons) <- "*"
# I will reduce the enhancers because some of them overlap between them
hg19_exons <- GenomicRanges::reduce(hg19_exons)#315668

exbygene <- exonsBy(human_txdb, "gene")#60234
intergenicRegions <- gaps(unlist(range(exbygene)))# 50142
intergenicRegions <- unique(intergenicRegions)#they are unique but have strand information
strand(intergenicRegions) <- "*"
intergenicRegions <- GenomicRanges::reduce(intergenicRegions)#14696

#do intergenic regions overlap exonic regions ??
x <- as.data.frame(findOverlaps(hg19_exons, intergenicRegions))#they do
#I guess I should remove exons from intergenic regions as they are not exactly intergenic
intergenicRegions <- GenomicRanges::setdiff(intergenicRegions, hg19_exons)#287926
x <- as.data.frame(findOverlaps(hg19_exons, intergenicRegions))#Now they don't overlap

#enhancers
recent <- read.table(file = "/srv/scratch/z5291980/replication_timing/data/human/villar/Hsap_H3K27ac_humanspEnhancers"
                     , header =T, stringsAsFactors = F, sep = '\t')#10434    22
non.recent <- read.table(file = "/srv/scratch/z5291980/replication_timing/data/human/villar/Hsap_H3K4me3.H3K27Ac_overlap_H3K27Aconly_nonRecent"
                         , header =F, stringsAsFactors = F, sep = '\t')#18743     4
recent <- recent[,1:3]
non.recent <- non.recent[,1:3]
colnames(non.recent) <- colnames(recent)
enhancers <- unique(rbind(recent, non.recent))#29177     3
enhancers$Chrom <- paste0("chr", enhancers$Chrom)#all unique enhancers :)
enhancers_gr <- with(enhancers, GRanges( Chrom , IRanges( Start+1, End )))
#reduce enhancers too :(
enhancers_gr <- GenomicRanges::reduce(enhancers_gr)#28295

hg19_exons.df <- as.data.frame(hg19_exons)#315668      5
intergenicRegions.df <- as.data.frame(intergenicRegions)#287926      5
enhancers.df  <- as.data.frame(enhancers_gr)#28295     5

#add ID to every region
hg19_exons.df$ID <- paste0('ID_',seq(1:nrow(hg19_exons.df)))
intergenicRegions.df$ID <- paste0('ID_',seq(1:nrow(intergenicRegions.df)))
enhancers.df$ID <- paste0('ID_',seq(1:nrow(enhancers.df)))

#HERE: I need to keep all sequences, not only the ones with mutations
interg.mut.agg.cpg <- merge(interg.mut.agg.cpg, intergenicRegions.df
                            , by.x =0, by.y = "ID", all =T)#287926     12
interg.mut.agg.non.cpg <- merge(interg.mut.agg.non.cpg, intergenicRegions.df
                                , by.x =0, by.y = "ID", all =T)#287926     12

exon.mut.agg.cpg <- merge(exon.mut.agg.cpg, hg19_exons.df, by.x =0
                          , by.y = "ID", all =T)#315668     12
exon.mut.agg.non.cpg <- merge(exon.mut.agg.non.cpg, hg19_exons.df, by.x =0
                              , by.y = "ID", all =T)#315668     12

enh.mut.agg.cpg <- merge(enh.mut.agg.cpg, enhancers.df, by.x =0
                         , by.y = "ID", all =T)#28295    12
enh.mut.agg.non.cpg <- merge(enh.mut.agg.non.cpg, enhancers.df, by.x =0
                             , by.y = "ID", all =T)#28295    12

#As now I am taking all sequences and not only sequences with mutations, I need to
#replace NAs with zeros

interg.mut.agg.cpg[is.na(interg.mut.agg.cpg)] <- 0
interg.mut.agg.non.cpg[is.na(interg.mut.agg.non.cpg)] <- 0
exon.mut.agg.cpg[is.na(exon.mut.agg.cpg)] <- 0
exon.mut.agg.non.cpg[is.na(exon.mut.agg.non.cpg)] <- 0
enh.mut.agg.cpg[is.na(enh.mut.agg.cpg)] <- 0
enh.mut.agg.non.cpg[is.na(enh.mut.agg.non.cpg)] <- 0

#####
## making GR objects:
#if I sum 1 to start position, I reduce width by 1
interg.mut.agg.cpg_gr <- with(interg.mut.agg.cpg, GRanges(seqnames , IRanges(start,  end)))
interg.mut.agg.non.cpg_gr <- with(interg.mut.agg.non.cpg, GRanges(seqnames , IRanges(start,  end)))

exon.mut.agg.cpg_gr <- with(exon.mut.agg.cpg, GRanges(seqnames , IRanges(start,  end)))
exon.mut.agg.non.cpg_gr <- with(exon.mut.agg.non.cpg, GRanges(seqnames , IRanges(start,  end)))

enh.mut.agg.cpg_gr <- with(enh.mut.agg.cpg, GRanges(seqnames , IRanges(  start,  end)))
enh.mut.agg.non.cpg_gr <- with(enh.mut.agg.non.cpg, GRanges(seqnames , IRanges(  start,  end)))

rt_gr <- with(rt, GRanges( chr , IRanges( start+1, end)))

### OVERLAP REGIONS AND RT
#Note: quintiles can make duplicated entries that is why I am using unique()
x <- as.data.frame(findOverlaps(interg.mut.agg.cpg_gr, rt_gr))
df.inter.cpg <- cbind(interg.mut.agg.cpg[x$queryHits,], rt[x$subjectHits,"quintile",drop=F])
df.inter.cpg <- unique(df.inter.cpg)#211349     13

x <- as.data.frame(findOverlaps(interg.mut.agg.non.cpg_gr, rt_gr))
df.inter.non.cpg <- cbind(interg.mut.agg.non.cpg[x$queryHits,], rt[x$subjectHits,"quintile",drop=F])
df.inter.non.cpg <- unique(df.inter.non.cpg)# 211349     13

x <- as.data.frame(findOverlaps(exon.mut.agg.cpg_gr, rt_gr))
df.exons.cpg <- cbind(exon.mut.agg.cpg[x$queryHits,], rt[x$subjectHits,"quintile",drop=F])
df.exons.cpg <- unique(df.exons.cpg)#96348    13

x <- as.data.frame(findOverlaps(exon.mut.agg.non.cpg_gr, rt_gr))
df.exons.non.cpg <- cbind(exon.mut.agg.non.cpg[x$queryHits,], rt[x$subjectHits,"quintile",drop=F])
df.exons.non.cpg <- unique(df.exons.non.cpg)#96348    13

x <- as.data.frame(findOverlaps(enh.mut.agg.cpg_gr, rt_gr))
df.enh.cpg <- cbind(enh.mut.agg.cpg[x$queryHits,], rt[x$subjectHits,"quintile",drop=F])
df.enh.cpg <- unique(df.enh.cpg)#27503    13

x <- as.data.frame(findOverlaps(enh.mut.agg.non.cpg_gr, rt_gr))
df.enh.non.cpg <- cbind(enh.mut.agg.non.cpg[x$queryHits,], rt[x$subjectHits,"quintile",drop=F])
df.enh.non.cpg <- unique(df.enh.non.cpg)#27503    13

write.table(x = df.inter.cpg, file = "intergenic_cpg_with_RT", sep = '\t', quote = F)
write.table(x = df.inter.non.cpg, file = "inter_non_cpg_with_RT", sep = '\t', quote = F)
write.table(x = df.exons.cpg, file = "exons_cpg_with_RT", sep = '\t', quote = F)
write.table(x = df.exons.non.cpg, file = "exons_non_cpg_with_RT", sep = '\t', quote = F)
write.table(x = df.enh.cpg, file = "enh_cpg_with_RT", sep = '\t', quote = F)
write.table(x = df.enh.non.cpg, file = "enh_non_cpg_with_RT", sep = '\t', quote = F)


interg.nuc <- read.table(file = "hg19_intergenic_reduced.nuc", header =F
                         , stringsAsFactors=F, sep = '\t')
exon.nuc <- read.table(file = "hg19_exons_reduced.nuc", header =F
                       , stringsAsFactors=F, sep = '\t')
enh.nuc <- read.table(file = "hg19_villar_liver_enhancers_reduced.nuc", header =F
                      , stringsAsFactors=F, sep = '\t')

colnames(interg.nuc) <- c("usercol_1", "usercol_2", "usercol_3"
                          , "usercol_4", "usercol_5", "usercol_6", "pct_at"
                          , "pct_gc", "num_A", "num_C", "num_G", "num_T"
                          , "num_N", "num_oth", "seq_len")

colnames(exon.nuc) <- c("usercol_1", "usercol_2", "usercol_3", "usercol_4", "usercol_5"
                        , "usercol_6", "pct_at", "pct_gc", "num_A", "num_C"
                        , "num_G", "num_T", "num_N", "num_oth", "seq_len")

colnames(enh.nuc) <- c("usercol_1", "usercol_2", "usercol_3", "usercol_4"
                       , "usercol_5", "usercol_6", "pct_at", "pct_gc", "num_A", "num_C"
                       , "num_G", "num_T", "num_N", "num_oth", "seq_len")

#I need to correct the nucleotide composition so that I am using the ancestral genome
#background composition.
#For this, I will use the mutations mapped at both cpg and non cpg sites

interg.cpg <- read.table(file = "intergenic_human_anc_mut_CpG", sep = '\t'
                         , header =T, stringsAsFactors = F)#
interg.non.cpg <- read.table(file = "intergenic_human_anc_mut_non.CpG", sep = '\t'
                             , header =T, stringsAsFactors = F)#

exon.cpg <- read.table(file = "exon_human_anc_mut_CpG", sep = '\t'
                       , header =T, stringsAsFactors = F)#
exon.non.cpg <- read.table(file = "exon_human_anc_mut_non.CpG", sep = '\t'
                           , header =T, stringsAsFactors = F)#

enh.cpg <- read.table(file = "liverEnh_human_anc_mut_CpG", sep = '\t'
                      , header =T, stringsAsFactors = F)#
enh.non.cpg <- read.table(file = "liverEnh_human_anc_mut_non.CpG", sep = '\t'
                          , header =T, stringsAsFactors = F)#

#sum mutations at CpG and nonCpG sites by mutation type
interg.all <- merge(interg.cpg, interg.non.cpg, by=0, all=T)
exon.all <- merge(exon.cpg, exon.non.cpg, by=0, all=T)
enh.all <- merge(enh.cpg, enh.non.cpg, by=0, all=T)

#There are NAs in the tables
interg.all[is.na(interg.all)] <- 0
exon.all[is.na(exon.all)] <- 0
enh.all[is.na(enh.all)] <- 0

sum_mutations <- function(x) {
  x$CA <- with(x, CA.x + CA.y)
  x$CG <- with(x, CG.x + CG.y)
  x$CT <- with(x, CT.x + CT.y)
  x$TA <- with(x, TA.x + TA.y)
  x$TC <- with(x, TC.x + TC.y)
  x$TG <- with(x, TG.x + TG.y)
  return(x[,c("Row.names", "CA", "CG", "CT", "TA", "TC", "TG")])}

interg.all <- sum_mutations(interg.all)#203501      7
exon.all <- sum_mutations(exon.all)#100069      7
enh.all <- sum_mutations(enh.all)#26315     7

#number of GC nucleotides in the ancestor (lost in human)
#number of AT nucleotides in the ancestor (lost in human)

interg.all$GC_only.anc <- with(interg.all, CA+CT)#
interg.all$AT_only.anc <- with(interg.all, TC+TG)

exon.all$GC_only.anc <- with(exon.all, CA+CT)#
exon.all$AT_only.anc <- with(exon.all, TC+TG)

enh.all$GC_only.anc <- with(enh.all, CA+CT)#
enh.all$AT_only.anc <- with(enh.all, TC+TG)

interg.nuc <- merge(interg.nuc, interg.all[,c("Row.names", "GC_only.anc", "AT_only.anc")]
                    , by.x = "usercol_6", by.y="Row.names", all=T)#203501     17
# table(exon.all$Row.names %in% exon.nuc$usercol_6)#all T
exon.nuc <- merge(exon.nuc, exon.all[,c("Row.names", "GC_only.anc", "AT_only.anc")]
                  , by.x = "usercol_6", by.y="Row.names", all=T)#100069     17
# table(enh.all$Row.names %in% enh.nuc$usercol_6)#all T
enh.nuc <- merge(enh.nuc, enh.all[,c("Row.names", "GC_only.anc", "AT_only.anc")]
                 , by.x = "usercol_6", by.y="Row.names", all=T)#26315    17
#.all dataframes contain all the mutations but they don't contain the sequences without mutations
#that is why i need to keep them all and change NAs to zeros:
interg.nuc[is.na(interg.nuc)] <- 0
exon.nuc[is.na(exon.nuc)] <- 0
enh.nuc[is.na(enh.nuc)] <- 0

#calculate ancestral GC and AT frequencies
#GC ancestral frequency will be:
#human C freq + human G freq + ancestral GC lost in humans - ancestral AT that are GC in humans
interg.nuc$ancestral_GC_freq <- with(interg.nuc, num_C + num_G + GC_only.anc - AT_only.anc)
#AT ancestral frequency will be:
#human A freq + human T freq + ancestral AT lost in humans - ancestral GC that are AT in humans
interg.nuc$ancestral_AT_freq <- with(interg.nuc, num_A + num_T - GC_only.anc + AT_only.anc)

exon.nuc$ancestral_GC_freq <- with(exon.nuc, num_C + num_G + GC_only.anc - AT_only.anc)
exon.nuc$ancestral_AT_freq <- with(exon.nuc, num_A + num_T - GC_only.anc + AT_only.anc)

enh.nuc$ancestral_GC_freq <- with(enh.nuc, num_C + num_G + GC_only.anc - AT_only.anc)
enh.nuc$ancestral_AT_freq <- with(enh.nuc, num_A + num_T - GC_only.anc + AT_only.anc)

#save ancestral nucleotide content
write.table(x = interg.nuc, file = "intergenic_ancestral_nuc", sep = '\t', quote = F)
write.table(x = exon.nuc, file = "exons_ancestral_nuc", sep = '\t', quote = F)
write.table(x = enh.nuc, file = "enhancer_ancestral_nuc", sep = '\t', quote = F)
#continue:
df.inter.cpg <- read.table(file = "intergenic_cpg_with_RT", sep = '\t', header=T, stringsAsFactors = F)
df.inter.non.cpg <- read.table(file = "inter_non_cpg_with_RT", sep = '\t', header=T, stringsAsFactors = F)
df.exons.cpg <- read.table(file = "exons_cpg_with_RT", sep = '\t', header=T, stringsAsFactors = F)
df.exons.non.cpg <- read.table(file = "exons_non_cpg_with_RT", sep = '\t', header=T, stringsAsFactors = F)
df.enh.cpg <- read.table(file = "enh_cpg_with_RT", sep = '\t', header=T, stringsAsFactors = F)
df.enh.non.cpg <- read.table(file = "enh_non_cpg_with_RT", sep = '\t', header=T, stringsAsFactors = F)

interg.nuc <- read.table(file = "intergenic_ancestral_nuc", sep = '\t', header=T, stringsAsFactors = F)
exon.nuc <- read.table(file = "exons_ancestral_nuc", sep = '\t', header=T, stringsAsFactors = F)
enh.nuc <- read.table(file = "enhancer_ancestral_nuc", sep = '\t', header=T, stringsAsFactors = F)


df.inter.cpg <- merge(df.inter.cpg
                      , interg.nuc[,c("usercol_6", "ancestral_GC_freq", "ancestral_AT_freq")]
                      , by.x = "Row.names", by.y = "usercol_6", all.x=T)#211349     15
df.inter.non.cpg <- merge(df.inter.non.cpg
                          , interg.nuc[,c("usercol_6", "ancestral_GC_freq", "ancestral_AT_freq")]
                          , by.x = "Row.names", by.y = "usercol_6", all.x=T)#211349     15
df.exons.cpg <- merge(df.exons.cpg
                      , exon.nuc[,c("usercol_6", "ancestral_GC_freq", "ancestral_AT_freq")]
                      , by.x = "Row.names", by.y = "usercol_6", all.x=T)# 96348    15
df.exons.non.cpg <- merge(df.exons.non.cpg
                          , exon.nuc[,c("usercol_6", "ancestral_GC_freq", "ancestral_AT_freq")]
                          , by.x = "Row.names", by.y = "usercol_6", all.x=T)#96348    15
df.enh.cpg <- merge(df.enh.cpg
                    , enh.nuc[,c("usercol_6", "ancestral_GC_freq", "ancestral_AT_freq")]
                    , by.x = "Row.names", by.y = "usercol_6", all.x=T)#27503    15
df.enh.non.cpg <- merge(df.enh.non.cpg
                        , enh.nuc[,c("usercol_6", "ancestral_GC_freq", "ancestral_AT_freq")]
                        , by.x = "Row.names", by.y = "usercol_6", all.x=T)#27503    15


df.inter.cpg$has_Ns <- ifelse(with(df.inter.cpg, ancestral_GC_freq+ancestral_AT_freq) != df.inter.cpg$width, "yes", "no")
df.exons.cpg$has_Ns <- ifelse(with(df.exons.cpg, ancestral_GC_freq+ancestral_AT_freq) != df.exons.cpg$width, "yes", "no")
df.enh.cpg$has_Ns <- ifelse(with(df.enh.cpg, ancestral_GC_freq+ancestral_AT_freq) != df.enh.cpg$width, "yes", "no")
df.inter.non.cpg$has_Ns <- ifelse(with(df.inter.non.cpg, ancestral_GC_freq+ancestral_AT_freq) != df.inter.non.cpg$width, "yes", "no")
df.exons.non.cpg$has_Ns <- ifelse(with(df.exons.non.cpg, ancestral_GC_freq+ancestral_AT_freq) != df.exons.non.cpg$width, "yes", "no")
df.enh.non.cpg$has_Ns <- ifelse(with(df.enh.non.cpg, ancestral_GC_freq+ancestral_AT_freq) != df.enh.non.cpg$width, "yes", "no")

df.inter.cpg <- df.inter.cpg[df.inter.cpg$has_Ns == "no",]
df.exons.cpg <- df.exons.cpg[df.exons.cpg$has_Ns == "no",]
df.enh.cpg <- df.enh.cpg[df.enh.cpg$has_Ns == "no",]

df.inter.non.cpg <- df.inter.non.cpg[df.inter.non.cpg$has_Ns == "no",]
df.exons.non.cpg <- df.exons.non.cpg[df.exons.non.cpg$has_Ns == "no",]
df.enh.non.cpg <- df.enh.non.cpg[df.enh.non.cpg$has_Ns == "no",]

#normalise GC loss and gain mutations by ancestral background
norm_by_ancestral_nuc <- function(x){
  x$GC_loss <- with(x, CA + CT)  
  x$GC_gain <- with(x, TC + TG)
  #aggregate by quintile
  x$GCloss_norm <- with(x, GC_loss/ancestral_GC_freq)
  x$GCgain_norm <- with(x, GC_gain/ancestral_AT_freq)
 
  return(x)}

#### normalize mutation counts by ancestral background:
df.inter.cpg <- norm_by_ancestral_nuc(df.inter.cpg)
df.inter.non.cpg <- norm_by_ancestral_nuc(df.inter.non.cpg)

df.exons.cpg <- norm_by_ancestral_nuc(df.exons.cpg)
df.exons.non.cpg <- norm_by_ancestral_nuc(df.exons.non.cpg)

df.enh.cpg <- norm_by_ancestral_nuc(df.enh.cpg)
df.enh.non.cpg <- norm_by_ancestral_nuc(df.enh.non.cpg)

df.inter.cpg$region <- df.inter.non.cpg$region <- "intergenic"
df.exons.cpg$region <- df.exons.non.cpg$region <- "exons"
df.enh.cpg$region <- df.enh.non.cpg$region <- "enhancers"

df.inter.cpg$type <-
  df.exons.cpg$type <-
  df.enh.cpg$type <-  "CpG"

df.inter.non.cpg$type <-
  df.exons.non.cpg$type <-
  df.enh.non.cpg$type <-"nonCpG"

df <- rbind(df.inter.cpg[,c("GCloss_norm", "GCgain_norm", "quintile", "region", "type")]
            , df.inter.non.cpg[,c("GCloss_norm", "GCgain_norm", "quintile", "region", "type")]
            , df.exons.cpg[,c("GCloss_norm", "GCgain_norm", "quintile", "region", "type")]
            , df.exons.non.cpg[,c("GCloss_norm", "GCgain_norm", "quintile", "region", "type")]
            , df.enh.cpg[,c("GCloss_norm", "GCgain_norm", "quintile", "region", "type")]
            , df.enh.non.cpg[,c("GCloss_norm", "GCgain_norm", "quintile", "region", "type")])#669710      5


df$quintile <- as.factor(df$quintile)
df$type <- as.factor(df$type)

df <- melt(df)

table(is.na(df$value))
# FALSE    TRUE
# 1339404      16

df <- df[!is.na(df$value),]
summary(df$value)

se <- function(x) sqrt(var(x) / length(x))

df$log10_value <- log10(df$value + 1)
df.sub <- df[,c("quintile", "region", "type", "variable", "log10_value")]
df.se <- aggregate(log10_value~., df.sub, se)#60  5
colnames(df.se)[5] <- "se"
df <- aggregate(log10_value~., df.sub, mean)#
df <- merge(df, df.se, by=c("quintile", "region", "type", "variable"))
df$region <- as.factor(df$region)

pdf("human_bonobo_substitutions_CpG_vs_nonCpG_ancNuc_se.pdf")
ggplot(data=df, aes(x=quintile, y=log10_value, group =region, linetype=region)) +
  geom_line() + #geom_point() +
  geom_errorbar(aes(ymin=log10_value-se, ymax=log10_value+se), width=.2,
                position=position_dodge(0.05)) +
  theme_classic() + facet_wrap(.~variable+type)
dev.off()

