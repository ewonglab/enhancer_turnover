library(GenomicRanges)
library(ggplot2)

setwd("./replication_timing/mouse")
all_enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                      , header = F, stringsAsFactors = F, sep ='\t')
active_atLeast2sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast2sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
cons_enh <- rbind(active_atLeast2sp, poised_atLeast2sp)
cons_enh$V1 <- paste0("chr", cons_enh$V1)
unique(sub(":.*", "", cons_enh$V1))
# [1] "chr1"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17"
# [10] "chr18" "chr19" "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"
# [19] "chr9"  "chrX"

all_enh <- all_enh[all_enh$V6 == "Recent" | all_enh$V4 %in% cons_enh$V1, ] 

rt <- read.table(file = "mouse_sperm_pgc_meanRT.bedGraph"
                 , header = F, stringsAsFactors = F, sep = '\t')
rt <- rt[complete.cases(rt),]
colnames(rt)[4] <- "rt"

# making gr objects
all_enh_gr <- with(all_enh, GRanges(V1, IRanges(V2 + 1, V3)))
rt.gr <- with(rt, GRanges(V1, IRanges(V2 + 1, V3)))

x <- as.data.frame(findOverlaps(all_enh_gr, rt.gr))
all_enh <- cbind(all_enh[x$queryHits,]
                 , rt[x$subjectHits, "rt", drop =F])

# aggregate RT
all_enh <- aggregate(rt ~. , all_enh, mean)

## add relative RT and separate into early RT and late RT
colnames(all_enh) <- c("chr", "start", "end", "ID"
                       , "mark", "type", "tissue","rt")
all_enh$rel_rt <- ""
all_enh$rel_rt <- ifelse(all_enh$rt > 0.5, "early", all_enh$rel_rt)
all_enh$rel_rt <- ifelse(all_enh$rt < -0.5, "late", all_enh$rel_rt)


#keeping only early and late enhancers

all_enh <- all_enh[all_enh$rel_rt %in% c("early", "late"),]

# separate betwen tissue specific and non-tissue specific
all_enh$tissue_sp <- ifelse(grepl(",", all_enh$tissue)
                            , "Non_tissue_specffic", "tissue_specific")

table(all_enh$tissue_sp)
# Non_tissue_specffic     tissue_specific
#               19094               78188

as.data.frame(table(all_enh[,c("tissue_sp", "rel_rt")]))
#             tissue_sp rel_rt  Freq
# 1 Non_tissue_specffic  early 18323
# 2     tissue_specific  early 67976
# 3 Non_tissue_specffic   late   771
# 4     tissue_specific   late 10212

# fisher test
m <- matrix(data = c(67976, 18323, 10212, 771), nrow = 2, ncol = 2, byrow = T)
colnames(m) <- c("Tissue-specific", "Non-tissue-specific")
rownames(m) <- c("early", "Late")


x <- fisher.test(m)
x

# Fisher's Exact Test for Count Data
#
# data:  m
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.2595220 0.3019607
# sample estimates:
# odds ratio
#  0.2800732

x$p.value
# [1] 0
x$estimate
# odds ratio
# 0.2800732

m.df <- as.data.frame(m)
m.df$rt <- rownames(m.df)
colnames(m.df) <- gsub("-", "_", colnames(m.df))

# all enhancers

pdf("mouse_enh_non.tissue.spec_piechart_cons2sp.pdf")
ggplot(m.df, aes(x = "", y = Non_tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()

pdf("mouse_enh_tissue.spec_piechart_cons2sp.pdf")
ggplot(m.df, aes(x = "", y = Tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()


# separate between conserved and recent enhancers
cons <- all_enh[all_enh$type == "Conserved",] 
recent <- all_enh[all_enh$type == "Recent",] 

# CONSERVED
as.data.frame(table(cons[,c("tissue_sp", "rel_rt")]))
#       tissue_sp rel_rt  Freq
# 1 Non_tissue_specffic  early 12953
# 2     tissue_specific  early 37674
# 3 Non_tissue_specffic   late   390
# 4     tissue_specific   late  2469

m <- matrix(data = c(37674, 12953, 2469, 390), nrow = 2, ncol = 2, byrow = T)
colnames(m) <- c("Tissue-specific", "Non-tissue-specific")
rownames(m) <- c("early", "Late")


x <- fisher.test(m)
x
# Fisher's Exact Test for Count Data
#
# data:  m
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.4110689 0.5124200
# sample estimates:
# odds ratio
#  0.4594308
#  
x$p.value
# [1] 1.758916e-52
x$estimate
# odds ratio
# 0.4594308

m.df <- as.data.frame(m)
m.df$rt <- rownames(m.df)
colnames(m.df) <- gsub("-", "_", colnames(m.df))

pdf("mouse_enh_non.tissue.spec_CONS_piechart_cons2sp.pdf")
ggplot(m.df, aes(x = "", y = Non_tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()

pdf("mouse_enh_tissue.spec_CONS_piechart_cons2sp.pdf")
ggplot(m.df, aes(x = "", y = Tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()


## RECENT
as.data.frame(table(recent[,c("tissue_sp", "rel_rt")])) 
#             tissue_sp rel_rt  Freq
# 1 Non_tissue_specffic  early  5370
# 2     tissue_specific  early 30302
# 3 Non_tissue_specffic   late   381
# 4     tissue_specific   late  7743

m <- matrix(data = c(30302, 5370, 7743, 381), nrow = 2, ncol = 2, byrow = T)
colnames(m) <- c("Tissue-specific", "Non-tissue-specific")
rownames(m) <- c("early", "Late")

x <- fisher.test(m)
x
#
# Fisher's Exact Test for Count Data
#
# data:  m
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.2488702 0.3090943
# sample estimates:
# odds ratio
#  0.2776401
#  
x$p.value
# [1] 1.428166e-167
x$estimate
# odds ratio
# 0.2776401

m.df <- as.data.frame(m)
m.df$rt <- rownames(m.df)
colnames(m.df) <- gsub("-", "_", colnames(m.df))

pdf("mouse_enh_non.tissue.spec_RECENT_piechart.pdf")
ggplot(m.df, aes(x = "", y = Non_tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()

pdf("mouse_enh_tissue.spec_RECENT_piechart.pdf")
ggplot(m.df, aes(x = "", y = Tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()


# percentages

# all enhancers

#             tissue_sp rel_rt  Freq
# 1 Non_tissue_specffic  early 18323
# 2     tissue_specific  early 67976
# 3 Non_tissue_specffic   late   771
# 4     tissue_specific   late 10212


# non-tissue-spec
round((18323/(18323+771))*100, 2)
# 95.96
round((771/(18323+771))*100, 2)
# 4.04

# tissue-spec
round((67976/(67976+10212))*100, 2)
# 86.94
round((10212/(67976+10212))*100, 2)
# 13.06

# conserved
#       tissue_sp rel_rt  Freq
# 1 Non_tissue_specffic  early 12953
# 2     tissue_specific  early 37674
# 3 Non_tissue_specffic   late   390
# 4     tissue_specific   late  2469

# non-tissue-specific
round((12953/(12953+390))*100, 2) # early
# 97.08
round((390/(12953+390))*100, 2) # late
# 2.92

# tissue-specific
round((37674/(37674+2469))*100, 2) # early
# 93.85
round((2469/(37674+2469))*100, 2) # late
# 6.15

# recent
# 1 Non_tissue_specffic  early  5370
# 2     tissue_specific  early 30302
# 3 Non_tissue_specffic   late   381
# 4     tissue_specific   late  7743

# non-tissue-specific
round((5370/(5370+381))*100, 2) # early
# 93.38
round((381/(5370+381))*100, 2) # late
# 6.62

# tissue-specific
round((30302/(30302+7743))*100, 2) # early
# 79.65
round((7743/(30302+7743))*100, 2) # late
# 20.35
