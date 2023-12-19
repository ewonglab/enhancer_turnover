library(GenomicRanges)
library(ggplot2)

all_enh <- read.table(file = "mouse_all_enh_byMark_type_and_tissue.bed"
                      , header = F, stringsAsFactors = F, sep ='\t')#175011      7

rt <- read.table(file = "mouse_sperm_pgc_meanRT.bedGraph"
                 , header = F, stringsAsFactors = F, sep = '\t')
rt <- rt[complete.cases(rt),]#20809     4
colnames(rt)[4] <- "rt"

#making gr objects
all_enh_gr <- with(all_enh, GRanges(V1, IRanges(V2 + 1, V3)))#175011
rt.gr <- with(rt, GRanges(V1, IRanges(V2 + 1, V3)))#20809

x <- as.data.frame(findOverlaps(all_enh_gr, rt.gr))
all_enh <- cbind(all_enh[x$queryHits,]
                 , rt[x$subjectHits, "rt", drop =F])#[1] 142503      8
length(unique(all_enh$V4))
# [1] 140755

# aggregate RT
all_enh <- aggregate(rt ~. , all_enh, mean) # 140755      8

## add relative RT and separate into early RT and late RT
colnames(all_enh) <- c("chr", "start", "end", "ID"
                       , "mark", "type", "tissue","rt")
all_enh$rel_rt <- ""
all_enh$rel_rt <- ifelse(all_enh$rt > 0.5, "early", all_enh$rel_rt)
all_enh$rel_rt <- ifelse(all_enh$rt < -0.5, "late", all_enh$rel_rt)

table(all_enh$rel_rt)
#       early  late 
# 39646 89921 11188 

#keeping only early and late enhancers

all_enh <- all_enh[all_enh$rel_rt %in% c("early", "late"),]#101109      9
range(all_enh[all_enh$rel_rt=="early", "rt"])#0.5000533 1.9178578
range(all_enh[all_enh$rel_rt=="late", "rt"])#-2.238118 -0.500722

# separate betwen tissue specific and non-tissue specific
all_enh$tissue_sp <- ifelse(grepl(",", all_enh$tissue)
                       , "Non_tissue_specffic", "tissue_specific")

table(all_enh$tissue_sp)
# Non_tissue_specffic     tissue_specific 
#               20095               81014

as.data.frame(table(all_enh[,c("tissue_sp", "rel_rt")]))
#             tissue_sp rel_rt  Freq
# 1 Non_tissue_specffic  early 19290
# 2     tissue_specific  early 70631
# 3 Non_tissue_specffic   late   805
# 4     tissue_specific   late 10383

# fisher test
m <- matrix(data = c(70631, 19290, 10383, 805), nrow = 2, ncol = 2, byrow = T)
colnames(m) <- c("Tissue-specific", "Non-tissue-specific")
rownames(m) <- c("early", "Late")
m
#       Tissue-specific Non-tissue-specific
# early           70631               19290
# Late            10383                 805


x <- fisher.test(m)
x
# 
# Fisher's Exact Test for Count Data
# 
# data:  m
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.2634527 0.3055594
# sample estimates:
# odds ratio 
#  0.2838599 
#  
x$p.value
# [1] 0
x$estimate
# odds ratio 
# 0.2838599 

m.df <- as.data.frame(m)
m.df$rt <- rownames(m.df)
colnames(m.df) <- gsub("-", "_", colnames(m.df))

pdf("mouse_enh_non.tissue.spec_piechart.pdf")
ggplot(m.df, aes(x = "", y = Non_tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()

pdf("mouse_enh_tissue.spec_piechart.pdf")
ggplot(m.df, aes(x = "", y = Tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()


# separate between conserved and recent enhancers
cons <- all_enh[all_enh$type == "Conserved",] # 57313    10
recent <- all_enh[all_enh$type == "Recent",] # 43796    10

# CONSERVED
as.data.frame(table(cons[,c("tissue_sp", "rel_rt")]))
#             tissue_sp rel_rt  Freq
# 1 Non_tissue_specffic  early 13920
# 2     tissue_specific  early 40329
# 3 Non_tissue_specffic   late   424
# 4     tissue_specific   late  2640

m <- matrix(data = c(40329, 13920, 2640, 424), nrow = 2, ncol = 2, byrow = T)
colnames(m) <- c("Tissue-specific", "Non-tissue-specific")
rownames(m) <- c("early", "Late")
m
#       Tissue-specific Non-tissue-specific
# early           40329               13920
# Late             2640                 424


x <- fisher.test(m)
x
# Fisher's Exact Test for Count Data
# 
# data:  m
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.4182206 0.5167204
# sample estimates:
# odds ratio 
#  0.4653128 
#  
x$p.value
# [1] 1.088807e-54
x$estimate
# odds ratio 
# 0.4653128

m.df <- as.data.frame(m)
m.df$rt <- rownames(m.df)
colnames(m.df) <- gsub("-", "_", colnames(m.df))

pdf("mouse_enh_non.tissue.spec_CONS_piechart.pdf")
ggplot(m.df, aes(x = "", y = Non_tissue_specific, fill = rt)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void()
dev.off()

pdf("mouse_enh_tissue.spec_CONS_piechart.pdf")
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
m
#       Tissue-specific Non-tissue-specific
# early           30302                5370
# Late             7743                 381

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
