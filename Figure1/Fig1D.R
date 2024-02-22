
gilbert <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header =T
                      , stringsAsFactors = F, sep = '\t')

#calculating mean RT by sample cluster
sample.clust1 <- apply(gilbert[,c("X46C", "D3", "TT2", "iPSC", "iPSC_1D4"
                                  , "iPSC_2D4", "EPL")], 1, mean)
sample.clust2 <-  apply(gilbert[,c("EBM3", "EpiSC5", "EpiSC7")], 1, mean)
sample.clust3 <- apply(gilbert[,c("EBM6", "X46CNPC", "TT2NPC", "EBM9" )], 1, mean)
sample.clust4 <- apply(gilbert[,c("Mesoderm", "Endoderm")], 1, mean)
sample.clust5 <- apply(gilbert[,c("piPSC_1A2", "piPSC_1B3", "piPSC_V3")], 1, mean)
sample.clust6 <- apply(gilbert[,c("MEF_female", "MEF_male", "Myoblast")], 1, mean)
sample.clust7 <- gilbert$PGC.female.1
sample.clust8 <- gilbert$PGC.male.1
sample.clust9 <- apply(gilbert[,c("sperm.1", "sperm.2")], 1, mean)

samples.clust.annot <- data.frame(cluster1 = sample.clust1
                                  , cluster2 = sample.clust2
                                  , cluster3 = sample.clust3
                                  , cluster4 =sample.clust4
                                  , cluster5 = sample.clust5
                                  , cluster6 =sample.clust6
                                  , cluster7 = sample.clust7
                                  , cluster8 = sample.clust8
                                  , cluster9 = sample.clust9
                                  , stringsAsFactors = F)#8966    9

early <- apply(samples.clust.annot, 1, function(x) all(x > 0.5))
late <- apply(samples.clust.annot, 1, function(x) all(x < -0.5))

table(early)
# early
# FALSE  TRUE
# 7599  1367
              
round((1367*100)/nrow(gilbert), 2)# 15.25

table(late)
# late
# FALSE  TRUE
# 7397  1569

round((1569*100)/nrow(gilbert),2)# 17.5

gilbert$peak_id <- rownames(gilbert)
gilbert$chr <- gsub(":.*", "", gilbert$peak_id)
gilbert$start <- gsub("_.*", "", gsub(".*:", "", gilbert$peak_id))
gilbert$end <- gsub(".*_", "", gsub(".*:", "", gilbert$peak_id))
gilbert$start <- as.numeric(gilbert$start)
gilbert$end <- as.numeric(gilbert$end)
gilbert$width <- with(gilbert, end-start)

always_early <- gilbert[rownames(gilbert) %in% names(early[early==TRUE]), ]
always_late <- gilbert[rownames(gilbert) %in% names(late[late==TRUE]), ]

round((sum(always_early$width)*100)/sum(gilbert$width), 2)# 13.94
round((sum(always_late$width)*100)/sum(gilbert$width), 2)# 19.14

# constant regions
all_constantRT <- rbind(always_early, always_late) 
              
# define variable regions
variableRT <- gilbert[!rownames(gilbert) %in% rownames(all_constantRT), ] 


# read enhancers and calculate the proportion of those enhancers overlapping
#reading enhancers
mouse.enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                        , header = F, stringsAsFactors = F
                        , sep = "\t")#
# keep only recent enhancers #
recent_enh <- mouse.enh[mouse.enh$V6 == "Recent",]# 80904     7
recent_enh$V1 <- sub("chr", "", recent_enh$V1)
main_chr <- c("1", "10", "11", "12", "13", "14", "15", "16", "17", "18"
              , "19", "2", "3", "4", "5", "6", "7", "8", "9", "X", "Y")

recent_enh$V1 <- ifelse(recent_enh$V1 %in% main_chr, paste0("chr", recent_enh$V1) 
                        , recent_enh$V1)

# MAKE GR OBJECTS
library(GenomicRanges)

recent_enh.gr <- with(recent_enh, GRanges(V1, IRanges(V2+1, V3)))
all_constantRT.gr <- with(all_constantRT, GRanges(chr, IRanges(start+1, end)))
always_early.gr <- with(always_early, GRanges(chr, IRanges(start+1, end)))
always_late.gr <- with(always_late, GRanges(chr, IRanges(start+1, end)))
variableRT.gr <- with(variableRT, GRanges(chr, IRanges(start+1, end)))

# count number of enhancers overlapping each set of bins

# RECENT ENHANCERS IN CONSTANT BINS
x <- as.data.frame(findOverlaps(all_constantRT.gr, recent_enh.gr))
recent_ov_constant <- recent_enh[unique(x$subjectHits), "V4"]

# RECENT ENHANCERS IN CONSTANTLY EARLY BINS
x <- as.data.frame(findOverlaps(always_early.gr, recent_enh.gr))
recent_ov_constE <- recent_enh[unique(x$subjectHits), "V4"]

# RECENT ENHANCERS IN CONSTANTLY LATE BINS
x <- as.data.frame(findOverlaps(always_late.gr, recent_enh.gr))
recent_ov_constL <- recent_enh[unique(x$subjectHits), "V4"]

# RECENT ENHANCERS IN VARIABLE RT BINS
x <- as.data.frame(findOverlaps(variableRT.gr, recent_enh.gr))
recent_ov_variable <- recent_enh[unique(x$subjectHits), "V4"]

# for each category remove enhancers overlapping other categories
# only plot constantly E and T
recent_dup <- c(recent_ov_constE, recent_ov_constL, recent_ov_variable)
recent_dup <- recent_dup[duplicated(recent_dup)] #

recent_ov_constE <- recent_ov_constE[!recent_ov_constE %in% recent_dup]
recent_ov_constL <- recent_ov_constL[!recent_ov_constL %in% recent_dup]
recent_ov_variable <- recent_ov_variable[!recent_ov_variable %in% recent_dup]

## OVERLAP OF CONSERVED ENHANCERS
active_atLeast2sp <- read.table(file = "./roller/active_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast2sp <- read.table(file = "./roller/poised_conserved_atLeast2sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
cons_enh <- rbind(active_atLeast2sp, poised_atLeast2sp)

cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end)

unique(cons_enh$chr)
# [1] "chr1"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17"
# [10] "chr18" "chr19" "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"
# [19] "chr9"  "chrX"

cons_enh.gr <- with(cons_enh, GRanges(chr, IRanges(start+1, end)))#

# CONSERVED ENHANCERS IN CONSTANT BINS
x <- as.data.frame(findOverlaps(all_constantRT.gr, cons_enh.gr))
cons_ov_constant <- cons_enh[unique(x$subjectHits), "V1"]

# CONSERVED ENHANCERS IN CONSTANTLY EARLY BINS
x <- as.data.frame(findOverlaps(always_early.gr, cons_enh.gr))
cons_ov_constE <- cons_enh[unique(x$subjectHits), "V1"]

# CONSERVED ENHANCERS IN CONSTANTLY LATE BINS
x <- as.data.frame(findOverlaps(always_late.gr, cons_enh.gr))
cons_ov_constL <- cons_enh[unique(x$subjectHits), "V1"]

# CONSERVED ENHANCERS IN VARIABLE RT BINS
x <- as.data.frame(findOverlaps(variableRT.gr, cons_enh.gr))
cons_ov_variable <- cons_enh[unique(x$subjectHits), "V1"]

# only plot consE, consL and variable
cons_dup <- c(cons_ov_constE, cons_ov_constL, cons_ov_variable)
cons_dup <- cons_dup[duplicated(cons_dup)] #

cons_ov_constE <- cons_ov_constE[!cons_ov_constE %in% cons_dup]# 22042
cons_ov_constL <- cons_ov_constL[!cons_ov_constL %in% cons_dup]# 968
cons_ov_variable <- cons_ov_variable[!cons_ov_variable %in% cons_dup] # 48150

frequency_df <- data.frame(consE = c(13478, 22042)
                           , consL = c(3978, 968)
                           , Var = c(46005, 48150)
                           , type = c("recent", "conserved"))
library(reshape2)
library(viridis)
frequency_df <- melt(frequency_df, "type")

frequency_df$total <- ifelse(frequency_df$type == "conserved",
                             sum(frequency_df[frequency_df$type == "conserved", "value"])
                             , sum(frequency_df[frequency_df$type == "recent", "value"]))
frequency_df
#        type variable value total
# 1    recent    consE 13478 63461
# 2 conserved    consE 22042 71160
# 3    recent    consL  3978 63461
# 4 conserved    consL   968 71160
# 5    recent      Var 46005 63461
# 6 conserved      Var 48150 71160

# calculate percentage
frequency_df$perc <- with(frequency_df, (value/total)*100)
frequency_df
# type variable value total      perc
# 1    recent    consE 13478 63461 21.238241
# 2 conserved    consE 22042 71160 30.975267
# 3    recent    consL  3978 63461  6.268417
# 4 conserved    consL   968 71160  1.360315
# 5    recent      Var 46005 63461 72.493342
# 6 conserved      Var 48150 71160 67.664418

library(ggplot2)

pdf("Enh_RT_constant_v_rariable_nonDup.pdf")
ggplot(frequency_df, aes(fill=variable, y=perc, x=type)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) + theme_classic() +
  # theme_ipsum() +
  xlab("")
dev.off()
