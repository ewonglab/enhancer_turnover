library(GenomicRanges)

setwd("./replication_timing/mouse")
mouse_rt <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header =T
                       , stringsAsFactors = F, sep = '\t')
somatic <- c("X46C", "D3", "TT2", "iPSC", "iPSC_1D4"
             , "iPSC_2D4", "EPL", "EBM3", "EpiSC5", "EpiSC7"
             , "EBM6", "X46CNPC", "TT2NPC", "EBM9", "Mesoderm"
             , "Endoderm", "piPSC_1A2", "piPSC_1B3", "piPSC_V3"
             , "MEF_female", "MEF_male", "Myoblast")
germline <- c("sperm.1", "sperm.2", "PGC.male.1", "PGC.female.1")

#mean RT
mouse_rt$mean_germlineRT <- apply(mouse_rt[,germline], 1, mean)
mouse_rt$mean_somaticRT <- apply(mouse_rt[,somatic], 1, mean)

#reading enhancers
mouse.enh <- read.table(file = "./roller/mouse_all_enh_byMark_type_and_tissue.bed"
                        , header = F, stringsAsFactors = F
                        , sep = "\t")
recent_enh <- mouse.enh[mouse.enh$V6 == "Recent",]
# conserved in min 4 other species
active_atLeast4sp <- read.table(file = "./roller/active_conserved_atLeast4sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
poised_atLeast4sp <- read.table(file = "./roller/poised_conserved_atLeast4sp.txt"
                                , header = T, stringsAsFactors = F, sep = "\t")
cons_enh <- rbind(active_atLeast4sp, poised_atLeast4sp)

#add coordinates to rt data
mouse_rt$peakid <- rownames(mouse_rt)
mouse_rt$chr <- gsub(":.*", "", mouse_rt$peakid)
mouse_rt$start <- gsub("_.*", "", gsub(".*:", "", mouse_rt$peakid))
mouse_rt$end <- gsub(".*_", "", gsub(".*:", "", mouse_rt$peakid))

mouse_rt$start <- as.integer(mouse_rt$start)
mouse_rt$end <- as.integer(mouse_rt$end)

# separate conserved enhancers coordinates
cons_enh$chr <- paste0("chr", sub(":.*", "", cons_enh$V1))
cons_enh$start <- sub("-.*", "", sub(".*:", "", cons_enh$V1))
cons_enh$end <- sub(".*-", "", sub(".*:", "", cons_enh$V1))
cons_enh$start <- as.integer(cons_enh$start)
cons_enh$end <- as.integer(cons_enh$end)

#making gr objects
mouse_rt.gr <- with(mouse_rt, GRanges(chr, IRanges(start+1, end)))
cons_enh.gr <- with(cons_enh, GRanges(chr, IRanges(start+1, end)))
recent_enh.gr <- with(recent_enh, GRanges(V1, IRanges(V2+1, V3)))

#overlap enhancers with bins and calculate number of enhancers per bin
x <- as.data.frame(findOverlaps(mouse_rt.gr, cons_enh.gr))
bins_with_consEnh <- cbind(mouse_rt[x$queryHits,"peakid",drop =F]
                           , cons_enh[x$subjectHits,"V1", drop =F]) # V1 = ID

x <- as.data.frame(findOverlaps(mouse_rt.gr, recent_enh.gr))
bins_with_recentEnh <- cbind(mouse_rt[x$queryHits,"peakid",drop =F]
                             , recent_enh[x$subjectHits,"V4", drop =F])# V4 = ID

# add kmeans cluster
bins_with_consEnh <- merge(bins_with_consEnh, mouse_rt[,c("peakid", "kmeans.cluster")]
                           , by = "peakid", all.x =T)
bins_with_recentEnh <- merge(bins_with_recentEnh, mouse_rt[,c("peakid", "kmeans.cluster")]
                             , by = "peakid", all.x =T)

#add cluster number
#count unique number of enhancers overlapping every cluster
n_consEnh <- lapply(1:18, function(x) length(unique(bins_with_consEnh[bins_with_consEnh$kmeans.cluster == x, "V1"])))
n_recentEnh <- lapply(1:18, function(x) length(unique(bins_with_recentEnh[bins_with_recentEnh$kmeans.cluster == x, "V4"])))

# mean RT per kmeans cluster
mean_rt <- aggregate(.~kmeans.cluster
                     , mouse_rt[,c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT")], mean)
mean_rt <- mean_rt[with(mean_rt, order(kmeans.cluster)),]
mean_rt$n_conserved <- unlist(n_consEnh)
mean_rt$n_recent <- unlist(n_recentEnh)

mean_rt$recent_vs_cons_LogFC <- with(mean_rt, log(n_recent/n_conserved))

mean_rt
# kmeans.cluster mean_germlineRT mean_somaticRT n_conserved n_recent
# 1               1       0.2557718     -0.3377223        1769     3731
# 2               2      -1.5276388     -0.6201144          83      720
# 3               3       0.5344305      0.5811477        1921     3350
# 4               4       1.1588596      1.5392447        7700     7525
# 5               5      -0.2522994     -0.7056573        1116     2682
# 6               6       1.0143551      1.3136903        7410     7329
# 7               7      -1.3777489     -0.8931174         274     2946
# 8               8      -0.3946769     -0.7273448         897     2884
# 9               9       0.5564480      0.5738714        1684     2419
# 10             10       0.7491716      0.6339254        2302     2805
# 11             11       0.8332980      1.0137935        2739     4149
# 12             12      -0.6998528     -1.0728049         795     4109
# 13             13       0.3041689      0.1175160        1631     3135
# 14             14       0.2708385     -0.1796346        1675     2637
# 15             15       0.1021640     -0.8683291        1516     3045
# 16             16       0.8287007      1.0084502        3627     4985
# 17             17       0.4266030      0.2032035        2010     2953
# 18             18      -0.1720332     -0.3165335         692     2288
# recent_vs_cons_LogFC
# 1            0.74626188
# 2            2.16041060
# 3            0.55611446
# 4           -0.02298952
# 5            0.87681192
# 6           -0.01099136
# 7            2.37507549
# 8            1.16787764
# 9            0.36218232
# 10           0.19762523
# 11           0.41527445
# 12           1.64259285
# 13           0.65343585
# 14           0.45382874
# 15           0.69742561
# 16           0.31802754
# 17           0.38468688
# 18           1.19584740

library(ggplot2)
pdf("mouse_germlineRT_vs_LogFCrecent_v_cons4sp.pdf")
ggplot(data=mean_rt, aes_string(x="mean_germlineRT", y="recent_vs_cons_LogFC", group=1)) +
  geom_point() + theme_classic() + geom_smooth(method=lm)
dev.off()

pdf("mouse_somaticRT_vs_LogFCrecent_v_cons4sp.pdf")
ggplot(data=mean_rt, aes_string(x="mean_somaticRT", y="recent_vs_cons_LogFC", group=1)) +
  geom_point() + theme_classic() + geom_smooth(method=lm)
dev.off()

x <- cor.test(mean_rt$mean_germlineRT, mean_rt$recent_vs_cons_LogFC)
x$p.value
#[1] 4.196598e-12
x$estimate
# cor
# -0.9766027
round(x$estimate^2, 2)
# 0.95

x <- cor.test(mean_rt$mean_somaticRT, mean_rt$recent_vs_cons_LogFC)
# x$p.value
# [1] 0.0001298478
x$estimate
# cor
# -0.7811051
round(x$estimate^2, 2)
# cor
# 0.61
