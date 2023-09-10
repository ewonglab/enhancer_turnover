# Script used to produce Figure 2D and calculate R2 values between germ line and somatic mean DNA replication time and enhancer turnover across clusters in Fig. 2C

library(GenomicRanges)
library(ggplot2)

mouse_rt <- read.table(file = "hiratani_plus_germline_18kmeans.txt", header =T
                       , stringsAsFactors = F, sep = '\t')
somatic <- c("X46C", "D3", "TT2", "iPSC", "iPSC_1D4"
             , "iPSC_2D4", "EPL", "EBM3", "EpiSC5", "EpiSC7"
             , "EBM6", "X46CNPC", "TT2NPC", "EBM9", "Mesoderm"
             , "Endoderm", "piPSC_1A2", "piPSC_1B3", "piPSC_V3"
             , "MEF_female", "MEF_male", "Myoblast")
germline <- c("sperm.1", "sperm.2", "PGC.male.1", "PGC.female.1")

# mean RT
mouse_rt$mean_germlineRT <- apply(mouse_rt[,germline], 1, mean)
mouse_rt$mean_somaticRT <- apply(mouse_rt[,somatic], 1, mean)

#reading enhancers
mouse.enh <- read.table(file = "mouse_all_enh_byMark_type_and_tissue.bed", header = F, stringsAsFactors = F, sep = "\t")#
table(mouse.enh$V6)
# Conserved    Recent
# 94107     80904

#add coordinates to rt data
mouse_rt$peakid <- rownames(mouse_rt)
mouse_rt$chr <- gsub(":.*", "", mouse_rt$peakid)
mouse_rt$start <- gsub("_.*", "", gsub(".*:", "", mouse_rt$peakid))
mouse_rt$end <- gsub(".*_", "", gsub(".*:", "", mouse_rt$peakid))

mouse_rt$start <- as.integer(mouse_rt$start)
mouse_rt$end <- as.integer(mouse_rt$end)

#making gr objects
mouse_rt.gr <- with(mouse_rt, GRanges(chr, IRanges(start+1, end)))#8966

cons_enh <- mouse.enh[mouse.enh$V6=="Conserved", ]
recent_enh <- mouse.enh[mouse.enh$V6=="Recent", ]

cons_enh.gr <- with(cons_enh, GRanges(V1, IRanges(V2+1, V3)))#94107
recent_enh.gr <- with(recent_enh, GRanges(V1, IRanges(V2+1, V3)))#80904

#overlap enhancers with bins and calculate number of enhancers per bin
x <- as.data.frame(findOverlaps(mouse_rt.gr, cons_enh.gr))
bins_with_consEnh <- cbind(mouse_rt[x$queryHits,"peakid",drop =F]
                           , cons_enh[x$subjectHits,"V4", drop =F])

x <- as.data.frame(findOverlaps(mouse_rt.gr, recent_enh.gr))
bins_with_recentEnh <- cbind(mouse_rt[x$queryHits,"peakid",drop =F]
                             , recent_enh[x$subjectHits,"V4", drop =F])

bins_with_consEnh <- merge(bins_with_consEnh, mouse_rt[,c("peakid", "kmeans.cluster")], by = "peakid", all.x =T)
bins_with_recentEnh <- merge(bins_with_recentEnh, mouse_rt[,c("peakid", "kmeans.cluster")], by = "peakid", all.x =T)

#add cluster number
#count unique number of enhancers overlapping every cluster
n_consEnh <- lapply(1:18, function(x) length(unique(bins_with_consEnh[bins_with_consEnh$kmeans.cluster == x, "V4"])))
n_recentEnh <- lapply(1:18, function(x) length(unique(bins_with_recentEnh[bins_with_recentEnh$kmeans.cluster == x, "V4"])))

# mean RT and number of enhancers per cluster
mean_rt <- aggregate(.~kmeans.cluster
                     , mouse_rt[,c("kmeans.cluster", "mean_germlineRT", "mean_somaticRT")], mean)
mean_rt <- mean_rt[with(mean_rt, order(kmeans.cluster)),]
mean_rt$n_conserved <- unlist(n_consEnh)
mean_rt$n_recent <- unlist(n_recentEnh)

#turnover calculated as logFC CONSERVED VS RECENT FOR LINE PLOTS
mean_rt$cons_vs_recent_LogFC <- with(mean_rt, log(n_conserved/n_recent))
mean_rt
# kmeans.cluster mean_germlineRT mean_somaticRT n_conserved n_recent
# 1               1       0.2557718     -0.3377223        3513     3731
# 2               2      -1.5276388     -0.6201144         218      720
# 3               3       0.5344305      0.5811477        3809     3350
# 4               4       1.1588596      1.5392447       13333     7525
# 5               5      -0.2522994     -0.7056573        2288     2682
# 6               6       1.0143551      1.3136903       13555     7329
# 7               7      -1.3777489     -0.8931174         771     2946
# 8               8      -0.3946769     -0.7273448        1872     2884
# 9               9       0.5564480      0.5738714        3326     2419
# 10             10       0.7491716      0.6339254        4416     2805
# 11             11       0.8332980      1.0137935        5447     4149
# 12             12      -0.6998528     -1.0728049        2020     4109
# 13             13       0.3041689      0.1175160        3212     3135
# 14             14       0.2708385     -0.1796346        3429     2637
# 15             15       0.1021640     -0.8683291        3060     3045
# 16             16       0.8287007      1.0084502        6912     4985
# 17             17       0.4266030      0.2032035        3905     2953
# 18             18      -0.1720332     -0.3165335        1548     2288
# cons_vs_recent_LogFC
# 1          -0.060205921
# 2          -1.194756149
# 3           0.128406342
# 4           0.572011354
# 5          -0.158884711
# 6           0.614916402
# 7          -1.340515223
# 8          -0.432170841
# 9           0.318416149
# 10          0.453830770
# 11          0.272197657
# 12         -0.710082178
# 13          0.024264622
# 14          0.262626766
# 15          0.004914015
# 16          0.326825628
# 17          0.279436181
# 18         -0.390714298

#order of clusters:

#order of clusters:
mean_rt$kmeans.cluster <- factor(mean_rt$kmeans.cluster
                                 , levels=c(4,6,10,11,16,3,9,17,13,14,
                                            1,15,18,5,8,12,2,7))


pdf("mouse_germlineRT_lineplot_mm10_TEST.pdf")
ggplot(data=mean_rt, aes_string(y="mean_germlineRT"
                                , x="kmeans.cluster", group=1)) +
  geom_line() + theme_classic()
dev.off()

x <- cor.test(mean_rt$mean_germlineRT, mean_rt$cons_vs_recent_LogFC)
x$p.value
# [1] 3.788733e-12
# 0.9769026 * 0.9769026
# [1] 0.9543387
x$estimate
#0.9769026  #round((0.9769026)^2, 2)

x <- cor.test(mean_rt$mean_somaticRT, mean_rt$cons_vs_recent_LogFC)
x$p.value
# [1] 0.0001157407
x$estimate
#0.7845459
# 0.7845459  * 0.7845459
# [1] 0.6155123
