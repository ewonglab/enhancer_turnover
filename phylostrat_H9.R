#the table jaspar_tfs_and_TFBSs_RT.txt with TFs replication time was produced using the script jaspar_human_liver_RT_correlat.R
# that table is used along with the data from Litman and Stein 2019 (Obtaining estimates for the ages of all the protein-coding genes and most of the ontology-identified noncoding genes of the human genome, assigned to 19 phylostrata) to make violin plots of the TFs divided by replication time quintiles.

#Note: the TFs replication time was calculated for all the TFs in jaspar 2020 that overlap with H9 rep time bedgraph.
#The TFs were NOT filtered to keep only those from Human
# Tfs with duplicated names were removed (like those profiled in more than one species) 

#code to produce violin plots
phylo <- read.csv(file = "fixed_up_no_19.3s_phylo.csv", header =T
                  , stringsAsFactors = F, skip =2)
phylo <- phylo[!is.na(phylo$modal.value),]#19653    35

quintile <- function(obs){
  qn <- quantile(obs, probs = (0:5)/5)
  result <- as.numeric(cut(obs, qn, include.lowest = T))
  return(result)}


tf_rt <- read.table(file = "jaspar_tfs_and_TFBSs_RT.txt"
                     , stringsAsFactors = F, sep = '\t')
tf_rt <- unique(tf_rt[, c("Name", "gene_id", "earlyvlate")])# 520   3

tf_rt$earlyvlate_bin <- quintile(tf_rt$earlyvlate)

table(tf_rt$Name %in% phylo$HGNC)
# FALSE  TRUE 
# 87   433 
tf_rt <- merge(tf_rt, unique(phylo[,c("HGNC", "modal.value")])
               , by.x = "Name", by.y = "HGNC")#433   5
tf_rt$earlyvlate_bin <- factor(tf_rt$earlyvlate_bin, levels = c(1,2,3,4,5))

library(ggplot2)
pdf("phylostrat_jaspar.pdf")
ggplot(tf_rt, aes(x=earlyvlate_bin, y=modal.value)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black")
dev.off()

#######

#phylostrat
library(ggplot2)

phylo <- read.csv(file = "fixed_up_no_19.3s_phylo.csv", header =T
                  , stringsAsFactors = F, skip =2)
phylo <- phylo[!is.na(phylo$modal.value),]#19653    35

quintile <- function(obs){
  qn <- quantile(obs, probs = (0:5)/5)
  result <- as.numeric(cut(obs, qn, include.lowest = T))
  return(result)}


tf_rt <- read.table(file = "all_jaspar_tfs_and_TFBSs_RT.txt"
                    , stringsAsFactors = F, sep = '\t')
tf_rt <- unique(tf_rt[, c("Name", "gene_id", "earlyvlate")])# 681   4
tf_rt$earlyvlate_bin <- quintile(tf_rt$earlyvlate)
#number of TFs in every replication time bin
table(tf_rt$earlyvlate_bin)
# 1   2   3   4   5 
# 102 101 101 101 102 

table(tf_rt$Name %in% phylo$HGNC)
# FALSE  TRUE 
# 148   533 
tf_rt <- merge(tf_rt, unique(phylo[,c("HGNC", "modal.value")])
               , by.x = "Name", by.y = "HGNC")#533   5
tf_rt$earlyvlate_bin <- factor(tf_rt$earlyvlate_bin, levels = c(1,2,3,4,5))
pdf("phylostrat_all_jaspar.pdf")
ggplot(tf_rt, aes(x=earlyvlate_bin, y=modal.value)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black")
dev.off()


tf_rt <- read.table(file = "human_jaspar_tfs_and_TFBSs_RT.txt"
                    , stringsAsFactors = F, sep = '\t')
tf_rt <- unique(tf_rt[, c("Name", "gene_id", "earlyvlate")])# 507   3
tf_rt$earlyvlate_bin <- quintile(tf_rt$earlyvlate)
table(tf_rt$Name %in% phylo$HGNC)
# TRUE 
# 507
tf_rt <- merge(tf_rt, unique(phylo[,c("HGNC", "modal.value")])
               , by.x = "Name", by.y = "HGNC")#533   5
tf_rt$earlyvlate_bin <- factor(tf_rt$earlyvlate_bin, levels = c(1,2,3,4,5))

pdf("phylostrat_human_jaspar.pdf")
ggplot(tf_rt, aes(x=earlyvlate_bin, y=modal.value)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black") +
  scale_y_continuous(limits = c(-5, 25), breaks = seq(-5, 25, by = 1))
dev.off()

for(i in 1:4){
  print(paste("quintiles", i, "vs", (i+1)))
  x <- wilcox.test(tf_rt[tf_rt$earlyvlate_bin==i,"modal.value"]
                   , tf_rt[tf_rt$earlyvlate_bin==(i+1),"modal.value"], alternative = "greater")
  print(x$p.value)
}

# [1] "quintiles 1 vs 2"
# [1] 0.1379856
# [1] "quintiles 2 vs 3"
# [1] 0.01441811
# [1] "quintiles 3 vs 4"
# [1] 0.7328193
# [1] "quintiles 4 vs 5"
# [1] 0.3191132
