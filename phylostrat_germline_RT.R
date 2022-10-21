library(ggplot2)

# reading phylostrat data
setwd("/g/data/zk16/cc3704/replication_timing/human/cres_with_replicon")
phylo <- read.csv(file = "/g/data/zk16/cc3704/replication_timing/human/fixed_up_no_19.3s_phylo.csv", header =T
                  , stringsAsFactors = F, skip =2)
phylo <- phylo[!is.na(phylo$modal.value),]

quintile <- function(obs){
  qn <- quantile(obs, probs = (0:5)/5)
  result <- as.numeric(cut(obs, qn, include.lowest = T))
  return(result)}

# read TFs with replication time 
# DNA replication time predicted from ovary ATAC-seq

tf_ovary.rt <- read.table(file = "human_jaspar_tfs_and_TFBSs_ovary.RT.txt"
                    , stringsAsFactors = F, sep = '\t')
tf_ovary.rt <- unique(tf_ovary.rt[, c("Name", "gene_id", "earlyvlate")])

# defining quintiles of TF binding motifs replication time  
tf_ovary.rt$earlyvlate_bin <- quintile(tf_ovary.rt$earlyvlate)
table(tf_ovary.rt$earlyvlate_bin)
# 1   2   3   4   5 
# 102 101 101 101 101

tf_ovary.rt <- merge(tf_ovary.rt, unique(phylo[,c("HGNC", "modal.value")])
               , by.x = "Name", by.y = "HGNC")
tf_ovary.rt$earlyvlate_bin <- factor(tf_ovary.rt$earlyvlate_bin, levels = c(1,2,3,4,5))

# NUmber of TFs per DNA replication quintiles
table(tf_ovary.rt$earlyvlate_bin)
# 1   2   3   4   5 
# 102 101 101 101 101 

pdf("phylostrat_human.jaspar_ovaryRT.pdf")
ggplot(tf_ovary.rt, aes(x=earlyvlate_bin, y=modal.value)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black") +
  scale_y_continuous(breaks = seq(0, 20, by = 1))
dev.off()

# Fisher's test for the difference in phylostrat

for(i in 1:4){
  print(paste("quintiles", i, "vs", (i+1)))
  x <- wilcox.test(tf_ovary.rt[tf_ovary.rt$earlyvlate_bin==i,"modal.value"]
                   , tf_ovary.rt[tf_ovary.rt$earlyvlate_bin==(i+1),"modal.value"]
                   , alternative = "greater")
  print(x$p.value)
}

# [1] "quintiles 1 vs 2"
# [1] 0.206608

# [1] "quintiles 2 vs 3"
# [1] 0.003449166

# [1] "quintiles 3 vs 4"
# [1] 0.7253045

# [1] "quintiles 4 vs 5"
# [1] 0.5930718

# DNA replication time predicted from testis ATAC-seq

phylo <- read.csv(file = "/g/data/zk16/cc3704/replication_timing/human/fixed_up_no_19.3s_phylo.csv", header =T
                  , stringsAsFactors = F, skip =2)
phylo <- phylo[!is.na(phylo$modal.value),]


tf_testis.rt <- read.table(file = "human_jaspar_tfs_and_TFBSs_testis.RT.txt"
                          , stringsAsFactors = F, sep = '\t')#506  24
tf_testis.rt <- unique(tf_testis.rt[, c("Name", "gene_id", "earlyvlate")])# 506   3

tf_testis.rt$earlyvlate_bin <- quintile(tf_testis.rt$earlyvlate)
table(tf_testis.rt$earlyvlate_bin)
# 1   2   3   4   5 
# 103 102 111  93 103

tf_testis.rt <- merge(tf_testis.rt, unique(phylo[,c("HGNC", "modal.value")])
                     , by.x = "Name", by.y = "HGNC")
tf_testis.rt$earlyvlate_bin <- factor(tf_testis.rt$earlyvlate_bin, levels = c(1,2,3,4,5))

# Number of TFs per DNA replication time quintile

table(tf_testis.rt$earlyvlate_bin)
# 1   2   3   4   5 
# 103 102 111  93 103

pdf("phylostrat_human.jaspar_testisRT.pdf")
ggplot(tf_testis.rt, aes(x=earlyvlate_bin, y=modal.value)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black") + 
  scale_y_continuous(breaks = seq(0, 20, by = 1))
dev.off()

# Fisher's exact test for the difference in phylostrat

for(i in 1:4){
  print(paste("quintiles", i, "vs", (i+1)))
  x <- wilcox.test(tf_testis.rt[tf_testis.rt$earlyvlate_bin==i,"modal.value"]
                   , tf_testis.rt[tf_testis.rt$earlyvlate_bin==(i+1),"modal.value"], alternative = "greater")
  print(x$p.value)}

# [1] "quintiles 1 vs 2"
# [1] 0.1149558

# [1] "quintiles 2 vs 3"
# [1] 0.06995738

# [1] "quintiles 3 vs 4"
# [1] 0.4641701

# [1] "quintiles 4 vs 5"
# [1] 0.8793054
