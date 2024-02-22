setwd("/replication_timing/human")
human.jaspar <- read.table(file = "human_jaspar_tfs_and_TFBSs_RT.txt", header =T
                           , stringsAsFactors = F, sep = '\t')
human.jaspar <- unique(human.jaspar[, c("Name", "gene_id", "earlyvlate")])

quintile <- function(obs){
  qn <- quantile(obs, probs = (0:5)/5)
  result <- as.numeric(cut(obs, qn, include.lowest = T))
  return(result)}


human.jaspar$earlyvlate_bin <- quintile(human.jaspar$earlyvlate)
table( human.jaspar$earlyvlate_bin)
# 1   2   3   4   5
# 102 101 101 101 102
human.expr1 <- read.table(file = "HumFagerbergTScomparisonTable_9_27.txt", header =T
                          , stringsAsFactors = F)

table(human.jaspar$gene_id %in% human.expr1$Ensembl.Gene.ID)

setwd("./replication_timing/human")
human.jaspar.tau <- merge(human.jaspar, human.expr1[,c("Ensembl.Gene.ID", "Tau")]
                          , by.x = "gene_id", by.y = "Ensembl.Gene.ID")


library(ggplot2)
human.jaspar.tau$earlyvlate_bin <- factor(human.jaspar.tau$earlyvlate_bin, levels = c(1,2,3,4,5))
pdf("HumFagerberg_tau_vs_H9_RT_humanJaspar_motifs.pdf")
ggplot(human.jaspar.tau, aes(x=earlyvlate_bin, y=Tau)) +
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black")
dev.off()
length(unique(genes_rt.tau1$gene_id))

table(human.jaspar.tau$earlyvlate_bin)
#  1   2   3   4   5
# 91  93  96  97 100

for(i in 1:4){
  print(paste("test: quintile", i, "vs quintile", i+1))
  x <- wilcox.test(human.jaspar.tau[human.jaspar.tau$earlyvlate_bin ==i, "Tau"]
                   , human.jaspar.tau[human.jaspar.tau$earlyvlate_bin ==(i+1), "Tau"]
                   , alternative = "greater")
  print(x$p.value)
}
# [1] "test: quintile 1 vs quintile 2"
# [1] 0.3454562
# [1] "test: quintile 2 vs quintile 3"
# [1] 7.146121e-06
# [1] "test: quintile 3 vs quintile 4"
# [1] 0.002677551
# [1] "test: quintile 4 vs quintile 5"
# [1] 0.9171333
