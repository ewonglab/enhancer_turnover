library(ggplot2)

setwd("/g/data/zk16/cc3704/replication_timing/human")

# reading human TFs 
human.jaspar <- read.table(file = "human_jaspar_tfs_and_TFBSs_RT.txt", header =T
                           , stringsAsFactors = F, sep = '\t')
                           
# reading mouse tau scores of gene expression specificity
setwd("/g/data/zk16/cc3704/replication_timing/tau/Expression_data")
mouse_mic <- read.table(file = "MusThorrezTScomparisonTable_9_19mT.txt", header =T
                        , stringsAsFactors = F)
                        
#orthologous downloaded from ensembl biomart
orth <- read.table(file = "/g/data/zk16/cc3704/useful/orthologs/mart_human_mouse.txt"
                   , header =T, stringsAsFactors = F, sep = '\t')

orth <- orth[orth$Gene.stable.ID %in% human.jaspar$gene_id,]

#keeping only genes with orthologous in mouse
mouse_mic <- merge(mouse_mic, unique(orth[,c("Gene.stable.ID", "Mouse.gene.stable.ID")])
                   , by.x = "Ensembl.Gene.ID", by.y = "Mouse.gene.stable.ID", all.x = T)

mouse_mic <- mouse_mic[mouse_mic$Gene.stable.ID %in% human.jaspar$gene_id,]
length(unique(mouse_mic$Ensembl.Gene.ID))#379 #two genes in mouse map to one gene in human #keep both 
length(unique(mouse_mic$Gene.stable.ID))#378

# add TF name to gene ids
mouse_mic.full <- merge(mouse_mic, human.jaspar[,c("gene_id", "Row.names")]
                        , by.x = "Gene.stable.ID", by.y="gene_id")#379  33

setwd("/g/data/zk16/cc3704/replication_timing/mouse")
# germline DNA replication time of mouse TF binding sites
mouse_comp <- read.table(file = "mouse_all_tissues_germ.rt_humanJaspar_heat"
                         , header =T, stringsAsFactors = F, sep = '\t')
mouse_mic.full <- merge(mouse_mic.full, mouse_comp[,"earlyvlate",drop=F]
                        , by.x="Row.names", by.y=0)

quintile <- function(obs){
  qn <- quantile(obs, probs = (0:5)/5)
  result <- as.numeric(cut(obs, qn, include.lowest = T))
  return(result)}

# Define quintiles of TFBS DNA replication time
mouse_mic.full$earlyvlate_bin <- quintile(mouse_mic.full$earlyvlate)

# Number of TFs per DNA replication time quintile
table( mouse_mic.full$earlyvlate_bin)
# 1  2  3  4  5 
# 72 72 72 72 72 

mouse_mic.full$earlyvlate_bin <- factor(mouse_mic.full$earlyvlate_bin
                                         , levels = c(1,2,3,4,5))
mouse_mic.full <- unique(mouse_mic.full)

pdf("MusThorrez_tau_vs_germline_RT_humanJaspar_motifs.pdf")
ggplot(mouse_mic.full, aes(x=earlyvlate_bin, y=Tau)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black")
dev.off()

# Fisher's exact test for the difference in tau scores

for(i in 1:4){
  print(paste("test: quintile", i, "vs quintile", i+1))
  x <- wilcox.test(mouse_mic.full[mouse_mic.full$earlyvlate_bin ==i, "Tau"]
                   , mouse_mic.full[mouse_mic.full$earlyvlate_bin ==(i+1), "Tau"]
                   , alternative = "greater")
  print(round(x$p.value, 2))
}

# [1] "test: quintile 1 vs quintile 2"
# [1] 0.27

# [1] "test: quintile 2 vs quintile 3"
# [1] 0.04

# [1] "test: quintile 3 vs quintile 4"
# [1] 0.23

# [1] "test: quintile 4 vs quintile 5"
# [1] 0.09
