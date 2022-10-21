library(GenomicRanges)
library(GenomicFeatures, lib = '/scratch/zk16/cc3704/R_packages')

setwd("/g/data/zk16/cc3704/replication_timing/human")
#get hg19 genes
human_txdb <- makeTxDbFromGFF(file="/g/data/zk16/cc3704/human_data/hg19.ensGene.gtf"
                              , format="gtf",organism="Homo sapiens")

hg19_genes <- genes(human_txdb)
hg19_genes.df <- as.data.frame(hg19_genes)

# read DNA replication time predicted from ovary ATAC-seq
ovary_rt <- fread('RT_H9_ESC_Ext29405702_hg19.bedgraph', header=F)
# create GenomicRanges object
ovary_rt_gr <- with(ovary_rt, GRanges( V1 , IRanges( V2+1, V3 )))
names(ovary_rt) <- c('chr','start','end','rt')

# read DNA replication time predicted from testis ATAC-seq
testis_rt <- fread('RT_H9_ESC_Ext29405702_hg19.bedgraph', header=F)
# create GenomicRanges object
testis_rt_gr <- with(testis_rt, GRanges( V1 , IRanges( V2+1, V3 )))
names(testis_rt) <- c('chr','start','end','rt')

#mean replication time per gene
x <- as.data.frame(findOverlaps(hg19_genes, ovary_rt_gr))
genes_ovary.rt <- cbind(hg19_genes.df[x$queryHits,]
                        , ovary_rt[x$subjectHits,"rt", drop =F])
genes_ovary.rt <- aggregate(rt~., genes_ovary.rt, mean)

x <- as.data.frame(findOverlaps(hg19_genes, testis_rt_gr))
genes_testis.rt <- cbind(hg19_genes.df[x$queryHits,]
                         , testis_rt[x$subjectHits,"rt", drop =F])
genes_testis.rt <- aggregate(rt~., genes_testis.rt, mean)

quintile <- function(obs){
  qn <- quantile(obs, probs = (0:5)/5)
  result <- as.numeric(cut(obs, qn, include.lowest = T))
  return(result)}

# assign DNA replication time quintile
genes_ovary.rt$rt_quintile <- quintile(genes_ovary.rt$rt)
genes_testis.rt$rt_quintile <- quintile(genes_testis.rt$rt)

#save
write.table(x = genes_ovary.rt, file = "hg19_genes_with_ovaryReplicon_RT.txt"
            , quote = F, sep = '\t', row.names = F)
write.table(x = genes_testis.rt, file = "hg19_genes_with_testisReplicon_RT.txt"
            , quote = F, sep = '\t', row.names = F)


# Plot tau scores (tissue specificity scores) separated by predicted germline DNA replication time quintiles

human.jaspar <- read.table(file = "human_jaspar_tfs_and_TFBSs_RT.txt", header =T
                           , stringsAsFactors = F, sep = '\t')
                           
setwd("/g/data/zk16/cc3704/replication_timing/tau/Expression_data")

# read tau scores
human.expr1 <- read.table(file = "HumFagerbergTScomparisonTable_9_27.txt", header =T
                          , stringsAsFactors = F)

# reading genes with DNA replication time
setwd("/g/data/zk16/cc3704/replication_timing/human")
genes_ovary.rt <- read.table(file = "hg19_genes_with_ovaryReplicon_RT.txt", header =T
                       , stringsAsFactors = F, sep = '\t')
genes_testis.rt <- read.table(file = "hg19_genes_with_testisReplicon_RT.txt", header =T
                       , stringsAsFactors = F, sep = '\t')

#subset genes to human transcription factors
genes_ovary.rt <- genes_ovary.rt[genes_ovary.rt$gene_id %in% human.jaspar$gene_id,]
genes_testis.rt <- genes_testis.rt[genes_testis.rt$gene_id %in% human.jaspar$gene_id,]

genes_ovary.rt.tau1 <- merge(genes_ovary.rt, human.expr1[,c("Ensembl.Gene.ID", "Tau")]
                       , by.x = "gene_id", by.y = "Ensembl.Gene.ID")
genes_testis.rt.tau1 <- merge(genes_testis.rt, human.expr1[,c("Ensembl.Gene.ID", "Tau")]
                       , by.x = "gene_id", by.y = "Ensembl.Gene.ID")

library(ggplot2)
genes_ovary.rt.tau1$rt_quintile <- factor(genes_ovary.rt.tau1$rt_quintile, levels = c(1,2,3,4,5))
genes_testis.rt.tau1$rt_quintile <- factor(genes_testis.rt.tau1$rt_quintile, levels = c(1,2,3,4,5))

pdf("HumFagerberg_tau_vs_ovaryReplicon.RT_humanJaspar_TFs.pdf")
ggplot(genes_ovary.rt.tau1, aes(x=rt_quintile, y=Tau)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black")
dev.off()

pdf("HumFagerberg_tau_vs_testisReplicon.RT_humanJaspar_TFs.pdf")
ggplot(genes_testis.rt.tau1, aes(x=rt_quintile, y=Tau)) + 
  geom_violin(trim=FALSE) + theme_classic() +
  stat_summary(fun=median, geom="point", size=2, color="black")
dev.off()

#Number of TFs in every dna replication time quintile

table(genes_ovary.rt.tau1$rt_quintile)
# 1   2   3   4   5 
# 54  67 113 113 130

table(genes_testis.rt.tau1$rt_quintile)
# 1   2   3   4   5 
# 54  67 113 113 130 
