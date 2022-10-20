library(data.table)
library(stringr)

setwd("/g/data/zk16/cc3704/replication_timing/human")

# read predicted replication time by chromosome
replicon_aml_files <- list.files(path = "./AML_hg19_IPLS"
                                     , pattern = "(.*)timing$", full.names = T)
replicon_aml <- lapply(replicon_aml_files, fread)
names(replicon_aml) <- basename(replicon_aml_files)

for(i in 1:length(replicon_aml)){
  replicon_aml[[i]]$chr <- names(replicon_aml[i])}

replicon_aml <- do.call("rbind", replicon_aml)#6186215       3

replicon_aml$chr <- gsub("\\.timing", "", replicon_aml$chr)
replicon_aml$end <- replicon_aml$V1 + 500
replicon_aml <- replicon_aml[,c("chr", "V1", "end", "V2")]

#multiply by -1 to match experimental RT values direction
replicon_aml$V2 <- (replicon_aml$V2) * (-1)
replicon_aml <- replicon_aml[with(replicon_aml, order(chr, V1)),]

write.table(x = replicon_aml, file = "replicon_pHSC_hg19.bed"
            , col.names = F, row.names = F, sep = '\t', quote = F)

