library(ggplot2)

setwd("./replication_timing/human")
liver <- read.table(file = "human_liver_all_enh_rt_humanJaspar_heat"
                    , header = T, stringsAsFactors = F, sep = '\t')
annot <- read.table(file = "jaspar_motifs_TFs_annot_all.txt", header =T
                    , stringsAsFactors = F, sep = '\t')
liver <- merge(liver, annot[,c("ID", "Name")], by.x = 0, by.y = "ID"
               , all.x =T)

liver$Row.names <- paste(liver$Row.names, liver$Name, sep = '_')
rownames(liver) <- liver$Row.names
liver.sub <- liver[,c("earlyvlate", "Class", "GC.proportion")]
liver.sub <- liver.sub[with(liver.sub, order(Class, earlyvlate)),]

# mean earlyvlate by class
mean.rt.by.class <- aggregate(earlyvlate~Class, liver.sub, mean)
mean.rt.by.class <- mean.rt.by.class[with(mean.rt.by.class, order(earlyvlate)),]
colnames(mean.rt.by.class) <- c("Class", "mean_earlyvlate")
liver.sub$motif_id <- rownames(liver.sub)
liver.sub <- merge(liver.sub, mean.rt.by.class, by ="Class")
liver.sub <- liver.sub[with(liver.sub
                            , order(mean_earlyvlate, earlyvlate)),]

rownames(liver.sub) <- liver.sub$motif_id
liver.sub$motif_id <- NULL


#removing classes with less than 10 motifs
x <- as.data.frame(table(liver.sub$Class))
x <- x[x$Freq>10,]
liver.sub <- liver.sub[liver.sub$Class %in% x$Var1,]

liver.sub$id <- 1:nrow(liver.sub)

pdf("human_liver_earlyvlate_humanJaspar_by_TFClass_BAR.pdf")
ggplot(data = liver.sub, aes(x = id, y = GC.proportion, fill = Class)) +
  geom_bar(stat="identity", width=0.5) + theme_classic()
dev.off()
