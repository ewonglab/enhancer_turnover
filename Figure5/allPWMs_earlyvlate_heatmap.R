library(pheatmap)                  
library(RColorBrewer)

color.palette <- function(steps, n.steps.between=NULL, ...){
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}

liver <- read.table(file = "human_liver_all_enh_rt_humanJaspar_heat"
                    , header = T, stringsAsFactors = F, sep = '\t')

annot <- read.table(file = "jaspar_motifs_TFs_annot_all.txt", header =T
                    , stringsAsFactors = F, sep = '\t')

liver <- merge(liver, annot[,c("ID", "Name")], by.x = 0, by.y = "ID"
               , all.x =T)

liver$Row.names <- paste(liver$Row.names, liver$Name, sep = '_')
rownames(liver) <- liver$Row.names
liver.sub <- liver[,c("earlyvlate", "Class")]
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

# showing some motifs as Row names
motifs <- c("MA0442.1_SOX10", "MA0866.1_SOX21", "MA0769.2_TCF7"
            , "MA1421.1_TCF7L1", "MA0077.1_SOX9", "MA0442.2_SOX10"
            , "MA0143.4_SOX2", "MA0868.2_SOX8", "MA0867.2_SOX4"
            , "MA1563.1_SOX18", "MA1120.1_SOX13", "MA1152.1_SOX15"
            , "MA1562.1_SOX14", "MA1478.1_DMRTA2", "MA0610.1_DMRT3"
            , "MA1479.1_DMRTC2", "MA0498.2_MEIS1", "MA1507.1_HOXD4"
            , "MA1504.1_HOXC4", "MA0597.1_THAP1", "MA1573.1_THAP11"
            , "MA1561.1_SOX12", "MA0084.1_SRY")

liver.sub$tmp <- 1:nrow(liver.sub)

for(i in 1:nrow(liver.sub)){
  if(rownames(liver.sub)[i] %in% motifs){
    liver.sub$tmp[i] <- rownames(liver.sub)[i]}}

rownames(liver.sub) <- liver.sub$tmp
liver.sub$tmp <- NULL

#removing classes with less than 10 motifs

x <- as.data.frame(table(liver.sub$Class))
x <- x[x$Freq>10,]#13
liver.sub <- liver.sub[liver.sub$Class %in% x$Var1,] 

Class_df <- data.frame("TF Class" = liver.sub$Class)
rownames(Class_df) <- rownames(liver.sub)

steps <- rev(brewer.pal(n = 8, name = "RdYlBu"))
pal <- color.palette(rev(steps), c(20, 20,5,5,5, 20,20), space="rgb")

pdf("human_liver_earlyvlate_humanJaspar_by_TFClass.5.pdf")
pheatmap(liver.sub[,2, drop =F], scale='column', annotation_row = Class_df
         , show_colnames=TRUE, show_rownames = TRUE
         ,  main="Human liver earlyvlate by TF class",
         , color = pal(100)
         , cluster_rows=FALSE, cluster_cols=FALSE, fontsize=5)
dev.off()
