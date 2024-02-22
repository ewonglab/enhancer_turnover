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

mouse <- read.table(file = "mouse_allTissues_homeodomain_topQuintile_humanJAS_heat"
                    , header = T, stringsAsFactors = F, sep = '\t')#754  10

annot <- read.table(file = "jaspar_motifs_TFs_annot_all.txt"
                    , header = T, stringsAsFactors = F, sep = '\t')

mouse <- merge(mouse, annot[,c("ID", "Name")], by.x = 0, by.y = "ID"
               , all.x =T)

mouse$Row.names <- paste(mouse$Row.names, mouse$Name, sep = '_')
rownames(mouse) <- mouse$Row.names
mouse.sub <- mouse[,c("earlyvlate", "Class")]
mouse.sub <- mouse.sub[with(mouse.sub, order(Class, earlyvlate)),]

# mean earlyvlate by class
mean.rt.by.class <- aggregate(earlyvlate~Class, mouse.sub, mean)
mean.rt.by.class <- mean.rt.by.class[with(mean.rt.by.class, order(earlyvlate)),]
colnames(mean.rt.by.class) <- c("Class", "mean_earlyvlate")
mouse.sub$motif_id <- rownames(mouse.sub)
mouse.sub <- merge(mouse.sub, mean.rt.by.class, by ="Class")
mouse.sub <- mouse.sub[with(mouse.sub
                            , order(mean_earlyvlate, earlyvlate)),]

rownames(mouse.sub) <- mouse.sub$motif_id
mouse.sub$motif_id <- NULL


Class_df <- data.frame("TF Class" = mouse.sub$Class)
rownames(Class_df) <- rownames(mouse.sub)

steps <- rev(brewer.pal(n = 8, name = "RdYlBu"))
pal <- color.palette(rev(steps), c(20, 20,5,5,5, 20,20), space="rgb")

pdf("mouse_earlyvlate_homeodomain_topQ_humanJAS_by_TFClass.pdf")
pheatmap(mouse.sub[,2, drop =F], scale='column', annotation_row = Class_df
         , show_colnames=TRUE, show_rownames = TRUE
         ,  main="Human mouse earlyvlate by TF class",
         , color = pal(100)
         , cluster_rows=FALSE, cluster_cols=FALSE, fontsize=5)
dev.off()
