data.4.heatmap_jas <- function(comp.df, enh_type){
  enh_cols <- paste(c('early_motif_non_zeros', 'early_zero',  'late_motif_non_zeros', 'late_zero')
                    , enh_type, sep = "_")
  df<- comp.df[,enh_cols]
  row.names(df) <- comp.df$pwm
  df <- df[complete.cases(df), ]
  print(nrow(df))
  comp.df$earlyvlate <- log2((comp.df[,enh_cols[1]] / comp.df[,enh_cols[2]]) /
                               (comp.df[,enh_cols[3]] / comp.df[,enh_cols[4]]) )
  comp.df <- merge(comp.df, df, by.x='pwm', by.y=0)#
  #plot earlyvlate with % AT %GC
  df <- comp.df[,c('earlyvlate','C.proportion','G.proportion'
                   ,'GC.proportion','A.proportion','T.proportion','AT.proportion'
                   , "Family", "Class", "Species")]
  row.names(df) <- comp.df$pwm
  df <- df[is.finite(rowSums(df[,c(1:7)])),]#
  print(nrow(df))
  df <- df[order(df$earlyvlate), ]
  return(df)
}

comp <- read.table(file = "motif_comp_mouse_enh_ov_germRT_allTissues_humanJAS_cons2sp.txt"
                   , header = T, stringsAsFactors = F, sep = '\t')
### ALL
all.heat <-  data.4.heatmap_jas(comp, "all")


write.table(x = all.heat
            , file = "mouse_all_tissues_germ.rt_humanJAS_cons2sp_heat"
            , quote = F, sep = '\t')

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

setwd("./replication_timing/mouse")
mouse <- read.table(file = "mouse_all_tissues_germ.rt_humanJAS_cons2sp_heat"
                    , header = T, stringsAsFactors = F, sep = '\t')# already ordered by earlyvlate

steps <- rev(brewer.pal(n = 8, name = "PuOr"))
pal <- color.palette(steps, c(20, 20,5,5,5, 20,20), space="rgb")

# all enhancers
pdf('heatmap.earlylate.mouse_allTissues_humanJAS_cons2sp.pdf')#
pheatmap(mouse[,c(2:7)], scale='column'
         , show_rownames = FALSE
         , color = pal(500)
         , show_colnames=TRUE
         ,  main="Mouse all tissues earlyvlate",
         cluster_rows=FALSE, cluster_cols=F, fontsize=5)
dev.off()
