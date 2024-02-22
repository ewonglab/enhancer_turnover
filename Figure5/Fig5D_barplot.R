library(ggplot2)
setwd("./replication_timing/human")

heatmap_data <- read.table(file = "human_liver_all_enh_rt_humanJaspar_heat"
                           , header = T, stringsAsFactors = F, sep ='\t') 

# confirm order
heatmap_data <- heatmap_data[with(heatmap_data, order(earlyvlate)),]
# unique id
heatmap_data$id <- 1:nrow(heatmap_data)

range(heatmap_data$GC.proportion)
# [1] 0.0000000 0.9090909

pdf("barplot.earlylate.humanLiver_allEnh_humanJas.pdf")
ggplot(data = heatmap_data, aes(x = id, y = GC.proportion)) +
  geom_bar(stat="identity", width=0.5) + theme_classic()
dev.off()
