library(data.table)
library(GenomicRanges)
library("Formula")
library("Hmisc")
library(reshape2)

ln_LNCaP_RT_cons <- read.delim('mutations_conserved_LNCap_with_LNCap_RT.txt')
ln_LNCaP_RT_losses <- read.delim('mutations_losses_LNCap_with_LNCap_RT.txt')
ln_LNCaP_RT_gains <- read.delim("mutations_gains_LNCap_with_LNCaP_RT.txt")

ln_LNCaP_RT_cons$study <- "conserved"
ln_LNCaP_RT_losses$study <- "losses"
ln_LNCaP_RT_gains$study <- "gains"

ln_LNCaP_RT_cons$width <- with(ln_LNCaP_RT_cons, end - start)
ln_LNCaP_RT_losses$width <- with(ln_LNCaP_RT_losses, end - start)
ln_LNCaP_RT_gains$width <- with(ln_LNCaP_RT_gains, end - start)

df <- rbind(ln_LNCaP_RT_cons[,c("ID", "mean_rt", "study", "n","width")]
            , ln_LNCaP_RT_losses[,c("ID", "mean_rt", "study", "n","width")]
            , ln_LNCaP_RT_gains[,c("ID", "mean_rt", "study", "n","width")])
df$study <- factor(df$study, levels = c("gains", "conserved", "losses"))

df <- unique(df)

pdf('boxplot_ln_LNCap_rt_cons_gains_losses_mut_adjustedBy_width_log.e.pdf')
ggplot(df, aes(x = study, y = log(n/width), fill = study)) +
  geom_boxplot() + theme_classic() + scale_fill_brewer(palette="Set1")
dev.off()

df$log_n_by_width <- with(df, log(n/width))

x <- wilcox.test(df[df$study=='gains',"log_n_by_width"]
                 ,df[df$study=='conserved',"log_n_by_width"])
# Wilcoxon rank sum test with continuity correction
# 
# data:  df[df$study == "gains", "log_n_by_width"] and df[df$study == "conserved", "log_n_by_width"]
# W = 91143312, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

# x$p.value
# [1] 8.197485e-172

x <- wilcox.test(df[df$study=='losses',"log_n_by_width"]
                 ,df[df$study=='conserved',"log_n_by_width"])
# Wilcoxon rank sum test with continuity correction
# 
# data:  df[df$study == "losses", "log_n_by_width"] and df[df$study == "conserved", "log_n_by_width"]
# W = 25810822, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

#x$p.value
#6.637546e-147

table(df$study)
# gains conserved   losses 
# 24672      5996   6814 
