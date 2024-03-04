library(ggplot2)
library(RColorBrewer)
library(data.table)
library(tidyr)
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

df <- rbind(ln_LNCaP_RT_cons[,c("ID", "rt_quintile", "study", "n","width")]
            , ln_LNCaP_RT_losses[,c("ID", "rt_quintile", "study", "n","width")]
            , ln_LNCaP_RT_gains[,c("ID", "rt_quintile", "study", "n","width")])
df$study <- factor(df$study, levels = c("gains", "conserved","losses"))

df <- unique(df)
df$log_n_by_width <- with(df, log(n/width))
df$rt_quintile <- factor(df$rt_quintile, levels=c(1,2,3,4,5))

pdf('boxplot_ln_LNCap_RTquintiles_cons_gains_losses_mut_adjustedBy_width_log.e.pdf')
ggplot(df, aes(x = rt_quintile, y = log_n_by_width, fill = study)) +
  geom_boxplot() + theme_classic() + scale_fill_brewer(palette="Set1") +
  facet_wrap(~study)
dev.off()

pdf('boxplot_ln_LNCap_RTquintiles_cons_gains_losses_mut_adjustedBy_width_log.e.1.pdf')
ggplot(df, aes(x = study, y = log_n_by_width, fill = study)) +
  geom_boxplot() + theme_classic() + scale_fill_brewer(palette="Set1") +
  facet_wrap(~rt_quintile)
dev.off()

### testing the difference between unchanged enhancers and gains/losses
for(i in 1:5){
  print(paste("quintile", i))
  x <- wilcox.test(df[df$rt_quintile==i & df$study == "conserved", "log_n_by_width"]
                   , df[df$rt_quintile==i & df$study == "gains", "log_n_by_width"])
  print(paste("gains vs unchanged:", paste("p-value=", x$p.value)))
  x <- wilcox.test(df[df$rt_quintile==i & df$study == "conserved", "log_n_by_width"]
                   , df[df$rt_quintile==i & df$study == "losses", "log_n_by_width"])
  print(paste("losses vs unchanged:", paste("p-value=", x$p.value)))}


# [1] "quintile 1"
# [1] "gains vs unchanged: p-value= 2.45694850969587e-07"
# [1] "losses vs unchanged: p-value= 0.147820818660658"
# [1] "quintile 2"
# [1] "gains vs unchanged: p-value= 6.40061463570428e-07"
# [1] "losses vs unchanged: p-value= 0.000279819325145233"
# [1] "quintile 3"
# [1] "gains vs unchanged: p-value= 2.02563146737779e-21"
# [1] "losses vs unchanged: p-value= 1.57697550150131e-22"
# [1] "quintile 4"
# [1] "gains vs unchanged: p-value= 7.78613075691398e-41"
# [1] "losses vs unchanged: p-value= 4.0858824903647e-44"
# [1] "quintile 5"
# [1] "gains vs unchanged: p-value= 2.2717176704556e-61"
# [1] "losses vs unchanged: p-value= 1.01888533717475e-80"


### line plots with error bars
df_medians <- aggregate(log_n_by_width~(rt_quintile + study), df, median)
colnames(df_medians) <- c("rt_quintile", "study", "median")
df_medians$rt_quintile <- factor(df_medians$rt_quintile, levels = c(1,2,3,4,5))

se <- function(x) sqrt(var(x) / length(x))
df_se <- aggregate(log_n_by_width~(rt_quintile + study), df, se)
colnames(df_se) <- c("rt_quintile", "study", "se")
df_se$rt_quintile <- factor(df_se$rt_quintile, levels = c(1,2,3,4,5))

df_medians <- merge(df_medians, df_se, by=c("rt_quintile", "study"))
(df_medians)

#    rt_quintile     study    median          se
# 1            1 conserved -7.090077 0.076147531
# 2            1     gains -6.745236 0.015031948
# 3            1    losses -7.029620 0.043405182
# 4            2 conserved -7.244228 0.047353373
# 5            2     gains -6.956545 0.010254045
# 6            2    losses -7.090077 0.031431091
# 7            3 conserved -7.377759 0.030191630
# 8            3     gains -7.003065 0.009518009
# 9            3    losses -6.907755 0.027538504
# 10           4 conserved -7.495542 0.022575442
# 11           4     gains -7.047517 0.008814062
# 12           4    losses -6.972294 0.022872535
# 13           5 conserved -7.377759 0.015885296
# 14           5     gains -6.956545 0.008106921
# 15           5    losses -6.907755 0.015669477

pdf("prostate_enh_median_n_by_enhWidth_LNCaP_RT.1.pdf")
ggplot(df_medians, aes(x=rt_quintile, y=(median), group=study, color=study)) +
  geom_line() + theme_classic() +
  geom_errorbar(aes(ymin=(median-se), ymax=(median+se)), width=.2,
                position=position_dodge(0.05))
dev.off()

(table(df[, c("study", "rt_quintile")]))
# rt_quintile
# study          1    2    3    4    5
# gains     2657 3884 4784 5814 7547
# conserved  120  367  857 1517 3138
# losses     370  699 1046 1489 3212


