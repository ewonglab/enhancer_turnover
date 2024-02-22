library(data.table)
library("Hmisc")
library("GenomicFeatures")
library(dplyr)

setwd("./replication_timing/human")
GC_all <- read.table(file = "prom_enh_and_random_gc_byRTquintile_sept2022", header = T, sep = "\t")

# same results when quintiles is a factor
kw_enh <- kruskal.test(pct_gc ~ quintile, data = GC_all[GC_all$type == "enhancers", ])
kw_enh

# Kruskal-Wallis rank sum test
#
# data:  pct_gc by quintile
# Kruskal-Wallis chi-squared = 2708.6, df = 4, p-value < 2.2e-16

kw_prom <- kruskal.test(pct_gc ~ quintile, data = GC_all[GC_all$type == "promoters", ])
kw_prom
# Kruskal-Wallis rank sum test
#
# data:  pct_gc by quintile
# Kruskal-Wallis chi-squared = 688.25, df = 4, p-value < 2.2e-16

kw_random <- kruskal.test(pct_gc ~ quintile, data = GC_all[GC_all$type == "random", ])
kw_random
# Kruskal-Wallis rank sum test
#
# data:  pct_gc by quintile
# Kruskal-Wallis chi-squared = 902.98, df = 4, p-value < 2.2e-16

kw_prom$p.value
# [1] 1.219427e-147
kw_random$p.value
# [1] 3.764577e-194
kw_enh$p.value
# [1] 0

