library("gkmSVM", lib = "/g/data/zk16/software/Rpackages_paola/R_4.0.0")
library(Biostrings)
library("BSgenome.Hsapiens.UCSC.hg19.masked")

bed_path <- "/human_data/rt/liver/feb2022/"

setwd("./replication_timing/human")

set.seed(1)
genNullSeqs(
  paste0(bed_path, "liver_non.recent_early_quintile_H9RT_ext.bed"),
  genomeVersion='hg19',
  outputBedFN = 'conserved_early_quintile_H9RT_ext_NULL.bed',
  outputPosFastaFN = 'posSet.fa',
  outputNegFastaFN = 'conserved_early_quintile_H9RT_ext_NULL.fa',
  xfold = 1,
  repeat_match_tol = 0.02,
  GC_match_tol = 0.02,
  length_match_tol = 0.02,
  batchsize = 5000,
  nMaxTrials = 50)

set.seed(1)
genNullSeqs(
  inputBedFN=  paste0(bed_path, "liver_non.recent_late_quintile_H9RT_ext.bed"),
  genomeVersion='hg19',
  outputBedFN = 'conserved_late_quintile_H9RT_ext_NULL.bed',
  outputPosFastaFN = 'posSet.fa',
  outputNegFastaFN = 'conserved_late_quintile_H9RT_ext_NULL.fa',
  xfold = 1,
  repeat_match_tol = 0.02,
  GC_match_tol = 0.02,
  length_match_tol = 0.02,
  batchsize = 5000,
  nMaxTrials = 50)

set.seed(1)
genNullSeqs(
  paste0(bed_path, "liver_recent_early_quintile_H9RT_ext.bed"),
  genomeVersion='hg19',
  outputBedFN = 'recent_early_quintile_H9RT_ext_NULL.bed',
  outputPosFastaFN = 'posSet.fa',
  outputNegFastaFN = 'recent_early_quintile_H9RT_ext_NULL.fa',
  xfold = 1,
  repeat_match_tol = 0.02,
  GC_match_tol = 0.02,
  length_match_tol = 0.02,
  batchsize = 5000,
  nMaxTrials = 50)

set.seed(1)
genNullSeqs(
  paste0(bed_path, "liver_recent_late_quintile_H9RT_ext.bed"),
  genomeVersion='hg19',
  outputBedFN = 'recent_late_quintile_H9RT_ext_NULL.bed',
  outputPosFastaFN = 'posSet.fa',
  outputNegFastaFN = 'recent_late_quintile_H9RT_ext_NULL.fa',
  xfold = 1,
  repeat_match_tol = 0.02,
  GC_match_tol = 0.02,
  length_match_tol = 0.02,
  batchsize = 5000,
  nMaxTrials = 50)
