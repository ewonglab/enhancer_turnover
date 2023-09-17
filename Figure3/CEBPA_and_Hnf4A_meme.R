# Script used to generate motifs in meme format with human and mouse CEBPA and Hnf4A motifs to use with FIMO

library("universalmotif")

# HUMAN MOTIFS

# Reading Cis-BP 2.0 human motifs
motif_file <- "Homo_sapiens.meme"
motifs <- read_meme(motif_file)

# Get motif IDs
motif_ids <- sapply(1:length(motifs),  function(x) motifs[[x]]@altname)

motifs[[grep("CEBPA", motif_ids, ignore.case = T)]]

# Motif name:   M04045_2.00
# Alternate name:   CEBPA
# Alphabet:   DNA
# Type:   PPM
# Strands:   +-
#   Total IC:   9.85
# Pseudocount:   1
# Consensus:   NNTKRCGCAAYN
# Target sites:   1
# E-value:   0
# Extra info:   [eval.string] 0
#
# N    N    T    K    R    C    G    C    A    A    Y    N
# A 0.17 0.37 0.01 0.00 0.46 0.04 0.19 0.13 0.89 0.98 0.03 0.23
# C 0.30 0.24 0.01 0.01 0.01 0.73 0.07 0.80 0.10 0.00 0.39 0.33
# G 0.32 0.37 0.01 0.26 0.43 0.06 0.70 0.01 0.00 0.01 0.19 0.26
# T 0.21 0.03 0.98 0.73 0.10 0.18 0.04 0.06 0.00 0.00 0.39 0.17


write_meme(motifs[[grep("CEBPA", motif_ids, ignore.case = T)]]
           , "human_Cebpa_CISBP_2.00.meme")


motifs[[grep("Hnf4a", motif_ids, ignore.case = T)]]


# Motif name:   M03357_2.00
# Alternate name:   HNF4A
# Alphabet:   DNA
# Type:   PPM
# Strands:   +-
#   Total IC:   19.82
# Pseudocount:   1
# Consensus:   RWTGGACTTTGGACYN
# Target sites:   1
# E-value:   0
# Extra info:   [eval.string] 0
#
# R    W    T    G    G A C    T    T    T    G    G    A    C    Y    N
# A 0.30 0.47 0.13 0.03 0.19 1 0 0.00 0.01 0.00 0.06 0.08 0.98 0.05 0.01 0.08
# C 0.00 0.05 0.10 0.00 0.09 0 1 0.04 0.01 0.03 0.00 0.19 0.01 0.94 0.43 0.41
# G 0.66 0.15 0.02 0.97 0.68 0 0 0.00 0.00 0.00 0.91 0.63 0.00 0.00 0.03 0.31
# T 0.04 0.33 0.75 0.00 0.04 0 0 0.96 0.98 0.97 0.03 0.10 0.00 0.01 0.53 0.20

write_meme(motifs[[grep("Hnf4a", motif_ids, ignore.case = T)]]
           , "human_Hnf4a_CISBP_2.00.meme")

                    
                    
# MOUSE MOTIFS

motif_file <- "Mus_musculus.meme"
motifs <- read_meme(motif_file)

                    # get motif ids
motif_ids <- sapply(1:length(motifs),  function(x) motifs[[x]]@altname)
motifs[[grep("CEBPA", motif_ids, ignore.case = T)]]

# Motif name:   M00993_2.00
# Alternate name:   Cebpa
# Alphabet:   DNA
# Type:   PPM
# Strands:   +-
#   Total IC:   12.47
# Pseudocount:   1
# Consensus:   NTTRCGCAA
# Target sites:   1
# E-value:   0
# Extra info:   [eval.string] 0
#
# N T    T    R C    G    C A    A
# A 0.27 0 0.00 0.30 0 0.03 0.06 1 0.98
# C 0.24 0 0.00 0.00 1 0.00 0.82 0 0.01
# G 0.24 0 0.07 0.51 0 0.97 0.00 0 0.01
# T 0.24 1 0.93 0.19 0 0.00 0.11 0 0.01

write_meme(motifs[[grep("CEBPA", motif_ids, ignore.case = T)]]
           , "mouse_Cebpa_CISBP_2.00.meme")
motifs[[grep("Hnf4a", motif_ids, ignore.case = T)]]

# Motif name:   M00182_2.00
# Alternate name:   Hnf4a
# Alphabet:   DNA
# Type:   PPM
# Strands:   +-
#   Total IC:   4.2
# Pseudocount:   1
# Consensus:   NNRGKYCRNN
# Target sites:   1
# E-value:   0
# Extra info:   [eval.string] 0
#
# N    N    R    G    K    Y    C    R    N    N
# A 0.34 0.31 0.38 0.02 0.07 0.13 0.15 0.59 0.31 0.27
# C 0.23 0.15 0.05 0.00 0.00 0.38 0.65 0.09 0.26 0.24
# G 0.22 0.41 0.52 0.94 0.42 0.10 0.10 0.25 0.20 0.23
# T 0.22 0.13 0.04 0.05 0.51 0.39 0.09 0.07 0.23 0.26

write_meme(motifs[[grep("Hnf4a", motif_ids, ignore.case = T)]]
           , "mouse_Hnf4a_CISBP_2.00.meme")
