liver_active <- read.table(file = "mouse_liver_activeEnh_conservationInOthers"
                           , header = T, stringsAsFactors = F)
liver_primed <- read.table(file = "mouse_liver_primedEnh_conservationInOthers"
                           , header = T, stringsAsFactors = F)
testis_active <- read.table(file = "mouse_testis_activeEnh_conservationInOthers"
                            , header = T, stringsAsFactors = F)
testis_primed <- read.table(file = "mouse_testis_primedEnh_conservationInOthers"
                            , header = T, stringsAsFactors = F)
muscle_active <- read.table(file = "mouse_muscle_activeEnh_conservationInOthers"
                            , header = T, stringsAsFactors = F)
muscle_primed <- read.table(file = "mouse_muscle_primedEnh_conservationInOthers"
                            , header = T, stringsAsFactors = F)
brain_active <- read.table(file = "mouse_brain_activeEnh_conservationInOthers"
                           , header = T, stringsAsFactors = F)
brain_primed <- read.table(file = "mouse_brain_primedEnh_conservationInOthers"
                           , header = T, stringsAsFactors = F)

# add type and tissue
liver_active$type <- testis_active$type <-  
  muscle_active$type <- brain_active$type <- "active"

liver_primed$type <- testis_primed$type <-
  muscle_primed$type <- brain_primed$type  <- "poised"

liver_active$tissue <- liver_primed$tissue <- "liver"
testis_active$tissue <- testis_primed$tissue <- "testis"
muscle_active$tissue <- muscle_primed$tissue <- "muscle"
brain_active$tissue <- brain_primed$tissue <- "brain"

all_enh <- rbind(liver_active, liver_primed, testis_active
                 , testis_primed, muscle_active, muscle_primed
                 , brain_active, brain_primed)

# define conservation in each type and tissue
library(stringr)

conserved_n <- function(enh, n){
  enh$n_species <- str_count(enh$species, ";") + 1
  enh <- enh[enh$n_species >= n, ]
  return(enh)
}

# conservation in at least other 2 sp
liver_active_atLeast2sp <- conserved_n(liver_active, 2) 
testis_active_atLeast2sp <- conserved_n(testis_active, 2) 
muscle_active_atLeast2sp <- conserved_n(muscle_active, 2) 
brain_active_atLeast2sp <- conserved_n(brain_active, 2) 


liver_primed_atLeast2sp <- conserved_n(liver_primed, 2)
testis_primed_atLeast2sp <- conserved_n(testis_primed, 2)
muscle_primed_atLeast2sp <- conserved_n(muscle_primed, 2)
brain_primed_atLeast2sp <- conserved_n(brain_primed, 2)

active_atLeast2sp <- rbind(liver_active_atLeast2sp
                           , testis_active_atLeast2sp
                           , muscle_active_atLeast2sp
                           , brain_active_atLeast2sp)
poised_atLeast2sp <- rbind(liver_primed_atLeast2sp
                           , testis_primed_atLeast2sp
                           , muscle_primed_atLeast2sp
                           , brain_primed_atLeast2sp)

# collapse tissue
active_atLeast2sp <- aggregate(tissue ~., active_atLeast2sp, paste, collapse=",")
poised_atLeast2sp <- aggregate(tissue ~., poised_atLeast2sp, paste, collapse=",")

write.table(x = active_atLeast2sp, file = "active_conserved_atLeast2sp.txt"
            , row.names = F, sep = "\t", quote = F)
write.table(x = poised_atLeast2sp, file = "poised_conserved_atLeast2sp.txt"
            , row.names = F, sep = "\t", quote = F)


## conservation in at least other 4 sp
liver_active_atLeast4sp <- conserved_n(liver_active, 4)
testis_active_atLeast4sp <- conserved_n(testis_active, 4)
muscle_active_atLeast4sp <- conserved_n(muscle_active, 4)
brain_active_atLeast4sp <- conserved_n(brain_active, 4)

liver_primed_atLeast4sp <- conserved_n(liver_primed, 4)
testis_primed_atLeast4sp <- conserved_n(testis_primed, 4)
muscle_primed_atLeast4sp <- conserved_n(muscle_primed, 4)
brain_primed_atLeast4sp <- conserved_n(brain_primed, 4)

active_atLeast4sp <- rbind(liver_active_atLeast4sp
                           , testis_active_atLeast4sp
                           , muscle_active_atLeast4sp
                           , brain_active_atLeast4sp)
poised_atLeast4sp <- rbind(liver_primed_atLeast4sp
                           , testis_primed_atLeast4sp
                           , muscle_primed_atLeast4sp
                           , brain_primed_atLeast4sp)

active_atLeast4sp <- aggregate(tissue ~., active_atLeast4sp, paste, collapse=",")
poised_atLeast4sp <- aggregate(tissue ~., poised_atLeast4sp, paste, collapse=",")

write.table(x = active_atLeast4sp, file = "active_conserved_atLeast4sp.txt"
            , row.names = F, sep = "\t", quote = F)
write.table(x = poised_atLeast4sp, file = "poised_conserved_atLeast4sp.txt"
            , row.names = F, sep = "\t", quote = F)
