cd ./replication_timing/mouse/
 
genome_path=Mus_musculus.GRCm38.dna_sm.primary_assembly.fa

annotatePeaks.pl \ 
mouse_enh_meanRT_by_germ_line_cons2sp_EARLY.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf Mus_musculus.GRCm38.92.gtf -size given \
-mbed mouse_recent_and_cons2sp_EARLY_human.jaspar.bed  \
> mouse_recent_and_cons2sp_EARLY_human.jaspar.tsv

/g/data/zk16/software/homer/bin/annotatePeaks.pl \
mouse_enh_meanRT_by_germ_line_cons2sp_LATE.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf Mus_musculus.GRCm38.92.gtf -size given \
-mbed mouse_recent_and_cons2sp_LATE_human.jaspar.bed  \
> mouse_recent_and_cons2sp_LATE_human.jaspar.tsv
