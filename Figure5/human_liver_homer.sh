genome_path=/human_data/hg19.fa
bed_path=/rt/liver/feb2022/

annotatePeaks.pl \
${bed_path}liver_non.recent_early_quintile_H9RT_ext.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf hg19.ensGene.gtf -size given \
-mbed ./homer/liver_non.recent_early_quintile_H9RT_ext.humanJAS.bed  \
> ./homer/liver_non.recent_early_quintile_H9RT_ext.humanJAS.tsv

annotatePeaks.pl \
${bed_path}liver_non.recent_late_quintile_H9RT_ext.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf hg19.ensGene.gtf -size given \
-mbed ./homer/liver_non.recent_late_quintile_H9RT_ext.humanJAS.bed  \
> ./homer/liver_non.recent_late_quintile_H9RT_ext.humanJAS.tsv

annotatePeaks.pl \
${bed_path}liver_recent_early_quintile_H9RT_ext.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf hg19.ensGene.gtf -size given \
-mbed ./homer/liver_recent_early_quintile_H9RT_ext.humanJAS.bed  \
> ./homer/liver_recent_early_quintile_H9RT_ext.humanJAS.tsv

annotatePeaks.pl \
${bed_path}liver_recent_late_quintile_H9RT_ext.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf hg19.ensGene.gtf -size given \
-mbed ./homer/liver_recent_late_quintile_H9RT_ext.humanJAS.bed  \
> ./homer/liver_recent_late_quintile_H9RT_ext.humanJAS.tsv
