genome_path=/human_data/hg19.fa
bed_path=/replication_timing/human/

annotatePeaks.pl \
${bed_path}conserved_early_quintile_H9RT_ext_NULL.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf hg19.ensGene.gtf -size given \
-mbed ./homer/conserved_early_quintile_H9RT_ext_NULL_humanJAS.bed  \
> ./homer/conserved_early_quintile_H9RT_ext_NULL_humanJAS.tsv

annotatePeaks.pl \
${bed_path}conserved_late_quintile_H9RT_ext_NULL.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf hg19.ensGene.gtf -size given \
-mbed ./homer/conserved_late_quintile_H9RT_ext_NULL_humanJAS.bed  \
> ./homer/conserved_late_quintile_H9RT_ext_NULL_humanJAS.tsv

annotatePeaks.pl \
${bed_path}recent_early_quintile_H9RT_ext_NULL.bed \
${genome_path} \
-m /g/data/zk16/cc3704/jaspar_data/all_human_jaspar2020.motif -cpu 10 \
-gtf /g/data/zk16/cc3704/human_data/hg19.ensGene.gtf -size given \
-mbed ./homer/recent_early_quintile_H9RT_ext_NULL_humanJAS.bed  \
> ./homer/recent_early_quintile_H9RT_ext_NULL_humanJAS.tsv

annotatePeaks.pl \
${bed_path}recent_late_quintile_H9RT_ext_NULL.bed \
${genome_path} \
-m all_human_jaspar2020.motif -cpu 10 \
-gtf hg19.ensGene.gtf -size given \
-mbed ./homer/recent_late_quintile_H9RT_ext_NULL_humanJAS.bed  \
> ./homer/recent_late_quintile_H9RT_ext_NULL_humanJAS.tsv
