#Create landscape file for chromosome 1
java -jar /g/data/zk16/cc3704/tools/RepliconWrench.jar BedToIPLS \
-i ENCFF728ZKQ_testis_DNAseq_mm10_small_norm.bed \
-o testis_DNAseq_mm10_norm

#Run Replicon
cd /g/data/zk16/cc3704/replication_timing
export PATH=/g/data/zk16/xzhang/bin:$PATH
mouse_path="/g/data/zk16/cc3704/replication_timing/mouse/testis_DNAseq_mm10_norm/chr1/"
human_path="/g/data/zk16/cc3704/replication_timing/human/ENCFF289LWU_ATAC_testis_hg19_IPLS/chr1/"
replicon -threads 8 ${mouse_path}chr1_ENCFF728ZKQ_testis_DNAseq_mm10_small_norm chr1_mm10_norm

# CGHnTiming to .timing
rep_path="/g/data/zk16/cc3704/tools/"
java -jar ${rep_path}RepliconWrench.jar CghToTiming -i chr1_mm10_norm.CGHnTiming.csv > chr1_mm10_norm.timing
