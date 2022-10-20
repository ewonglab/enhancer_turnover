# liftOver ATAC-seq data for replication time prediction

cd /g/data/zk16/cc3704/replication_timing/human

#testis
/g/data/zk16/cc3704/tools/liftOver ENCFF289LWU_ATAC_testis_hg38_scaled.bed \
/g/data/zk16/cc3704/useful/chain_files/hg38ToHg19.over.chain.gz \
ENCFF289LWU_ATAC_testis_hg19_scaled.bed ENCFF289LWU_ATAC_testis_unmap_scaled.bed

#ovary
/g/data/zk16/cc3704/tools/liftOver ENCFF586QDS_ATAC_ovary_hg38_scaled.bed \
/g/data/zk16/cc3704/useful/chain_files/hg38ToHg19.over.chain.gz \
ENCFF586QDS_ATAC_ovary_hg19_scaled.bed ENCFF586QDS_ATAC_ovary_unmap_scaled.bed

#create landscape files
java -jar /g/data/zk16/cc3704/tools/RepliconWrench.jar BedToIPLS -i ENCFF289LWU_ATAC_testis_hg19_scaled.bed -o testis_hg19_IPLS
java -jar /g/data/zk16/cc3704/tools/RepliconWrench.jar BedToIPLS -i ENCFF586QDS_ATAC_ovary_hg19_scaled.bed -o ovary_hg19_IPLS

cd /g/data/zk16/cc3704/replication_timing
export PATH=/g/data/zk16/xzhang/bin:$PATH

ovary_path="/g/data/zk16/cc3704/replication_timing/human/ovary_hg19_IPLS/"
testis_path="/g/data/zk16/cc3704/replication_timing/human/testis_hg19_IPLS/"

ovary_base="ENCFF586QDS_ATAC_ovary_hg19_scaled"
testis_base="ENCFF289LWU_ATAC_testis_hg19_scaled"

cd ${ovary_path}
replicon -threads 8 ./chr1/chr1_${ovary_base} chr1
replicon -threads 8 ./chr15/chr15_${ovary_base} chr15
replicon -threads 8 ./chr4/chr4_${ovary_base} chr4
replicon -threads 8 ./chrX/chrX_${ovary_base} chrX
replicon -threads 8 ./chr10/chr10_${ovary_base} chr10
replicon -threads 8 ./chr16/chr16_${ovary_base} chr16
replicon -threads 8 ./chr2/chr2_${ovary_base} chr2
replicon -threads 8 ./chr9/chr9_${ovary_base} chr9
replicon -threads 8 ./chr11/chr11_${ovary_base} chr11
replicon -threads 8 ./chr17/chr17_${ovary_base} chr17
replicon -threads 8 ./chr20/chr20_${ovary_base} chr20
replicon -threads 8 ./chr5/chr5_${ovary_base} chr5
replicon -threads 8 ./chr12/chr12_${ovary_base} chr12
replicon -threads 8 ./chr21/chr21_${ovary_base} chr21
replicon -threads 8 ./chr6/chr6_${ovary_base} chr6
replicon -threads 8 ./chr13/chr13_${ovary_base} chr13
replicon -threads 8 ./chr18/chr18_${ovary_base} chr18
replicon -threads 8 ./chr22/chr22_${ovary_base} chr22
replicon -threads 8 ./chr7/chr7_${ovary_base} chr7
replicon -threads 8 ./chr14/chr14_${ovary_base} chr14
replicon -threads 8 ./chr19/chr19_${ovary_base} chr19
replicon -threads 8 ./chr3/chr3_${ovary_base} chr3
replicon -threads 8 ./chr8/chr8_${ovary_base} chr8

cd ${testis_path}
replicon ./chr1/chr1_${testis_base} chr1 #done
replicon -threads 8 ./chr15/chr15_${testis_base} chr15
replicon -threads 8 ./chr4/chr4_${testis_base} chr4
replicon -threads 8 ./chrX/chrX_${testis_base} chrX
replicon -threads 8 ./chr10/chr10_${testis_base} chr10
replicon -threads 8 ./chr16/chr16_${testis_base} chr16
replicon -threads 8 ./chr2/chr2_${testis_base} chr2
replicon -threads 8 ./chr9/chr9_${testis_base} chr9
replicon -threads 8 ./chrY/chrY_${testis_base} chrY
replicon -threads 8 ./chr11/chr11_${testis_base} chr11
replicon -threads 8 ./chr17/chr17_${testis_base} chr17
replicon -threads 8 ./chr20/chr20_${testis_base} chr20
replicon -threads 8 ./chr5/chr5_${testis_base} chr5
replicon -threads 8 ./chr12/chr12_${testis_base} chr12
replicon -threads 8 ./chr21/chr21_${testis_base} chr21
replicon -threads 8 ./chr6/chr6_${testis_base} chr6
replicon -threads 8 ./chr13/chr13_${testis_base} chr13
replicon -threads 8 ./chr18/chr18_${testis_base} chr18
replicon -threads 8 ./chr22/chr22_${testis_base} chr22
replicon -threads 8 ./chr7/chr7_${testis_base} chr7
replicon -threads 8 ./chr14/chr14_${testis_base} chr14
replicon -threads 8 ./chr19/chr19_${testis_base} chr19
replicon -threads 8 ./chr3/chr3_${testis_base} chr3
replicon -threads 8 ./chr8/chr8_${testis_base} chr8


# CGHnTiming files to .timing files
rep_path="/g/data/zk16/cc3704/tools/"
cd /g/data/zk16/cc3704/replication_timing/human/testis_hg19_IPLS
for cgh in *.CGHnTiming.csv
do
new_name="${cgh/CGHnTiming.csv/}"
java -jar ${rep_path}RepliconWrench.jar CghToTiming -i ${cgh} > ${new_name}timing
done

cd /g/data/zk16/cc3704/replication_timing/human/ovary_hg19_IPLS
for cgh in *.CGHnTiming.csv
do
new_name="${cgh/CGHnTiming.csv/}"
java -jar ${rep_path}RepliconWrench.jar CghToTiming -i ${cgh} > ${new_name}timing
done
