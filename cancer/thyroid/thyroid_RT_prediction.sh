# liftOver ATAC-seq to hg19
cd /g/data/zk16/cc3704/replication_timing/human
/g/data/zk16/cc3704/tools/liftOver ENCFF682VUY_thyroid_ATAC_hg38_scaled.bed \
/g/data/zk16/cc3704/useful/chain_files/hg38ToHg19.over.chain.gz \
ENCFF682VUY_thyroid_ATAC_hg19_scaled.bed ENCFF682VUY_thyroid_ATAC_hg38_scaled_unmap.bed

#create landscape files
java -jar /g/data/zk16/cc3704/tools/RepliconWrench.jar BedToIPLS -i ENCFF682VUY_thyroid_ATAC_hg19_scaled.bed -o thyroid_hg19_IPLS

#run Replicon for every chromosome
cd /g/data/zk16/cc3704/replication_timing
export PATH=/g/data/zk16/xzhang/bin:$PATH

thyroid_path="/g/data/zk16/cc3704/replication_timing/human/thyroid_hg19_IPLS/"
thyroid_base="ENCFF682VUY_thyroid_ATAC_hg19_scaled"

cd ${thyroid_path}
# replicon -threads 8 ./chr1/chr1_${thyroid_base} chr1
# replicon -threads 8 ./chr15/chr15_${thyroid_base} chr15
# replicon -threads 8 ./chr4/chr4_${thyroid_base} chr4
# replicon -threads 8 ./chrX/chrX_${thyroid_base} chrX
# replicon -threads 8 ./chr10/chr10_${thyroid_base} chr10
# replicon -threads 8 ./chr16/chr16_${thyroid_base} chr16
# replicon -threads 8 ./chr2/chr2_${thyroid_base} chr2
# replicon -threads 8 ./chr9/chr9_${thyroid_base} chr9
# replicon -threads 8 ./chrY/chrY_${thyroid_base} chrY
# replicon -threads 8 ./chr11/chr11_${thyroid_base} chr11
# replicon -threads 8 ./chr17/chr17_${thyroid_base} chr17
# replicon -threads 8 ./chr20/chr20_${thyroid_base} chr20
# replicon -threads 8 ./chr5/chr5_${thyroid_base} chr5
# replicon -threads 8 ./chr12/chr12_${thyroid_base} chr12
# replicon -threads 8 ./chr21/chr21_${thyroid_base} chr21
# replicon -threads 8 ./chr6/chr6_${thyroid_base} chr6
# replicon -threads 8 ./chr13/chr13_${thyroid_base} chr13
# replicon -threads 8 ./chr18/chr18_${thyroid_base} chr18
# replicon -threads 8 ./chr22/chr22_${thyroid_base} chr22
# replicon -threads 8 ./chr7/chr7_${thyroid_base} chr7
# replicon -threads 8 ./chr14/chr14_${thyroid_base} chr14
# replicon -threads 8 ./chr19/chr19_${thyroid_base} chr19
# replicon -threads 8 ./chr3/chr3_${thyroid_base} chr3
# replicon -threads 8 ./chr8/chr8_${thyroid_base} chr8

# replicon -threads 8 ./chr17_gl000204_random/chr17_gl000204_random_${thyroid_base} chr17_gl000204_random
# replicon -threads 8 ./chrUn_gl000229/chrUn_gl000229_${thyroid_base} chrUn_gl000229
# replicon -threads 8 ./chrUn_gl000211/chrUn_gl000211_${thyroid_base} chrUn_gl000211
# replicon -threads 8 ./chrUn_gl000237/chrUn_gl000237_${thyroid_base} chrUn_gl000237
# replicon -threads 8 ./chrUn_gl000212/chrUn_gl000212_${thyroid_base} chrUn_gl000212
# replicon -threads 8 ./chr11_gl000202_random/chr11_gl000202_random_${thyroid_base} chr11_gl000202_random
# replicon -threads 8 ./chr1_gl000192_random/chr1_gl000192_random_${thyroid_base} chr1_gl000192_random
# replicon -threads 8 ./chrUn_gl000220/chrUn_gl000220_${thyroid_base} chrUn_gl000220
# replicon -threads 8 ./chr8_gl000196_random/chr8_gl000196_random_${thyroid_base} chr8_gl000196_random
# replicon -threads 8 ./chrUn_gl000228/chrUn_gl000228_${thyroid_base} chrUn_gl000228

#CGHnTiming to timing

rep_path="/g/data/zk16/cc3704/tools/"
cd /g/data/zk16/cc3704/replication_timing/human/thyroid_hg19_IPLS
for cgh in *.CGHnTiming.csv
do
new_name="${cgh/CGHnTiming.csv/}"
java -jar ${rep_path}RepliconWrench.jar CghToTiming -i ${cgh} > ${new_name}timing
done
