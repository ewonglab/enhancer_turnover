# Generate IPSL files
cd /g/data/zk16/cc3704/replication_timing/human
java -jar /g/data/zk16/cc3704/tools/RepliconWrench.jar BedToIPLS -i ./AML/AML_pHSC_mean_signal.scaled.bed -o AML_pHSC_hg19_IPLS

# Run Replicon
cd /g/data/zk16/cc3704/replication_timing
export PATH=/g/data/zk16/xzhang/bin:$PATH

aml_path="/g/data/zk16/cc3704/replication_timing/human/AML_pHSC_hg19_IPLS/"
aml_base="AML_pHSC_mean_signal"

replicon -threads 8 ./chr1/chr1_${aml_base} chr1
replicon -threads 8 ./chr15/chr15_${aml_base} chr15
replicon -threads 8 ./chr4/chr4_${aml_base} chr4
replicon -threads 8 ./chrX/chrX_${aml_base} chrX
replicon -threads 8 ./chr10/chr10_${aml_base} chr10
replicon -threads 8 ./chr16/chr16_${aml_base} chr16
replicon -threads 8 ./chr2/chr2_${aml_base} chr2
replicon -threads 8 ./chr9/chr9_${aml_base} chr9
replicon -threads 8 ./chrY/chrY_${aml_base} chrY
replicon -threads 8 ./chr11/chr11_${aml_base} chr11
replicon -threads 8 ./chr17/chr17_${aml_base} chr17
replicon -threads 8 ./chr20/chr20_${aml_base} chr20
replicon -threads 8 ./chr5/chr5_${aml_base} chr5
replicon -threads 8 ./chr12/chr12_${aml_base} chr12
replicon -threads 8 ./chr21/chr21_${aml_base} chr21
replicon -threads 8 ./chr6/chr6_${aml_base} chr6
replicon -threads 8 ./chr13/chr13_${aml_base} chr13
replicon -threads 8 ./chr18/chr18_${aml_base} chr18
replicon -threads 8 ./chr22/chr22_${aml_base} chr22
replicon -threads 8 ./chr7/chr7_${aml_base} chr7
replicon -threads 8 ./chr14/chr14_${aml_base} chr14
replicon -threads 8 ./chr19/chr19_${aml_base} chr19
replicon -threads 8 ./chr3/chr3_${aml_base} chr3
replicon -threads 8 ./chr8/chr8_${aml_base} chr8

# CGHnTiming.csv files to .timing files
rep_path="/g/data/zk16/cc3704/tools/"
cd /g/data/zk16/cc3704/replication_timing/human/AML_pHSC_hg19_IPLS
for cgh in *.CGHnTiming.csv
do
new_name="${cgh/CGHnTiming.csv/}"
java -jar ${rep_path}RepliconWrench.jar CghToTiming -i ${cgh} > ${new_name}timing
done
