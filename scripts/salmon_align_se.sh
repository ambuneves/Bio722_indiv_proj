####Salmon Counts w/ decoy ####
#This is a for loop that Ian sent me:
## SALMON counts with index

counts_dir=/2/scratch/amandaN/cornell_seqdata/salmon_counts

trim_dir=/2/scratch/amandaN/cornell_seqdata/dgrp_seq/trim_seq/bbduk_trim

index_dir=/2/scratch/amandaN/salmon_index_trial/salmon_index

files=(${trim_dir}/*.fastq)


#Real deal for loop
for file in ${files[@]}
do
  name=${file}
  sample_name=`basename ${name} .fastq`
salmon quant -i ${index_dir} -l A \
  -r ${trim_dir}/${sample_name}.fastq \
  -p 16 --validateMappings --rangeFactorizationBins 4 \
  --seqBias --gcBias --recoverOrphans \
  -o ${counts_dir}/${sample_name}_quant
done #I modified this to be SE instead of PE data
