#This code is adapted from what Tyler sent me, plus I played around with the settings a bit

raw_dir=/2/scratch/amandaN/cornell_seqdata/dgrp_seq/txt_files/fastq
trim_dir=/2/scratch/amandaN/cornell_seqdata/dgrp_seq/trim_seq/bbduk_trim

files=(${raw_dir}/*.fastq)

#For loop over every file
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .fastq`
/usr/local/BBmap/bbduk.sh \
in=${raw_dir}/${base}.fastq \
out=${trim_dir}/${base}_trim-bbduk.fastq \
/2/scratch/amandaN/cornell_seqdata/BBMap_adapters.fa \
threads=8 ktrim=r k=23 mink=10 hdist=1 tpe tbo \
qtrim=rl trimq=15 minlength=36 2> /2/scratch/amandaN/cornell_seqdata/dgrp_seq/trim_seq/bbduk_trim/log//${base}.log

done
