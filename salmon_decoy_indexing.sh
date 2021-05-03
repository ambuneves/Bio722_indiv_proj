#Download the data

wget ftp://ftp.flybase.net/releases/FB2021_01/dmel_r6.38/fasta/dmel-all-chromosome-r6.38.fasta.gz

wget ftp://ftp.flybase.net/releases/FB2021_01/dmel_r6.38/fasta/dmel-all-transcript-r6.38.fasta.gz

wget ftp://ftp.flybase.net/releases/FB2021_01/dmel_r6.38/gtf/dmel-all-r6.38.gtf.gz

wget ftp://ftp.flybase.net/releases/current/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2021_01.tsv.gz

#Unzip the files
gunzip *.gz

#Create a text file containining the decoys
grep "^>" dmel-all-chromosome-r6.38.fasta | cut -d " " -f 1 > decoys.txt

#Put the chromosome and trasncriptome files into one new file called gentrome.fa
cat dmel-all-transcript-r6.38.fasta dmel-all-chromosome-r6.38.fasta  > gentrome.fa

#Check that the newly created file makes sense by comparing word counts
grep -c "^>" dmel-all-transcript-r6.38.fasta
grep -c "^>" dmel-all-chromosome-r6.38.fasta
grep -c "^>" gentrome.fa

#Run the index
salmon index -t gentrome.fa -d decoys.txt -p 12 -i salmon_index --gencode

#Check that the new salmon index with decoys works
index_dir=/2/scratch/amandaN/salmon_index_trial/salmon_index/
sample_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads
sample_name=samw_wings_starved_R3_GCCAAT_L004
out_dir=/2/scratch/amandaN/salmon_index_trial/salmon_counts_w_decoys
salmon quant -i ${index_dir} -l A \
  -1 ${sample_dir}/${sample_name}_R1_PE.fastq \
  -2 ${sample_dir}/${sample_name}_R2_PE.fastq \
  -p 16 --validateMappings --rangeFactorizationBins 4 \
  --seqBias --gcBias --recoverOrphans \
  -o ${out_dir}/${sample_name}_quant

#Ian compared the counts he obtained with the decoy method to counts from normal indexing, nd it seems like a gene eh knows to be highly expressed isn't showing up..... Why's that? Genome version info?
