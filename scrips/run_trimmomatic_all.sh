for sample in control drought_1h heat_1h heat+drought_6h
do
  for rep in 1 2
  do
    /home/jiajinbu/bu/soft/sysoft/jre1.8.0_101/bin/java -jar /home/jiajinbu/bu/soft/bio/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/${sample}_rep${rep}_1.fastq.gz /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/${sample}_rep${rep}_2.fastq.gz /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/clean_data/${sample}_rep${rep}_1_paired.fq.gz /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/clean_data/${sample}_rep${rep}_1_unpaired.fq.gz /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/clean_data/${sample}_rep${rep}_2_paired.fq.gz /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/clean_data/${sample}_rep${rep}_2_unpaired.fq.gz ILLUMINACLIP:/home/jiajinbu/bu/soft/bio/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  done
done
