for sample in control drought_1h heat_1h heat+drought_6h
do
  for rep in 1 2
  do
    hisat2 -p 12 -x /home/jiajinbu/nas3/lzz/hisat2_index_Ta/hisat2_index_Ta -1 /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/clean_data/${sample}_rep${rep}_1_paired.fq.gz -2 /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/clean_data/${sample}_rep${rep}_2_paired.fq.gz -S /home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/${sample}_rep${rep}.sam
  done
done
