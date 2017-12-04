for sample in control drought_1h heat_1h heat+drought_6h
do
	fastqc -o /home/jiajinbu/nas3/lzz/10_16_rna-seq/fastqc/fastqc_result -f fastq -t 12 --nogroup /home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/${sample}*.gz
done
