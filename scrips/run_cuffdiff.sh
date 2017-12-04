cuffdiff -p 12 -o /home/jiajinbu/nas3/lzz/10_16_rna-seq/cuffdiff_result \
-L control,heat_1h,drought_1h,heat+drought_6h \
/home/jiajinbu/nas3/lzz/Triticum_aestivum.TGACv1.36.gtf \
/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/control_rep1_sorted.bam,/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/control_rep2.bam \
/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/heat_1h_rep1.bam,/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/heat_1h_rep2.bam \
/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/drought_1h_rep1.bam,/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/drought_1h_rep2.bam \
/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/heat+drought_6h_rep1.bam,/home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/heat+drought_6h_rep2.bam
