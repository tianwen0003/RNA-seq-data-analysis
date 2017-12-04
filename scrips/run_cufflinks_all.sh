for sample in control heat_1h drought_1h heat+drought_6h
do
   for rep in 1 2
   do
      cufflinks -p 12 -o /home/jiajinbu/nas3/lzz/10_16_rna-seq/cufflinks_result/${sample}_rep${rep}/ -G /home/jiajinbu/nas3/lzz/Triticum_aestivum.TGACv1.36.gtf /home/jiajinbu/nas3/lzz/10_16_rna-seq/hisat2_result/${sample}_rep${rep}.bam
   done
done
