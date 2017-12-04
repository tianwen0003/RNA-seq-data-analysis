# path to reference and GTF file
GTF="/home/jiajinbu/nas3/lzz/Triticum_aestivum.TGACv1.36.gtf"
REF="/home/jiajinbu/nas3/lzz/Triticum_aestivum.TGACv1.dna.nonchromosomal.fa"

# path to raw row_data and clean_data
#rowData="/home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir"
#cleanData="/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cleanData_snakemake"

# path to java,trimmomatic and trimmomatic_adapter_file
java="/home/jiajinbu/bu/soft/sysoft/jre1.8.0_101/bin/java"
trimmomatic="/home/jiajinbu/bu/soft/bio/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar"
adapter="/home/jiajinbu/bu/soft/bio/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10"

# the classify of samples
CLASS1 = 'control'.split()
CLASS2 = "heat_1h".split()
CLASS3 = "drought_1h".split()
CLASS4 = "heat+drought_6h".split()
SAMPLE = CLASS1 + CLASS2 + CLASS3 + CLASS4
REP=[1,2]

# define the sample group will be used in cuffdiff connand
CLASS1_BAM = expand('/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.bam',sample=CLASS1,rep=[1,2])
CLASS2_BAM = expand('/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.bam',sample=CLASS2,rep=[1,2])
CLASS3_BAM = expand('/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.bam',sample=CLASS3,rep=[1,2])
CLASS4_BAM = expand('/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.bam',sample=CLASS4,rep=[1,2])

# define the lable will be used in cuffdiff command
LABELS=['control','heat_1h','drought_1h','heat+drought_6h']

rule all:
    input:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cuffdiff_snakemake/gene_exp.diff",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cuffdiff_snakemake/isoform_exp.diff",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/merged/merged.gtf"

# filter the raw data by trimmomatic
rule trimmomatic:
    input:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/{sample}_rep{rep}_1.fastq.gz",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/row_data/conver_dir/{sample}_rep{rep}_2.fastq.gz"
    output:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cleanData_snakemake/{sample}_rep{rep}_1_paired.fq.gz",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cleanData_snakemake/{sample}_rep{rep}_1_unpaired.fq.gz",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cleanData_snakemake/{sample}_rep{rep}_2_paired.fq.gz",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cleanData_snakemake/{sample}_rep{rep}_2_unpaired.fq.gz"
    threads: 8
    shell:
        "{java} -jar {trimmomatic} PE -threads {threads} {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:{adapter} LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# reads mapping by hisat2
rule hisat2_map:
    input:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cleanData_snakemake/{sample}_rep{rep}_1_paired.fq.gz",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cleanData_snakemake/{sample}_rep{rep}_2_paired.fq.gz"
    output:
        temp("/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.sam")
    threads: 12
    shell:
        "hisat2 -p {threads} -x {GTF} -1 {input[0]} -2 {input[1]} -S {output}"

# transform the sam file to sorted bam file
rule sam_to_sorted_bam:
    input:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.sam"
    output:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.bam"
    threads: 8
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

# estimate the abundance of gene expression and assemble the new transcripts
rule cufflinks:
    input:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/hisat_snakemake/{sample}_rep{rep}.bam"
    output:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/{sample}_rep{rep}/transcripts.gtf",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/{sample}_rep{rep}/genes.fpkm_tracking",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/{sample}_rep{rep}/isoforms.fpkm_tracking",
        dir="/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/{sample}_rep{rep}"
    threads: 12
    shell:
        "cufflinks -p {threads} -o {output.dir} -G {GTF} {input}"

# creat the gtf files list
rule compse_merge_file:
    input:
        expand("/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/{sample}_rep{rep}/transcripts.gtf", sample=SAMPLE, rep=[1,2])
    output:
        txt='/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/assemblies.txt'
    run:
        with open(output.txt,'w') as out:
                print(*input, sep="\n", file=out)

# generate the merged transcripts file
rule cuffmerge:
    input:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/assemblies.txt"
    output:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/merged/merged.gtf",
        dir="/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/merged"
    threads: 8
    shell:
        "cuffmerge -p {threads} -g {GTF} -s {REF} -o {output.dir} {input}"

# comparision
rule compare_assemblies:
    input:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/merged/merged.gtf"
    output:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/comparison/all.stats",
        dir="/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/comparison"
    shell:
        "cuffcompare -o {output.dir}all -s {REF} -r {GTF} {input}"

# find the differential expression gene
rule cuffdiff:
    input:
        class1=CLASS1_BAM,
        class2=CLASS2_BAM,
        class3=CLASS3_BAM,
        class4=CLASS4_BAM,
#        gtf="/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cufflinks_snakemake/merged/merged.gtf"
    output:
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cuffdiff_snakemake/gene_exp.diff",
        "/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cuffdiff_snakemake/isoform_exp.diff",
        dir="/home/jiajinbu/nas3/lzz/10_16_rna-seq/snakemake/cuffdiff_snakemake"
    params:
        class1=",".join(CLASS1_BAM),
        class2=",".join(CLASS2_BAM),
        class3=",".join(CLASS3_BAM),
        class4=",".join(CLASS4_BAM),
        labels=",".join(LABELS)
    threads: 8
    shell:
        "cuffdiff -p {threads} -o {output.dir} -L {params.labels} {GTF} {params.class1} {params.class2} {params.class3} {params.class4}"
