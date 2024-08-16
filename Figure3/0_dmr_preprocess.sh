###dowload data
#paired-end reads GSM1296803: Homo sapiens; Bisulfite-Seq
fastq-dump --split-3 --gzip SRR1057536
fastq-dump --split-3 --gzip SRR1057537
fastq-dump --split-3 --gzip SRR1057558
fastq-dump --split-3 --gzip SRR1057559

###trim adapter and low quliaty data
#WI38
cutadapt -u 18 -o SRR1057537_1.out.fastq.gz SRR1057537_1.fastq.gz
cutadapt -u 18 -o SRR1057537_2.out.fastq.gz SRR1057537_2.fastq.gz
trim_galore --paired --gzip SRR1057537_1.out.fastq.gz SRR1057537_2.out.fastq.gz
echo "trim WI38 rep1 done!"
cutadapt -u 18 -o SRR1057536_1.out.fastq.gz SRR1057536_1.fastq.gz
cutadapt -u 18 -o SRR1057536_2.out.fastq.gz SRR1057536_2.fastq.gz
trim_galore --paired --gzip SRR1057536_1.out.fastq.gz SRR1057536_2.out.fastq.gz
echo "trim WI38 rep2 done!"
#H9
cutadapt -u 18 -o SRR1057558_1.out.fastq.gz SRR1057558_1.fastq.gz
cutadapt -u 18 -o SRR1057558_2.out.fastq.gz SRR1057558_2.fastq.gz
trim_galore --paired --gzip SRR1057558_1.out.fastq.gz SRR1057558_2.out.fastq.gz
echo "trim H9 rep1 done!"
cutadapt -u 18 -o SRR1057559_1.out.fastq.gz SRR1057559_1.fastq.gz
cutadapt -u 18 -o SRR1057559_2.out.fastq.gz SRR1057559_2.fastq.gz
trim_galore --paired --gzip SRR1057559_1.out.fastq.gz SRR1057559_2.out.fastq.gz
 echo "trim H9 rep2 done!"

###mapping reads to reference, then sort
#mapping
/data/biosoft/software/bismark_v0.22.1/bismark --genome_folder /data/reference/human/hg38_bismark --bowtie2 --path_to_bowtie2 /data/biosoft/software/bowtie2-2.3.5.1-linux-x86_64 --bam --non_directional -N 1 -L 20 -1 SRR1057536_1.out_val_1.fq.gz -2 SRR1057536_2.out_val_2.fq.gz -o ./
/data/biosoft/software/bismark_v0.22.1/bismark --genome_folder /data/reference/human/hg38_bismark --bowtie2 --path_to_bowtie2 /data/biosoft/software/bowtie2-2.3.5.1-linux-x86_64 --bam --non_directional -N 1 -L 20 -1 SRR1057537_1.out_val_1.fq.gz -2 SRR1057537_2.out_val_2.fq.gz -o ./
/data/biosoft/software/bismark_v0.22.1/bismark --genome_folder /data/reference/human/hg38_bismark --bowtie2 --path_to_bowtie2 /data/biosoft/software/bowtie2-2.3.5.1-linux-x86_64 --bam --non_directional -N 1 -L 20 -1 SRR1057559_1.out_val_1.fq.gz -2 SRR1057559_2.out_val_2.fq.gz -o ./
/data/biosoft/software/bismark_v0.22.1/bismark --genome_folder /data/reference/human/hg38_bismark --bowtie2 --path_to_bowtie2 /data/biosoft/software/bowtie2-2.3.5.1-linux-x86_64 --bam --non_directional -N 1 -L 20 -1 SRR1057558_1.out_val_1.fq.gz -2 SRR1057558_2.out_val_2.fq.gz -o ./
#sort
samtools sort -T ./ -o h9_p2.bam SRR1057559_1.out_val_1_bismark_bt2_pe.bam
samtools sort -T ./ -o h9_p1.bam SRR1057558_1.out_val_1_bismark_bt2_pe.bam
samtools sort -T ./ -o wi38_p1.bam SRR1057537_1.out_val_1_bismark_bt2_pe.bam
samtools sort -T ./ -o wi38_p2.bam SRR1057536_1.out_val_1_bismark_bt2_pe.bam




