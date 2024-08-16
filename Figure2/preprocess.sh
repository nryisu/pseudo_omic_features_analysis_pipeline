#sra to fastq file format
for i in *.sra; do fastq-dump --gzip --split-3 $i; done
#filter low quality reads and remove adapters
trim_galore --clip_R1 4 --clip_R2 4 --three_prime_clip_R1 4 --three_prime_clip_R2 4 --paired --gzip SRR10325453_1.fastq.gz SRR10325453_2.fastq.gz
#mapping 
sample=$1
qc_R1=${sample}_1_val_1.fq.gz
qc_R2=${sample}_2_val_2.fq.gz

output=$2
outdir_map=/data/nier/scATAC/HCT116/bulkATAC/rawdata
ref="/data/nier/ref/hg19/hg19"
ref_size="/data/nier/ref/hg19/hg19.size"
picard="/data/nier/bin/picard.jar"

bowtie2 -p 6 -X 1000 --dovetail --no-unal --no-mixed --no-discordant --very-sensitive -I 0 -x $ref -1 $qc_R1 -2 $qc_R2  -S  $outdir_map/$output.sam
awk -F"\t" '{if($3 !~ /chrM/ && $3 !~ /chrY/) print}' $outdir_map/$output.sam |samtools view -S -b -f 0x2 -q 30 - |samtools sort - -o $outdir_map/$output.pe.q30.sort.bam
samtools view -Sb $outdir_map/$output.sam  > $outdir_map/$output.bam
#java -Xmx20G -jar $picard MarkDuplicates INPUT=$outdir_map/$output.pe.q30.sort.bam OUTPUT=$outdir_map/$output.pe.q30.sort.rmdup.bam METRICS_FILE=$outdir_map/$output.Picard_Metrics_unfiltered_bam.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
samtools merge G1.pe.q30.sort.bam G1_rep1.pe.q30.sort.bam G1_rep2.pe.q30.sort.bam
#remover duplicates
ref="/data/nier/ref/hg19/hg19"
ref_size="/data/nier/ref/hg19/hg19.size"
picard="/data/nier/bin/picard.jar"
output="G1"
java -Xmx20G -jar $picard MarkDuplicates INPUT=$output.pe.q30.sort.bam OUTPUT=$output.pe.q30.sort.rmdup.bam METRICS_FILE=$output.Picard_Metrics_unfiltered_bam.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true

#simulation:G1：S：G2=5:3:2. The G1 reads is 12035894,S reads is 30887308,G2 reads is 32986784. The total reads is equal to G1 reads.
samtools view -bs 10000.2338027 S.pe.q30.sort.rmdup.bam > S.subsample.bam 
#10000 seed
samtools view -bs 10000.1459481 G2M.pe.q30.sort.rmdup.bam > G2M.subsample.bam
samtools merge mixture.simulate.bam G1.pe.q30.sort.rmdup.bam S.subsample.bam G2M.subsample.bam
#call peak
python2 /data/nier/bin/macs2 callpeak -n S --format BAMPE -g mm -q 0.01 --keep-dup all \
					--shift 0 -f BAM -B --nomodel \
					-t S.subsample.bam \
					--outdir ./S_peak

#filter peak
bedtools intersect -v -a G1_peaks.narrowPeak -b /data/nier/ref/hg19/ref/hg19_blacklist.JDB.bed | grep -P 'chr[0-9XY]+(?!_)' | cut -f1-3 > G1_filter.bed
