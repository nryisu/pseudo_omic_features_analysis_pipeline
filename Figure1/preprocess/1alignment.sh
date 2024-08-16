#mapping
REF="/data/nier/ref/mm9/mm9.fa"
bwa mem -t 8 -R "@RG\tID:esc\tLB:esc\tSM:esc\tPL:ILLUMINA" \
        $REF \
	DRR029032_1.fastq.gz DRR029032_2.fastq.gz| \
        samtools view -bt $REF.fai -o esc.bam -
samtools sort esc.bam -o esc.sort.bam
#remove duplicates
picard="/data/nier/bin/picard.jar"
PREFIX="esc"
samtools index $PREFIX.sort.bam
java -Xmx20g -jar $picard MarkDuplicates \
        INPUT=$PREFIX.sort.bam \
        OUTPUT=$PREFIX.sort.rmdup.bam \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true \
        METRICS_FILE=$PREFIX.metrics
java -Xmx20g -jar $picard BuildBamIndex \
        INPUT=$PREFIX.sort.rmdup.bam 

