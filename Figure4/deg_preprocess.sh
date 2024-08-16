####################################################################################
###bulk RNA 
for i in *.sra; do fastq-dump --gzip --split-3 $i; done
#trim adapters and remove low quality reads
trim_galore  SRR5420975_1.fastq.gz SRR5420975_2.fastq.gz -o ./ -q 20 --length 20 --paired --gzip
#STAR mapping
star='/data/biosoft/software/STAR-2.7.1a/bin/Linux_x86_64/STAR'
$star --genomeDir /data/reference/human/hg38_STAR \
	--readFilesCommand zcat \
     --readFilesIn SRR5420975_1_val_1.fq.gz SRR5420975_2_val_2.fq.gz\
     --runThreadN 4 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts 

#calculate the reads count of gene
dirlist=$(ls -t ./*.bam | tr '\n' ' ')
echo $dirlist
featureCounts -a /data/nier/ref/hg38/annotation/Homo_sapiens.GRCh38.94.gtf \
                -o final_counts.txt -g 'gene_name' $dirlist

#bam to bw file
bamCoverage -bs 10 -b esc.bam -o esc.bw
####################################################################################

###single cell RNA 
#Please refer to the "scRNA_DEGs_pipeline.md" documentation.

