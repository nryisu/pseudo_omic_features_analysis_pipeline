GTF='/data/reference/mouse/mm10_STAR/RNA/Mus_musculus.GRCm38.93.gtf'
cd bam
dirlist=$(ls -t *.bam | tr '\n' ' ')
echo $dirlist
featureCounts -a $GTF \
                -o mesc_npc.txt -g 'gene_name' $dirlist