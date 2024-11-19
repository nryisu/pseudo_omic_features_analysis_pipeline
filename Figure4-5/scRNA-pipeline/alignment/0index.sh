#software
STAR='/data/biosoft/software/STAR-2.7.1a/bin/Linux_x86_64/STAR'
#file
REF='/data/reference/mouse/mm10_STAR/RNA/'
REF_file='/data/reference/mouse/mm10_STAR/RNA/mm10.fa'
GTF='/data/reference/mouse/mm10_STAR/RNA/Mus_musculus.GRCm38.93.gtf'
#build the reference index in folder of 'mm10',eg mm10:
$STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $REF --genomeFastaFiles $REF_file --sjdbGTFfile $GTF --sjdbOverhang 100
