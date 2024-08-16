#build the reference index in folder of 'mm9' 
bwa index -a bwtsw mm9.fa
java -jar picard.jar CreateSequenceDictionary R=mm9.fa O=mm9.dict
samtools faidx mm9.fa
#build 10k window with 1bp step in folder of 'mm9' 
bedtools makewindows -g mm9.sizes -w 10000 -s 1 > mm9.10k.win