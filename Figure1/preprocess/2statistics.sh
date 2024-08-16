REF_DIR="/data/nier/ref/mm9"
#calculate Reads depth 
bedtools bamtobed -i $PREFIX.sort.rmdup.bam |bedtools intersect -a $REF_DIR/mm9.10k.win -b stdin -c > $PREFIX.mm9.wc

#Because the length of the replication timing domain(RTD) is ~50bp, only reserve the window, which has above 0.8 overlapped the replication timing domain.
#Besides, the unoverlaped windows with RTD are still reserved. 
#The file 'RT_46C_ESC_Int62150809_mm9.bedgraph' is download from https://www2.replicationdomain.com/. 
bedtools intersect -a $PREFIX.mm9.wc -b RT_46C_ESC_Int62150809_mm9.bedgraph -wo |awk '{if(($9/($7-$6))> 0.8)print}' |  sort -k1,1 -k2,2n -k3,3n -s -V > $PREFIX.rt.txt
#Remove duplicates windows
cat $PREFIX.rt.txt|awk '{FS=OFS="\t"}''{sum[$1"_"$2"_"$3"_"$4]+=$8;a[$1"_"$2"_"$3"_"$4]++}END{for(c in sum){printf("%s\t%f\n", c,sum[c]/a[c])}}' | tr "_" "\t" | sort -k1,1 -k2,2n -k3,3n -s -V  >$PREFIX.rt.uniq.txt
#obtain overlapped windows
bedtools intersect -a $PREFIX.mm9.wc -b $PREFIX.rt.uniq.txt -v -f 0.95 | awk '{print $0"\t0"}'> $PREFIX.rt.uniq.nonRT.txt
#merge all windows
cat $PREFIX.rt.uniq.txt $PREFIX.rt.uniq.nonRT.txt | sort -k1,1 -k2,2n -k3,3n -s -V > $PREFIX.rt.uniq.all.txt
