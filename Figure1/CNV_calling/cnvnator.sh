#cnv install in conda
conda create -n cnvnator python=3.8
conda activate cnvnator
conda install -c conda-forge root_base=6.20  
conda install -c bioconda cnvnator

#cnvnator
PREFIX="esc"
for chr in $(seq -f 'chr%g' 1 19)
do
#EXTRACTING READ MAPPING FROM BAM/SAM FILES
cnvnator -root ${PREFIX}.${chr}.root -genome mm9 -tree ${PREFIX}.${chr}.markdup.bam -chrom ${chr} -unique > ${PREFIX}.${chr}.extract.out
#GENERATING A HISTOGRAM
cnvnator -root ${PREFIX}.${chr}.root -his 200 -chrom ${chr} -d $REF > ${PREFIX}.${chr}.histogram.out
#CALCULATING STATISTICS
cnvnator -root ${PREFIX}.${chr}.root -chrom ${chr} -stat 200 > ${PREFIX}.${chr}.stats.out
#RD SIGNAL PARTITIONING
cnvnator -root ${PREFIX}.${chr}.root -chrom ${chr} -partition 200 > ${PREFIX}.${chr}.partition.out
#CNV CALLING
cnvnator -root ${PREFIX}.${chr}.root -chrom ${chr} -call 200 > ${PREFIX}.${chr}.cnv
done
echo "cnvnator done!"
#merge and filter
cat ${PREFIX}.${chr}.cnv >> {PREFIX}.cnv
awk -F"\t" '{split($2,a,":");split(a[2],b,"-");if($9<0.5)print a[1]"\t"b[1]"\t"b[2]"\t"$4"\t"$1}' {PREFIX}.cnv > {PREFIX}.bed