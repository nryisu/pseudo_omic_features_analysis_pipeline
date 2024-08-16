###filter CpG site: 1000bp contains at least 3 CpG sites.
bedtools intersect -a <(awk 'NR>1{print $1"\t"$2"\t"$3}' p.dm_annotated.txt) -b <(awk 'NR>1{print $1"\t"$2-1"\t"$2}' p_meth.txt) -wa -wb |uniq | cut -f1-3 | uniq -c |awk '{if($1>=3) print $2"\t"$3"\t"$4}' | awk '{FS=OFS="\t"}''NR==FNR{a[$1"_"$2"_"$3]=1}NR>FNR{if($1"_"$2"_"$3 in a )print $0}' - p.dm_annotated.txt  > p.dm_annotated1.txt
bedtools intersect -a <(awk 'NR>1{print $1"\t"$2"\t"$3}' g.dm_annotated.txt) -b <(awk 'NR>1{print $2"\t"$3-1"\t"$3}' g_meth.txt) -wa -wb |uniq | cut -f1-3 | uniq -c |awk '{if($1>=3) print $2"\t"$3"\t"$4}' | awk '{FS=OFS="\t"}''NR==FNR{a[$1"_"$2"_"$3]=1}NR>FNR{if($1"_"$2"_"$3 in a )print $0}' - g.dm_annotated.txt  > g.dm_annotated1.txt
bedtools intersect -a <(awk 'NR>1{print $1"\t"$2"\t"$3}' s.dm_annotated.txt) -b <(awk 'NR>1{print $2"\t"$3-1"\t"$3}' s_meth.txt) -wa -wb |uniq | cut -f1-3 | uniq -c |awk '{if($1>=3) print $2"\t"$3"\t"$4}' | awk '{FS=OFS="\t"}''NR==FNR{a[$1"_"$2"_"$3]=1}NR>FNR{if($1"_"$2"_"$3 in a )print $0}' - s.dm_annotated.txt  > s.dm_annotated1.txt

###solo-WCGW anglysis
#convert CpG site to bedgraph format
awk -F"\t" '{split($1,a,":");print a[1]"\t"a[2]"\t"a[2]"\t"$2}' p_wi38_meth_merge.txt > p_wi38.bedGraph
awk -F"\t" '{split($1,a,":");print a[1]"\t"a[2]"\t"a[2]"\t"$2}' p_esc_meth_merge.txt > p_esc.bedGraph
#obtain the three group site in 0 based
awk -F"\t" 'NR>1{print $1"\t"$2-3"\t"$3+2"\t"$7}' hg38.cpgsite.txt > hg38.cpgsite.dnmt.bed
bedtools getfasta -fi /data/nier/ref/hg38/hg38.fa -bed hg38.cpgsite.dnmt.bed > hg38.cpgsite.dnmt.fa 
paste <(awk 'NR%2==1' hg38.cpgsite.dnmt.fa) <(awk 'NR%2==0' hg38.cpgsite.dnmt.fa) |sed 's/^>//g' > hg38.cpg.dnmt.bed
awk -F"\t" 'NR>1{print $1"\t"$2-2"\t"$3+1"\t"$7}' hg38.cpgsite.txt > hg38.cpgsite.bed
bedtools getfasta -fi /data/nier/ref/hg38/hg38.fa -bed hg38.cpgsite.bed > hg38.cpgsite.fa 
paste <(awk 'NR%2==1' hg38.cpgsite.fa) <(awk 'NR%2==0' hg38.cpgsite.fa) |sed 's/^>//g' > hg38.cpg.bed
awk -F"\t" '{split($1,c,":");split(c[2],d,"-");a=substr($2,1,1);b=substr($2,4,4); if(a ~ /[AT]/ && b ~/[AT]/)print c[1]"\t"d[1]+1"\t"d[1]+2}' hg38.cpg.bed > hg38.wcgw.bed
awk -F"\t" '{split($1,c,":");split(c[2],d,"-");a=substr($2,1,1);b=substr($2,4,4); if(a ~ /[CG]/ && b ~/[CG]/)print c[1]"\t"d[1]+1"\t"d[1]+2}' hg38.cpg.bed > hg38.scgs.bed
awk -F"\t" '{split($1,c,":");split(c[2],d,"-");a=substr($2,1,1);b=substr($2,4,4); if((a ~ /[CG]/ && b ~/[AT]/)|| (b ~ /[CG]/ && a ~/[AT]/))print c[1]"\t"d[1]+1"\t"d[1]+2}' hg38.cpg.bed > hg38.wcgs.bed
#split group
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgw.bed -b p_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tmixture\twcgw\tESC"}' > esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgw.bed -b g_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tG1\twcgw\tESC"}' >> esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgw.bed -b s_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tS\twcgw\tESC"}' >> esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.scgs.bed -b p_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tmixture\tscgs\tESC"}' >> esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.scgs.bed -b g_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tG1\tscgs\tESC"}' >> esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.scgs.bed -b s_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tS\tscgs\tESC"}' >> esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgs.bed -b p_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tmixture\twcgs\tESC"}' >> esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgs.bed -b g_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tG1\twcgs\tESC"}' >> esc_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgs.bed -b s_esc_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tS\twcgs\tESC"}' >> esc_group_cg.txt

bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgw.bed -b p_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tmixture\twcgw\tESC"}' > wi38_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgw.bed -b g_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tG1\twcgw\tESC"}' >> wi38_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgw.bed -b s_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tS\twcgw\tESC"}' >> wi38_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.scgs.bed -b p_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tmixture\tscgs\tESC"}' >> wi38_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.scgs.bed -b g_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tG1\tscgs\tESC"}' >> wi38_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.scgs.bed -b s_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tS\tscgs\tESC"}' >> wi38_group_cg.txt bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgs.bed -b p_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tmixture\twcgs\tESC"}' >> wi38_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgs.bed -b g_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tG1\twcgs\tESC"}' >> wi38_group_cg.txt
bedtools intersect -a /data/nier/ref/hg38/cpg/hg38.wcgs.bed -b s_wi38_filter.bedGraph -wa -wb |awk -F"\t" '{print $1":"$5"\t"$7"\tS\twcgs\tESC"}' >> wi38_group_cg.txt
#0 flanking CpG site at 100 bp
awk 'NR>1' data_meth.txt | cut -f1 | tr ":" "\t" | awk -F"\t" '{print $1"\t"$2-1"\t"$2}' > data_meth.bed
awk -F"\t" '{print $1"\t"$3-50"\t"$3+50}' data_meth.bed |bedtools intersect -a stdin -b data_meth.bed -c | awk '{if($4==1) print $1":"$2+50}' > data_meth_onecg.pos
