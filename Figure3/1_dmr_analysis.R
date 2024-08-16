####convert bam to CpG site file
library(methylKit)
mySaveFolder="/data/nier/cycle/RRBS/data"
my.meth1=processBismarkAln(location = "h9_p1.bam",sample.id="h9_p_1", assembly="hg38",read.context="CpG",save.folder=mySaveFolder)
my.meth2=processBismarkAln(location = "h9_p2.bam",sample.id="h9_p_2", assembly="hg38",read.context="CpG",save.folder=mySaveFolder)
my.meth1=processBismarkAln(location = "wi38_p1.bam",sample.id="wi38_p_1", assembly="hg38",read.context="CpG",save.folder=mySaveFolder)
my.meth2=processBismarkAln(location = "wi38_p2.bam",sample.id="wi38_p_2", assembly="hg38",read.context="CpG",save.folder=mySaveFolder)

###merge replicate CpG site file 
library(methylKit)
myobj=methRead(list("h9_p_1_CpG.txt","h9_p_2_CpG.txt"),sample.id=list("h9_p_1","h9_p_2"),assembly="hg38",treatment=c(0,0),context="CpG")
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
getMethylationStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE) 
meth=unite(filtered.myobj, destrand=FALSE)
dim(meth)#2249349
write.table(as.data.frame(meth),file = "p_esc_meth.txt",quote = FALSE,sep = "\t",row.names=F)

myobj=methRead(list("wi38_p_1_CpG.txt","wi38_p_2_CpG.txt"),sample.id=list("wi38_p_1","wi38_p_2"),assembly="hg38",treatment=c(1,1),context="CpG")
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
getMethylationStats(filtered.myobj[[1]],plot=TRUE,both.strands=FALSE) 
meth=unite(filtered.myobj, destrand=FALSE)
dim(meth)#2126654
write.table(as.data.frame(meth),file = "p_wi38_meth.txt",quote = FALSE,sep = "\t",row.names=F)

p.esc<-read.table("p_esc_meth.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
p.wi38<-read.table("p_wi38_meth.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
p.esc$names<- paste(p.esc[,1],p.esc[,2],sep=":")
p.wi38$names<- paste(p.wi38[,1],p.wi38[,2],sep=":")
whole.esc<-data.frame("names"=p.esc$names,"esc.whole"=round((p.esc[,6]+p.esc[,9])/(p.esc[,5]+p.esc[,8]),3))
whole.wi38<-data.frame("names"=p.wi38$names,"wi38.whole"=round((p.wi38[,6]+p.wi38[,9])/(p.wi38[,5]+p.wi38[,8]),3))
write.table(whole.esc,file = "p_esc_meth_merge.txt",quote = FALSE,sep = "\t",row.names=F,col.names=F)
write.table(whole.wi38,file = "p_wi38_meth_merge.txt",quote = FALSE,sep = "\t",row.names=F,col.names=F)

###obtain DMR 
library(methylKit)
myobj=methRead(list("h9_p_1_CpG1.txt","h9_p_2_CpG1.txt","wi38_p_1_CpG1.txt","wi38_p_2_CpG1.txt"),sample.id=list("h9_p_1","h9_p_2","wi38_p_1","wi38_p_2"),assembly="hg38",treatment=c(0,0,1,1),context="CpG",mincov=10)
#filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000)
meth=unite(tiles, destrand=FALSE)
write.table(as.data.frame(meth),file = "p_meth.txt",quote = FALSE,sep = "\t",row.names=F)
myDiff=calculateDiffMeth(meth,num.cores=4)
#get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=30,qvalue=0.01)
write.table(x = as.data.frame(myDiff25p),file = "p_dmr_25p.txt", sep = '\t', quote = F,col.names = NA)

###DMRs annotation
library("annotatr",lib="/data/nier/bin/R3.6")
library("TxDb.Hsapiens.UCSC.hg38.knownGene",lib="/data/nier/bin/R3.6")
p_diff<-read.table("p_dmr_25p.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
p_diff$names<-paste(p_diff$chr,p_diff$start,p_diff$end,sep="-")
p<-read.table("p_meth_unit.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
p$names<-paste(p[,1],p[,2],p[,3],sep="-")
p_me<-data.frame("names"=p$names,"esc.whole"=round((p[,6]+p[,9])/(p[,5]+p[,8]),3),"wi38.whole"=round((p[,12]+p[,15])/(p[,11]+p[,14]),3))
dat<-merge(p_diff,p_me,by="names",all=TRUE)
dat<-na.omit(dat)
dat$DM_status<-with(dat,ifelse(meth.diff>0,"hyper",ifelse(meth.diff<0,"hypo","none")))
p<-dat[,c(3,4,5,12,7,6,9,10,11)]
extraCols = c(meth.diff = 'numeric', esc.whole = 'numeric', wi38.whole = 'numeric')
write.table(p,file = "p_anno.txt",quote = FALSE,sep = "\t",col.names=F,row.names=F)
p.dm_regions = read_regions(con = "p_anno.txt", genome = 'hg38', extraCols = extraCols, format = 'bed',rename_name = 'DM_status', rename_score = 'pval')
table(p.dm_regions$DM_status)
# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic',
           'hg38_genes_intronexonboundaries')
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)
# Intersect the regions we read in with the annotations
p.dm_annotated = annotate_regions(regions = p.dm_regions,annotations = annotations,ignore.strand = TRUE,quiet = FALSE,minoverlap = 500)
# A GRanges object is returned
#print(p.dm_annotated)
# Coerce to a data.frame
p.df_dm_annotated = data.frame(p.dm_annotated)
# See the GRanges column of dm_annotaed expanded
#print(head(p.df_dm_annotated))
length(unique(p.df_dm_annotated$annot.symbol))#10392
# Subset based on a gene symbol, in this case NOTCH1
notch1_subset = subset(p.df_dm_annotated, annot.symbol == 'NOTCH1')
p.dm_annsum = summarize_annotations(annotated_regions = p.dm_annotated,quiet = TRUE)
# annot.type   
write.table(as.data.frame(p.dm_annotated),file = "p.dm_annotated.txt",quote = FALSE,sep = "\t",row.names=F,col.names=T)

###predict TFBS by i-cisTarget,
#bed region upload in icis-Target web: https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget

###obtain CpG site in reference
# BiocManager::install("RnBeads")
# BiocManager::install("RnBeads.hg38")
library("RnBeads")
library("RnBeads.hg38")
cpg<-rnb.get.annotation('hg38',type="CpG")
cpg1<-lapply(cpg,as.data.frame)
cpg2<-do.call(rbind.data.frame, cpg1)
cpg3<-subset(cpg2,strand %in% "+")
write.table(file="hg38.cpgsite.txt",cpg3[,c(1,2,3,5,6,7,8)],quote=F,col.names=T,row.names=F,sep="\t") 

