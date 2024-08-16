#test_cnv_read_depth.txt is reads depth in windows
#test_cnv_cnvnator_call.txt is the result of CNVnator 
#This script is used to visualize the three samples.
x=read.table("../../data/test_cnv_read_depth.txt")
reg=read.table("../../data/test_cnv_cnvnator_call.txt",sep="\t")
reg1<-subset(reg,reg$V2=="chr19")
reg1$V6<-log(reg1$V5,2)
sample=c("mef_all","esc_30","esc_all")
x1<- data.frame(x[,1:2],log2(x[,(1+2)]/median(x[,(1+2)])),
log2(x[,(2+2)]/median(x[,(2+2)])),log2(x[,(3+2)]/median(x[,(3+2)])))
colnames(x1)<-colnames(x)
i=1
rt=reg1[reg1$V1==sample[i],c(2,3,4,6)]
colnames(rt)<-c("chromosome","start","end","segmean")
data<-x1[,c(1,2,i+2)]
colnames(data)<-c("chromosome","coordinate","cn")
library("GenVisR")
cnView(data, z = rt, chr = "chr19", genome = "mm9", ideogram_txtSize = 4,CNscale="relative")

i=2
rt=reg1[reg1$V1==sample[i],c(2,3,4,6)]
colnames(rt)<-c("chromosome","start","end","segmean")
data<-x1[,c(1,2,i+2)]
colnames(data)<-c("chromosome","coordinate","cn")
library("GenVisR")
cnView(data, z = rt, chr = "chr19", genome = "mm9", ideogram_txtSize = 4,CNscale="relative")

i=3
rt=reg1[reg1$V1==sample[i],c(2,3,4,6)]
colnames(rt)<-c("chromosome","start","end","segmean")
data<-x1[,c(1,2,i+2)]
colnames(data)<-c("chromosome","coordinate","cn")
library("GenVisR")
cnView(data, z = rt, chr = "chr19", genome = "mm9", ideogram_txtSize = 4,CNscale="relative")