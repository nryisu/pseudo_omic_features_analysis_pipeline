library("data.table")
data<-read.table("esc.mm9.rd.gc",sep="\t",stringsAsFactors = F)
data <- data[!(data$V1 %in% c("chrM", "chrX", "chrY")), ]
data$logrd <- log2(data[,4]/median(data[,4])+0.01)

lst<-list()
for (chr in unique(data$V1)){
  chrom <- data[V1==chr]
  gcCount.loess<-loess(logrd~V5,data=chrom,control =loess.control(surface = "direct"),degree=2)
  predictions1<-predict(gcCount.loess,chrom$V5)
  chrom$NRC <- chrom$logrd-predictions1
  #setkey(chrom,V5)
  #gcCount.loess<-loess(NRC~V5,data=chrom,control =loess.control(surface = "direct"),degree=2)
  #predictions<-predict(gcCount.loess,chrom$V5)
  lst[[chr]] <- chrom
}
datanrc <- rbindlist(lst)
saveRDS(datanrc,file="datanrc.rds")

datanrc<- readRDS("datanrc.rds")
dat <- datanrc[,c(1,2,3,7)]
colnames(dat)<-c("chrom","start","end","lasso")
head(dat)
dat$pos<-(dat[,3]+dat[,2])/2
write.table(dat,file="esc_mm9_rd_gclasso.txt",quote=FALSE, sep="\t",row.names = F)

class(dat$pos)
summary(dat$lasso)
library("DNAcopy")
CNA.corrtrd<- CNA(dat$lasso,
                  dat$chrom,dat$pos,
                  data.type="logratio",sampleid="correctedgc-Lasso")
smoothed.CNA.corrtrd <- smooth.CNA(CNA.corrtrd)
sdundo.CNA.corrtrd <- segment(smoothed.CNA.corrtrd,  undo.splits="sdundo")
png("esc_mm9_correctgclasso.png")
plot(sdundo.CNA.corrtrd,plot.type="w",ylim=c(-4,4))
dev.off()
write.table(segments.summary(sdundo.CNA.corrtrd),file="esc_mm9_gclasso.seg",quote=FALSE, sep="\t",row.names = F)


rm(list=ls())
esc<-read.table("esc_mm9_gclasso.seg.ratio1",header=F,sep="\t")
colnames(esc)<-c("chr","start","end","seg","rt")
summary(esc$seg)
data<-subset(esc,seg>0.25)
dat<-subset(esc,seg< -0.25)
library(LSD)
par(mar=c(0.5,0.5,0.5,0.5))
heatscatter(data$seg,data$rt,add.contour=TRUE,xlim=c(-1,2),ylim=c(-1,1))
par(new=T)
heatscatter(dat$seg,dat$rt,add.contour=TRUE,xlim=c(-1,1),ylim=c(-1,1))
abline(h=0,col="red",lwd=2)
abline(v=0,col="red",lwd=2)	 
head(data)
enrichscore<-function(data){
  #rt>0,cnv>0.25
  first <- length(subset(data,seg>0.25 & rt>0 & seg< 2)$chr)
  second <- length(subset(data,seg< -0.25 & rt>0 & seg> -2)$chr)
  third <- length(subset(data,seg< -0.25 & rt<0 & seg> -2)$chr)
  fourth <- length(subset(data,seg>0.25 & rt<0 & seg< 2)$chr)
  first.es<-round((first/second)/((first+fourth)/(second+third)),3)
  second.es<-round((second/first)/((second+third)/(first+fourth)),3)
  third.es<-round((third/fourth)/((second+third)/(first+fourth)),3)
  fourth.es<-round((fourth/third)/((first+fourth)/(second+third)),3)
  result.es<-c(first.es,second.es,third.es,fourth.es)
  names(result.es)<-c("first.es","second.es","third.es","fourth.es")
  result<-c(first,second,third,fourth)
  names(result)<-c("first","second","third","fourth")
  dat<-list(result,result.es)
  return(dat)
}
dat1<-enrichscore(esc)
#plot scatter and line
plot(data$V5,log(data$V4+0.01),cex=0.1,xlab="GC Content",ylab="RD",sep="")))
lines(RC_DT$GC,predictions1,col = "red")

chr19 <- dat[dat$chrom=="chr19",]
par(mfrow=c(1,1))
plot(chr19$lasso~chr19$pos,pch=".",cex=0.2,col="black",ylim=c(-2,2),main="Lasso_correct") 

plot(sdundo.CNA.corrtrd,plot.type="s",ylim=c(-2,2))


#sustract the influence of GC
resi <‐ log(RC_DT$RC+0.01)‐predictions1
RC_DT$RC <‐ resi
setkey(RC_DT,GC)

#plot scatter and line using Norm GC data
plot(RC_DT$GC,RC_DT$RC,cex=0.1,xlab="GC Content",ylab=expression("NRC"["G
C"]))
gcCount.loess <‐ loess(RC~GC,data=RC_DT,control = loess.control(surface ="direct"),degree=2)
predictions2 <‐ predict(gcCount.loess,RC_DT$GC)
lines(RC_DT$GC,predictions2,col="red")