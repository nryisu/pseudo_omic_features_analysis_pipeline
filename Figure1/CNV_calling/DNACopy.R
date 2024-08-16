dat<-read.table("test_depth_corrected.bed",sep="\t",header=F,stringsAsFactors = F)
colnames(dat)<-c("chrom","start","end","rd","logcorrtrd")
dat$pos<-(dat[,3]+dat[,2])/2
#CNV calling
library("DNAcopy")
CNA.corrtrd<- CNA(dat$logcorrtrd,
                  dat$chrom,dat$pos,
                  data.type="logratio",sampleid="correctedRT-RD")
smoothed.CNA.corrtrd <- smooth.CNA(CNA.corrtrd)
sdundo.CNA.corrtrd <- segment(smoothed.CNA.corrtrd,  undo.splits="sdundo")
png(paste(Args[6],"_correctRTRD.png",sep=""))
plot(sdundo.CNA.corrtrd,plot.type="w",ylim=c(-4,4))
dev.off()
write.table(segments.summary(sdundo.CNA.corrtrd),file=paste(Args[6],".seg",sep=""),quote=FALSE, sep="\t",row.names = F)