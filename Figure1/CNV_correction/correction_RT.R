Args <- commandArgs()
#mesc.rt.uniq.txt
normalizeRT <-function(dat){
  medians <- data.frame()
  overallcountsmedian <- median(dat$rd)
  thislength <- length(dat$rd)
  corr_rc <- vector(length=thislength)
  
  for(i in seq(0.00,1,0.01)){
    j <-which(round(dat$normrt, digits=2) == round(i, digits=2))
    numberofobs <- length(dat$rd[j])
    thismedian <- median(dat$rd[j])
    if (is.na(thismedian)){thismedian <- 0}
    if (thismedian == 0){weight <- 0}else{weight <- overallcountsmedian / thismedian}
    medians <- rbind(medians, t(c(i, thismedian, weight, numberofobs)))
    corr_rc[j] <- round(dat$rd[j] * weight)
  }
  
  dat$corrrtrd <- corr_rc
  return(dat)
}

rt<-read.table(Args[6],sep="\t",header=F,stringsAsFactors = F)
colnames(rt)<-c("chrom","start","end","rd","log2rt")
rt$normrt<-round(2^(rt$log2rt)/max(2^(rt$log2rt)),2)
dat<-normalizeRT(rt)
dat$pos<-(dat[,3]+dat[,2])/2
dat$logcorrtrd<-log2(dat[,7]/median(dat[,7]))
write.table(dat,paste(Args[6],"_depth_corrected.bed",sep=""),quote=FALSE, sep="\t",row.names = F)

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
#Then filter CNV from 'seg.mean' column.
