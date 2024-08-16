normalizeGC <-
function(dat){
   
  medians <- data.frame()
  overallcountsmedian <- median(dat$rd)
  overallgcmedian <- median(dat$GCpercent)
     
  thislength <- length(dat$rd)
    
  corr_rc <- vector(length=thislength)
  for(i in seq(0.00,1,0.01)){
    j <- (round(dat$GCpercent, digits=2) == round(i, digits=2))
    numberofobs <- length(dat$rd[j])
    thismedian <- median(dat$rd[j])
    if (is.na(thismedian)){thismedian <- 0}
    if (thismedian == 0){weight <- 0}else{weight <- overallcountsmedian / thismedian}
    medians <- rbind(medians, t(c(i, thismedian, weight, numberofobs)))
   corr_rc[j] <- round(dat$rd[j] * weight)
  }
  
  dat$Corrrd <- corr_rc
  return(dat)
}
#"esc.mm9.rd.gc is the reads depth file.
dat<-read.table("esc.mm9.rd.gc",sep="\t",header=F,stringsAsFactors = F)
colnames(dat)<-c("chrom","start","end","rd","GCpercent")
dat<-within(dat,GCpercent<-round(GCpercent,2))
dat1<-normalizeGC(dat)
write.table(dat1,file = "esc.mm9.correct.rd.txt", sep = '\t', quote = F,col.names = T,row.names=F)