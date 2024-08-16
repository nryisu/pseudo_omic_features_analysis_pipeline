#CNV in repeatmask
lst <-list()
file_list <- list.files("./", pattern = "\\.rm$")

for (file in file_list) {
  repat_mask <-list()
  esc<-read.table(file,header=F,sep="\t")
  colnames(esc)<-c("chr","start","end","seg","rt","rm")
  esc1<-subset(esc[,c(1,2,3,4,5)],abs(seg)>0.25 & abs(seg)<2)
  esc1 <- distinct(esc1)
  esc1 <- subset(esc1,(esc1[,3]-esc1[,2])>=50000)
  lst[[file]] <- esc1
}
lst  
names(lst)<-c("correctedGC","original","correctedRT","Lasso_GC")

dat22 <- do.call(rbind, lapply(seq_along(lst), function(i) {
  df <- lst[[i]]
  df$methods <- names(lst)[i]
  return(df)
}))
table(dat22$methods)
#plot
library(LSD)

pdf(file = "scatter_repeatmask_filter.pdf", width = 6, height = 6)
par(mfrow=c(2,2))
for (group in c("original","correctedGC","Lasso_GC","correctedRT")) {
  group_data <- subset(dat22, methods == group)
  data<-subset(group_data,seg>0.25 & seg< 2)
  dat<-subset(group_data,seg< -0.25 & seg> -2)
  par(mar=c(0.5,0.5,0.5,0.5))
  heatscatter(data$seg,data$rt, add.contour=TRUE,pch=".",xlim=c(-2,2),ylim=c(-2,2),main = group)
  par(new=T)
  heatscatter(dat$seg,dat$rt,add.contour=TRUE,pch=".",xlim=c(-2,2),ylim=c(-2,2), 
              main = group)
  abline(h=0,col="red",lwd=2)
  abline(v=0,col="red",lwd=2)	 
  cat(file,group)
  print(dim(group_data[,1:3][!duplicated(group_data[1:3,]), ])[1])
  }
dev.off()



lst <-list()
file_list <- list.files("./", pattern = "\\.rm$")

for (file in file_list) {
  esc<-read.table(file,header=F,sep="\t")
  colnames(esc)<-c("chr","start","end","seg","rt")
  esc1<-subset(esc[,c(1,2,3,4,5)],abs(seg)>0.25 & abs(seg)<2)
  esc1 <- distinct(esc1)
  esc1 <- subset(esc1,(esc1[,3]-esc1[,2])>=50000)
  dat1<-enrichscore(esc1)
  lst[[file]] <- dat1
}

lst  
dat2<-lapply(lst,function(x){x<-x[1][[1]]})
dat22<-do.call(rbind,dat2)
rownames(dat22)<-c("GC","Not","RT","Lasso")

p_cal1 <- function(row) {
  pvalue <- binom.test((as.numeric(row[1])+as.numeric(row[3])),
                        (as.numeric(row[1])+as.numeric(row[2])+as.numeric(row[4])+as.numeric(row[3])),0.5,alternative = "greater")
  return(pvalue)
}
apply(dat22,1,p_cal1)


###CNV in CpG Islands
lst<-list()
file_list <- list.files("./", pattern = "\\.cpg$")
for (file in file_list) {
  esc<-read.table(file,header=F,sep="\t")
  colnames(esc)<-c("chr","start","end","seg","rt")
  esc1<-subset(esc[,c(1,2,3,4,5)],abs(seg)>0.25 & abs(seg)<2)
  esc1 <- distinct(esc1[,c(1,2,3,4,5)])
  esc1 <- subset(esc1,(esc1[,3]-esc1[,2])>=50000)
  print(dim(esc1))
  lst[[file]] <- esc1
}

lst
names(lst)<-c("correctedGC","original","correctedRT","Lasso_GC")

dat22 <- do.call(rbind, lapply(seq_along(lst), function(i) {
  df <- lst[[i]]
  df$methods <- names(lst)[i]
  return(df)
}))

pdf(file = "scatter_cpg_filter.pdf", width = 6, height = 6)
par(mfrow=c(2,2))
for (group in c("original","correctedGC","Lasso_GC","correctedRT")) {
  group_data <- subset(dat22, methods == group)
  data<-subset(group_data,seg>0.25 & seg< 2)
  dat<-subset(group_data,seg< -0.25 & seg> -2)
  par(mar=c(0.5,0.5,0.5,0.5))
  heatscatter(data$seg,data$rt, add.contour=TRUE,pch=".",xlim=c(-2,2),ylim=c(-2,2),main = group)
  par(new=T)
  heatscatter(dat$seg,dat$rt,add.contour=TRUE,pch=".",xlim=c(-2,2),ylim=c(-2,2), 
              main = group)
  abline(h=0,col="red",lwd=2)
  abline(v=0,col="red",lwd=2)	 
  cat(file,group)
  print(dim(group_data[,1:3][!duplicated(group_data[1:3,]), ])[1])
  
}
dev.off()

lst<-list()

for (file in file_list) {
  esc<-read.table(file,header=F,sep="\t")
  colnames(esc)<-c("chr","start","end","seg","rt")
  esc1<-subset(esc[,c(1,2,3,4,5)],abs(seg)>0.25 & abs(seg)<2)
  esc1 <- distinct(esc1[,c(1,2,3,4,5)])
  esc1 <- subset(esc1,(esc1[,3]-esc1[,2])>=50000)
  dat1<-enrichscore(esc1)
  lst[[file]] <- dat1
}
lst  


##CNV number
dat2<-lapply(lst,function(x){x<-x[1][[1]]})
dat22<-do.call(rbind,dat2)
rownames(dat22)<-c("GC","Not","RT","Lasso")
apply(dat22,1,p_cal1)

library(ggplot2)
library(viridis)
library(reshape2)
dat22_long<-melt(dat22,measure.vars = c('first','second','third','fourth'),value.name='number')
dat22_long$Var1<-factor(dat22_long$Var1,levels=c("Not","GC","Lasso","RT"))

ggplot(dat22_long, aes(fill=Var2, y=number, x=Var1)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = T) +
  geom_text(
    aes(label = number),
    colour = "white", size = 3,
    vjust = 1, position = position_dodge(.9)
  ) +
  ggtitle("CpG number") +
  theme_classic()+
  xlab("")

#CNV enrich_score
dat2<-lapply(lst,function(x){x1<-x[2][[1]]})
dat2
dat22<-do.call(rbind,dat2)
rownames(dat22)<-c("GC","Not","RT","Lasso")
dat22_long<-melt(dat22,measure.vars = c('first.es','second.es','third.es','fourth.es'),
                 value.name='EnrichScore')
dat22_long$Var1<-factor(dat22_long$Var1,levels=c("Not","GC","Lasso","RT"))

ggplot(dat22_long, aes(fill=Var2, y=EnrichScore, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  geom_text(
    aes(label = EnrichScore),
    colour = "white", size = 3,
    vjust = 1, position = position_dodge(.9)
  ) +
  ggtitle("CpG EnrichScore") +
  theme_classic()+
  xlab("")
#visualization Figure 1g
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
library(dplyr)
##lasso
dir()
esc<-read.table("esc_mm9_gclasso.seg.ratio1",header=F,sep="\t")
colnames(esc)<-c("chr","start","end","seg","rt")
summary(esc$seg)
esc1<-subset(esc[,c(1,2,3,4,5)],abs(seg)>0.25 & abs(seg)<2)
esc1 <- distinct(esc1)
esc <- subset(esc1,(esc1[,3]-esc1[,2])>=50000)
pdf(file = "scatter_lasso.pdf", width = 6, height = 6)
data<-subset(esc,seg>0.25 & seg< 2)
dat<-subset(esc,seg< -0.25 & seg> -2)
library(LSD)
par(mar=c(0.5,0.5,0.5,0.5))
heatscatter(data$seg,data$rt,add.contour=TRUE,pch=".",xlim=c(-2,2),ylim=c(-2,2))
par(new=T)
heatscatter(dat$seg,dat$rt,add.contour=TRUE,pch=".",xlim=c(-2,2),ylim=c(-2,2))
abline(h=0,col="red",lwd=2)
abline(v=0,col="red",lwd=2)	 
dev.off()
dat1<-enrichscore(esc)

