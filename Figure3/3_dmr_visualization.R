library("GenomicRanges")
library("IRanges")
library("reshape2")
library("ggplot2")
library("BSgenome.Hsapiens.UCSC.hg38",character.only=TRUE)
library(tidyr)
library("RColorBrewer")
library("gridExtra")
display.brewer.all()
library("ggpubr")
library("ggpval")
library(rstatix)
library(plyr)

coul <- brewer.pal(3, "Pastel2") 


#Fig 3a
esc<-read.table("esc_group_cg.txt",sep="\t",header=F,stringsAsFactors =F )
colnames(esc)<-c("pos","methylation","group","type","cell")
esc2<-dcast(esc,pos+type ~ group,value.var="methylation")
esc3<-esc2[,3:5]
rownames(esc3)<-esc2$pos
cv <- apply(esc3,1,function(x){x<-sd(x)/mean(x)})
cv<-na.omit(cv)
esc5<-subset(esc, pos %in% names(cv))

wi38<-read.table("wi38_group_cg.txt",sep="\t",header=F,stringsAsFactors =F )
colnames(wi38)<-c("pos","methylation","group","type","cell")
wi38$cell<-gsub("ESC","WI38",wi38$cell)
wi382<-dcast(wi38,pos+type ~ group,value.var="methylation")
wi383<-wi382[,3:5]
rownames(wi383)<-wi382$pos
cv <- apply(wi383,1,function(x){x<-sd(x)/mean(x)})
cv<-na.omit(cv)
wi385<-subset(wi38, pos %in% names(cv))
data<-rbind(esc5,wi385)

data<-within(data,type<-factor(type,levels=c("scgs","wcgs","wcgw")))
data$groups<-paste0(data$cell,"_",data$group)
data<-within(data,groups<-factor(groups,levels=c("ESC_mixture","ESC_G1","ESC_S","WI38_mixture","WI38_G1","WI38_S")))
data1<-dcast(data,pos+type ~ groups,value.var="methylation")
write.table(data1,file="data_meth.txt",sep="\t",quote=F,row.names = F)
#The "data_meth.txt" convert to "data_meth_onecg.pos" referred to the '2_process.sh' file.
position<-read.table("data_meth_onecg.pos",header=F,sep="\t")$V1

data2<- subset(data, pos %in% position)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )  
  datac <- rename(datac, c("mean" = measurevar))  
  datac$se <- datac$sd / sqrt(datac$N)    
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult  
  return(datac)
}

colnames(delta_df)<- c("pos","type","esc_delta","wi38_delta")
delta_df1<- na.omit(gather(delta_df,category, value,3:4))
dfc <- summarySE(delta_df1, measurevar="value", groupvars=c("category","type"))
#ESC
esc<- subset(delta_df1, category %in% "esc_delta")
stat.test <- esc %>%
  t_test(value ~ type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
esc<- subset(dfc, category %in% "esc_delta")
esc$type <- toupper(esc$type)
p<- ggplot(esc, aes(x=type, y=value)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge())  +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2,
                position=position_dodge(.9)) + 
    labs(title="hESC",y = "difference methylation value") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p





