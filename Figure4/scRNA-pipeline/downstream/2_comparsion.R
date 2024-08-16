data<-load("data.RData")
input <- lapply(data,rownames)
names(input)<-c("G1","S","G2M","whole")
venn(input)
data<-list(G,S,M)
names(data)<-c("G1","S","G2M")
input <- lapply(data,rownames)
sapply(input,length)
names(input)<-c("G1","S","G2M")
venn(input)
deg<-unique(as.vector(c(input$G1,input$S,input$G2M))) 
length(deg)#2614
input <- list(deg,rownames(whole))
names(input)<-c("cycle_DEG","whole")
venn(input)
overlap<-Reduce(intersect, input) 
length(overlap) #2177
library(gplots)
library(VennDiagram)
dim(whole)#5398
venn.plot <- draw.pairwise.venn(area1=c(5398),area2 =c(2614),cross.area =c(2177),
				category = c("whole", "cycle"),fill = c("navy", "red"),
				lty = "blank",cex = 2, cat.cex = 2, cat.pos = c(180, 180),
				 cat.dist = 0.05,cat.just = list(c(0, 1), c(1, 1)))
grid.draw(venn.plot)
whole$gene<-rownames(whole)
whole_common<-subset(whole,gene %in% overlap)
dim(whole_common)#2177/5398 0.40 (60%的差异基因为假阳性。丢失了约437/5398=8.1%的基因)。
#对两者进行富集分析
whole_fp<-subset(whole,!(gene %in% overlap))
dim(whole_fp) #3221
save(whole_fp,file="whole_fp.RData")
save(whole_common,file="whole_common.RData")
save(deg,file="deg.RData")

library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
library(pathview)
whole.go <- enrichGO(gene = whole$gene,OrgDb = 'org.Mm.eg.db', 
	keyType = 'SYMBOL',ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.10,
	minGSSize = 3, maxGSSize = 500)
write.table(whole.go,file = "whole_go_bp.txt",quote = FALSE,sep = "\t")

whole_common.go <- enrichGO(gene = whole_common$gene,OrgDb = 'org.Mm.eg.db', 
		keyType = 'SYMBOL',ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.10,
		minGSSize = 3, maxGSSize = 500)
#dotplot(whole_common.go, showCategory=10)
write.table(whole_common.go,file = "whole_common_go_bp.txt",quote = FALSE,sep = "\t")

deg.go <- enrichGO(gene = deg,OrgDb = 'org.Mm.eg.db', 
		keyType = 'SYMBOL',ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.10,
		minGSSize = 3, maxGSSize = 500)
#dotplot(deg.go, showCategory=10)
write.table(deg.go,file = "deg_go_bp.txt",quote = FALSE,sep = "\t")

source("D:/18-11/mesc/com/go_slim.R")
library(GOstats) 
library(GO.db)
goterms <- Term(GOTERM) 
go_slim<-read.table("D:/18-11/mesc/com/goslim_term",header=F,stringsAsFactors=F)
go_slim<-as.vector(go_slim$V1)
go_term<- whole.go$ID
slim_stat_bp <- goSlim_fct(sample_go=go_term, gotype="BP")
write.table(slim_stat_bp,file = "whole_goslim.txt",quote = FALSE,sep = "\t")
go_term<- whole_common.go$ID

#######bar plot
dot<-read.table("cycle.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
dim(dot)
head(dot)
library(ggplot2) 
dot$GO<-factor(dot$Gene_set,levels=dot[order(dot$cluster,-dot$FDR_q.val), ]$Gene_set)
dot$fdr<--log10(dot$FDR_q.val)
ggplot(dot, aes(x=GO,fill=cluster,y=ifelse(test = cluster == "cycle",yes= -fdr,no = fdr)))+
	geom_bar(stat='identity',width=0.75)+
	scale_y_continuous(labels = abs, limits = max(dot$fdr) * c(-1,1)) +
	scale_fill_manual(labels = c("whole","cycle" ),values = c("cycle"="red","whole"="navy" )) +
   	labs(title="cycle processed related GO BP terms",y = "-log10(p.adjust.val)") +   
	theme(panel.grid.major=element_blank(),panel.background=element_blank(),
		panel.border=element_blank(),panel.grid.minor=element_blank())+
	theme(axis.line.x = element_line(color="black", size = .1),
       	 axis.line.y = element_line(color="black", size = .1)) +
	coord_flip()

setwd("E:/技能/百度云同步盘/cell_cycle/cell_cycle_de/rnaseqdata/mouse/mesc_npc")
dot<-read.table("cycle_simp.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
dim(dot)
head(dot)
library(ggplot2) 
dot<-dot[order(dot$cluster,-dot$FDR_q.val), ]
dot$fdr<--log10(dot$FDR_q.val)
ggplot(dot, aes(x=Gene_set,fill=cluster,y=ifelse(test = cluster == "cycle",yes= -fdr,no = fdr)))+
	geom_bar(stat='identity',width=0.75)+
	scale_y_continuous(labels = abs, limits = max(dot$fdr) * c(-1,1)) +
	scale_fill_manual(labels = c("whole","cycle" ),values = c("cycle"="red","whole"="navy" )) +
   	labs(title="cycle processed related GO BP terms",y = "-log10(p.adjust.val)") +   
	theme(panel.grid.major=element_blank(),panel.background=element_blank(),
		panel.border=element_blank(),panel.grid.minor=element_blank())+
	theme(axis.line.x = element_line(color="black", size = .1),
       	 axis.line.y = element_line(color="black", size = .1)) +
	coord_flip()	
	

dot<-read.table("development.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
dim(dot)
head(dot)
library(ggplot2) 
dot$GO<-factor(dot$Gene_set,levels=dot[order(dot$cluster,-dot$FDR_q.val), ]$Gene_set)
dot$fdr<--log10(dot$FDR_q.val)
ggplot(dot, aes(x=GO,fill=cluster,y=ifelse(test = cluster == "cycle",yes= -fdr,no = fdr)))+
	geom_bar(stat='identity',width=0.75)+
	scale_y_continuous(labels = abs, limits = max(dot$fdr) * c(-1,1)) +
	scale_fill_manual(labels = c("cycle", "whole"),values = c("cycle"="red","whole"="navy" )) +
   	labs(title="development processed related GO BP terms",y = "-log10(p.adjust.val)") +   
	theme(panel.grid.major=element_blank(),panel.background=element_blank(),
		panel.border=element_blank(),panel.grid.minor=element_blank())+
	theme(axis.line.x = element_line(color="black", size = .1),
       	 axis.line.y = element_line(color="black", size = .1)) +
	coord_flip()

dot<-read.table("development_simp.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
dim(dot)
library(ggplot2) 
dot<-dot[order(dot$cluster,-dot$FDR_q.val), ]
dot$fdr<--log10(dot$FDR_q.val)
dot$GO<- factor(dot$Gene_set, levels = dot$Gene_set) 
ggplot(dot, aes(x=Gene_set,fill=cluster,y=ifelse(test = cluster == "cycle",yes= -fdr,no = fdr)))+
	geom_bar(stat='identity',width=0.75)+
	scale_y_continuous(labels = abs, limits = max(dot$fdr) * c(-1,1)) +
	scale_fill_manual(labels = c("cycle", "whole"),values = c("cycle"="red","whole"="navy" )) +
   	labs(title="development processed related GO BP terms",y = "-log10(p.adjust.val)") +   
	theme(panel.grid.major=element_blank(),panel.background=element_blank(),
		panel.border=element_blank(),panel.grid.minor=element_blank())+
	theme(axis.line.x = element_line(color="black", size = .1),
       	 axis.line.y = element_line(color="black", size = .1)) +
	coord_flip()
	
	
dot<-read.table("development_simp.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
dim(dot)
library(ggplot2) 
dot$fdr<--log10(dot$FDR_q.val)
## extract duplicate elements
d <- duplicated(dot$Gene_set)
d1<-dot[d,]$Gene_set
dup<-subset(dot, Gene_set %in% d1 )
un<-subset(dot, !(Gene_set %in% d1 ))
dup<-dup[order(dup$cluster,-dup$FDR_q.val), ]
dup$GO<- factor(dup$Gene_set, levels = unique(dup$Gene_set) )
un$GO<-factor(un$Gene_set,levels=un[order(un$cluster,-un$FDR_q.val), ]$Gene_set)
dot<-rbind(dup,un)
ggplot(dot, aes(x=GO,fill=cluster,y=ifelse(test = cluster == "cycle",yes= -fdr,no = fdr)))+
  geom_bar(stat='identity',width=0.75)+
  scale_y_continuous(labels = abs, limits = max(dot$fdr) * c(-1,1)) +
  scale_fill_manual(labels = c("cycle", "whole"),values = c("cycle"="red","whole"="navy" )) +
  labs(title="development processed related GO BP terms",y = "-log10(p.adjust.val)") +   
  theme(panel.grid.major=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line.x = element_line(color="black", size = .1),
        axis.line.y = element_line(color="black", size = .1)) +
  coord_flip()

  
dot<-read.table("cycle_simp.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
dim(dot)
library(ggplot2) 
dot$fdr<--log10(dot$FDR_q.val)
## extract duplicate elements
d <- duplicated(dot$Gene_set)
d1<-dot[d,]$Gene_set
dup<-subset(dot, Gene_set %in% d1 )
un<-subset(dot, !(Gene_set %in% d1 ))
dup<-dup[order(dup$cluster,-dup$FDR_q.val), ]
dup$GO<- factor(dup$Gene_set, levels = unique(dup$Gene_set) )
un$GO<-factor(un$Gene_set,levels=un[order(un$cluster,-un$FDR_q.val), ]$Gene_set)
dot<-rbind(dup,un)
ggplot(dot, aes(x=GO,fill=cluster,y=ifelse(test = cluster == "cycle",yes= -fdr,no = fdr)))+
  geom_bar(stat='identity',width=0.75)+
  scale_y_continuous(labels = abs, limits = max(dot$fdr) * c(-1,1)) +
  scale_fill_manual(labels = c("whole","cycle" ),values = c("cycle"="red","whole"="navy" )) +
  labs(title="cycle processed related GO BP terms",y = "-log10(p.adjust.val)") +   
  theme(panel.grid.major=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.line.x = element_line(color="black", size = .1),
        axis.line.y = element_line(color="black", size = .1)) +
  coord_flip()	

library(data.table)
###bar BP num
dat<-data.table(c(84,NA,NA),c(122,44,58),c(61,44,52))
colnames(dat)<-c("common","whole","cycle")
barplot(as.matrix(dat), main="whole vs cycle BP NUM",
  		col=c("cornflowerblue","darkgoldenrod1","slategray2"),
  		legend = c("ALL","whole_cycle_same","common_same")
		, beside=TRUE)


###compare cycle or development process		
dat<-data.table(c(21,4),c(10,10))
colnames(dat)<-c("whole","cycle")
barplot(as.matrix(dat), main="whole vs cycle development BP",
  		col=c("darkblue","slategray2"),
  		legend = c("whole","common")
		, beside=TRUE)

dat<-data.table(c(11,0),c(0,0))
colnames(dat)<-c("whole","cycle")
barplot(as.matrix(dat), main="whole vs cycle cycle BP",
  		col=c("darkblue","slategray2"),
  		legend = c("whole","common")
		, beside=TRUE)

##slim
dat<-data.table(c(11,5,4),c(11,11,8),c(13,14,9))
colnames(dat)<-c("cell proliferation",
	"cell differentiation","system development")
barplot(as.matrix(dat), main="whole vs cycle BP NUM",
  		col=c("darkblue","red","darkgoldenrod1"),
  		legend = c("whole","cycle","same")
		, beside=TRUE)


all<-read.delim("D:/cell_cycle_de/rnaseqdata/classify/sen_table.txt",header=F,sep="\t")
colnames(all)<-c("stage","software","sample","sensitivity")
head(all)
all<-within(all,sensitivity<-round(sensitivity,3))
all<-within(all,stage<-factor(stage,levels=c("G1","S","G2M")))
library(ggplot2)
p<-ggplot(all, aes(x = software, y = sensitivity, fill = sample)) + 
	geom_bar(stat = 'identity', position = 'stack') + 
	facet_grid(~ stage)+ 
	geom_text(aes(label=sensitivity), position = position_stack(vjust = 0.5),color="white", size=3.5)+ 
	scale_color_grey() + 
	theme_classic()+
	scale_fill_manual(values=c("#3399FF","#FF66CC","#33CC99"))
p+theme(legend.position="top")+
	labs(x="software",title="sensitivity",fill= "sample")