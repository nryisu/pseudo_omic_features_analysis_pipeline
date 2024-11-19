library(Seurat)

cal<-function(data){
    lapply(data,function(x){
        x1<-table(x)
        G1.TPR=x1[1,1]/sum(x1[,1])
        G1.FPR=sum(x1[1,-1])/(sum(x1[,-1])+sum(x1[1,-1]))
        G1.F1=2*x1[1,1]/(sum(x1[,1])+sum(x1[1,]))
        #S
        M.TPR=x1[2,2]/sum(x1[,2])
        M.FPR=sum(x1[2,-2])/(sum(x1[,-2])+sum(x1[2,-2]))
        M.F1=2*x1[2,2]/(sum(x1[,2])+sum(x1[2,]))
        #S
        S.TPR=x1[3,3]/sum(x1[,3])
        S.FPR=sum(x1[3,-3])/(sum(x1[,-3])+sum(x1[3,-3]))
    # S.prescion=x1[3,3]/sum(x1[3,])
        S.F1=2*x1[3,3]/(sum(x1[,3])+sum(x1[3,]))
        #S.F11=2 * S.TPR * S.prescion/(S.TPR + S.prescion)
        x11<-data.frame(G1.TPR,G1.FPR,G1.F1,S.TPR,S.FPR,S.F1,M.TPR,M.FPR,M.F1)
        return(x11)
  })}

#cyclone-merge
raw<-list()
re.cyclone.merge<-list()
#90%: 35 *0.9=31.5
for(i in 1:20){
  raw[[i]]<-count[,sample(colnames(count), 32, replace=F)]
  sce <- SingleCellExperiment(assays = list(counts =as.matrix(raw[[i]])))
  assignments <- cyclone(sce,mm.pair, gene.names=rownames(counts(sce)))
  re.cyclone.merge[[i]]<-data.frame(assignments$phases,substr(colnames(raw[[i]]),4,5))
  colnames(re.cyclone.merge[[i]])<-c("pre","true")
}
cal.cyclone.merge<-do.call(rbind,cal(re.cyclone.merge))
library(matrixStats)
cal.cyclone.merge.ave<-round(colMeans(cal.cyclone.merge),3)
#cycleone-raw
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
raw<-list()
re.cyclone.raw<-list()
#90%: 35 *0.9=31.5
for(i in 1:20){
  raw[[i]]<-count[,sample(colnames(count), 32, replace=F)]
  sce <- SingleCellExperiment(assays = list(counts =as.matrix(raw[[i]])))
  assignments <- cyclone(sce,mm.pairs, gene.names=rownames(counts(sce)))
  re.cyclone.raw[[i]]<-data.frame(assignments$phases,substr(colnames(raw[[i]]),4,5))
  colnames(re.cyclone.raw[[i]])<-c("pre","true")
}
cal.cyclone.raw<-do.call(rbind,cal(re.cyclone.raw))
cal.cyclone.raw.ave<-round(colMeans(cal.cyclone.raw),3)
##reCAT
##in windows
setwd("reCAT-master/R")
count<-read.table("../data/mesc35.counts.predict" ,header=TRUE,sep="\t",stringsAsFactors = FALSE,row.names=1)
source("get_test_exp.R")
source("get_ordIndex.R")
source("get_score.R")
source("get_hmm.R")
source("plot.R")
load("../data/ola_mES_2i.RData")
library(edgeR)
dge <- DGEList(counts = count)
dge <- edgeR::calcNormFactors(dge)
cd <- log2(cpm(dge)+1)
cycle.genes<-rownames(t(test_exp))
cycle.genesok<-intersect(cycle.genes,rownames(cd))
length(cycle.genesok)#312
##sample
raw<-list()
re.recat<-list()
#90%: 35 *0.9=31.5
for(i in 1:20){
  raw[[i]]<-count[,sample(colnames(count), 32, replace=F)]
  dge <- DGEList(counts = raw[[i]])
  dge <- edgeR::calcNormFactors(dge)
  cd <- log2(cpm(dge)+1)
  cd.exp<-cd[cycle.genesok,]
  cd.exp<-t(cd.exp)
  dim(cd.exp)
  set.seed(122)
  cd.exp <- get_test_exp(t(cd.exp))
  ordIndex <- get_ordIndex(cd.exp, 1)
  score_result <- get_score(t(cd.exp))
  #plot_bayes(score_result$bayes_score, ordIndex)
  #table(apply(score_result$bayes_score,1,which.max))
  myord = c(1:32)
  rdata = t(data.frame(c(2,8),c(13,19),c(24,30)))
  colnames(rdata)<-c("start","end")
  rownames(rdata)<-NULL
  rdata
  hmm_result <- get_hmm_order(bayes_score=score_result$bayes_score, 
                              mean_score=score_result$mean_score,ordIndex = ordIndex,
                              cls_num = 3, myord = myord, rdata = rdata)
  #plot_bayes(score_result$bayes_score, ordIndex, cls_result = hmm_result, cls_ord = myord, colorbar = 1)
  #table(hmm_result)
  re.recat[[i]]<-data.frame(hmm_result,substr(colnames(cd[,ordIndex]),4,5))
  colnames(re.recat[[i]])<-c("pre","true")
}
saveRDS(re.recat,file="re.recat.rds")
cal.recat<-do.call(rbind,cal(re.recat))
cal.recat.ave<-round(colMeans(cal.recat),3)
saveRDS(cal.recat,file="cal.recat.rds")
#seurat
cd<-count
colnames(cd)<- gsub("ES_", "", colnames(cd), fixed = T)
#ensemblTogenenames
library("biomaRt")
mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
bm.query <- getBM(values=rownames(cd),attributes=c("ensembl_gene_id","external_gene_name"),filters=c("ensembl_gene_id"),mart=mart)
genes<-bm.query[match(rownames(cd), bm.query$ensembl_gene_id),]$external_gene_name
genes.na<-data.frame(genes,rownames(cd),stringsAsFactors=FALSE)
colnames(genes.na)<-c("names","id")
ind <- which(is.na(genes.na$names ))
genes.na[ind, "names"] <- genes.na[ind, "id"]
rownames(cd) <- make.unique(as.character(genes.na$names))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
length(s.genes)#43
length(g2m.genes) #54
s.genes<-convertHumanGeneList(s.genes)
s.genes<-intersect(s.genes,rownames(cd))
length(s.genes) #40
g2m.genes<-convertHumanGeneList(g2m.genes)
g2m.genes<-intersect(g2m.genes,rownames(cd))
length(g2m.genes) #52
##sample
raw<-list()
re.seurat<-list()
#90%: 35 *0.9=31.5
for(i in 1:20){
  raw[[i]]<-cd[,sample(colnames(cd), 32, replace=F)]
  marrow <- CreateSeuratObject(counts = raw[[i]])
  marrow <- NormalizeData(marrow)
  marrow <- FindVariableFeatures(marrow, selection.method = "vst")
  marrow <- ScaleData(marrow, features = rownames(marrow))
  marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  re.seurat[[i]]<-(marrow[[]][,6:7])
  colnames(re.seurat[[i]])<-c("pre","true")
}
cal.seurat<-do.call(rbind,cal(re.seurat))
cal.seurat.ave<-round(colMeans(cal.seurat),3)
#hela-2002
#normFishCts <- GAPDHnorm ( cd, gapdhThreshold = 0 )
setwd("../../")
load("mu_genes.RData")
obscycle<-lapply(genes,function(x){
  x<-intersect(x,rownames(cd))
})

getCorrelatedGenes <- function ( obsPhaseGenes , normCts ) {
  aves <- lapply ( obsPhaseGenes , function ( x ) { colMeans ( normCts [ x , ] ) } )
  cors <- lapply ( names ( obsPhaseGenes ), function ( phase ) {
    pGenes <- obsPhaseGenes [[ phase ]]
    pCors <- sapply ( pGenes , function ( x ) {
      cor ( normCts [ x , ] , aves [[ phase ]] )
    } )
    pCors [ ! is.finite ( pCors ) ] <- 0
    pCors
  } )
  names ( cors ) <- names ( obsPhaseGenes)
  lapply ( cors , function ( x ) { names ( x ) [ x >= 0.3 ] } )
}
raw<-list()
re.whi<-list()
#90%: 35 *0.9=31.5
for(i in 1:20){
  raw[[i]]<-as.matrix(cd[,sample(colnames(cd), 32, replace=F)])
  dsPhaseGenes <- getCorrelatedGenes ( obscycle , raw[[i]])
###### Assign sample score
  scoreMatrix <- assignSampleScore ( dsPhaseGenes , raw[[i]] )
###### Normalize scores
  normScores <- getNormalizedScores ( scoreMatrix )
###### Score cycle similarity
  refCors <- assignRefCors ( normScores )
  assignedPhase <- assignPhase ( refCors )
  re.whi[[i]]<-data.frame(assignedPhase,substr(names(assignedPhase),1,2))
  colnames(re.whi[[i]])<-c("pre","true")
}
#*S:S;*G2,G2/M,M/G1:G2M;*G1/S:G1
cal1<-function(data){
  lapply(data,function(x){
    x2<-table(x)
    x1<-as.matrix(rbind(x2[1,],colSums(x2[2:4,]),x2[5,]))
    G1.TPR=x1[1,1]/sum(x1[,1])
    G1.FPR=sum(x1[1,-1])/(sum(x1[,-1])+sum(x1[1,-1]))
    G1.F1=2*x1[1,1]/(sum(x1[,1])+sum(x1[1,]))
    #S
    M.TPR=x1[2,2]/sum(x1[,2])
    M.FPR=sum(x1[2,-2])/(sum(x1[,-2])+sum(x1[2,-2]))
    M.F1=2*x1[2,2]/(sum(x1[,2])+sum(x1[2,]))
    #S
    S.TPR=x1[3,3]/sum(x1[,3])
    S.FPR=sum(x1[3,-3])/(sum(x1[,-3])+sum(x1[3,-3]))
    # S.prescion=x1[3,3]/sum(x1[3,])
    S.F1=2*x1[3,3]/(sum(x1[,3])+sum(x1[3,]))
    #S.F11=2 * S.TPR * S.prescion/(S.TPR + S.prescion)
    x11<-data.frame(G1.TPR,G1.FPR,G1.F1,S.TPR,S.FPR,S.F1,M.TPR,M.FPR,M.F1)
    return(x11)
  })}
cal.whi<-do.call(rbind,cal1(re.whi))
cal.whi.ave<-round(colMeans(cal.whi),3)
re.recat<-readRDS("./auc/re.recat.rds")
cal.recat<-readRDS("./auc/cal.recat.rds")
cal.recat.ave<-round(colMeans(cal.recat),3)
mesc35.sample<-list("re.cyclone.merge"=re.cyclone.merge,"re.cyclone.raw"=re.cyclone.raw,
                    "re.recat"=re.recat,"re.seurat"=re.seurat,"re.whi"=re.whi)
mesc35.cal<-list("cal.cyclone.merge"=cal.cyclone.merge,"cal.cyclone.raw"=cal.cyclone.raw,
                 "cal.recat"=cal.recat,"cal.seurat"=cal.seurat,"cal.whi"=cal.whi)
saveRDS(mesc35.cal,file="mesc35.cal.rds")
saveRDS(mesc35.sample,file="mesc35.sample.rds")
mesc35.cal.ave<-rbind(cal.cyclone.merge.ave,cal.cyclone.raw.ave,cal.recat.ave,cal.seurat.ave,cal.whi.ave)
write.table(mesc35.cal.ave,file="./auc/mesc35_32.f1_score.txt",sep="\t",quote=F,row.names = T,col.names = T)

#plot heatmap
library(viridis)
library(ggplot2)
library(hrbrthemes)
library(plotly)
rawfile=dir()
file<-list()
for(i in 1:length(rawfile)){
 file[[i]]<-read.table(file=rawfile[i],sep="\t",header=T,row.names = 1)
}
names(file) <-c("hESC-SMART","mESC-SMART","mESC-Quartz")
data<-do.call(cbind,file)
dat<-data[,c(3,6,9,12,15,18,21,24,27)]
library(reshape2)
dat$names<-rownames(dat)
df<-melt(dat)
library(tidyr)
dfc<-separate(df,variable ,c("sample","satge","F1"),sep="\\.")
table(dfc[,1:2])
dfc$X<-rep(1:3,each=5)
library(viridis)
library(ggplot2)
library(hrbrthemes)
library(plotly)
dfc<-within(dfc,names<-factor(names,levels=rev(unique(dfc$names))))
ggplot(dfc, aes(X, names, fill= value,text=satge)) + 
  geom_tile() + facet_wrap( ~sample) +
  scale_fill_viridis(discrete=FALSE) +
  theme_ipsum()


ggplot(dfc, aes(X, names, fill= value,text=satge)) + 
  geom_tile() + facet_wrap( ~sample) +
  scale_fill_viridis(discrete=FALSE) +
  theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text = element_text(size = 12), 
          plot.title = element_text(size = 16, face = 'bold'), 
          axis.title = element_text(face = 'bold'))