library(SingleCellExperiment)
library(scater)
library(edgeR)
library(MAST)
library(ggplot2)
library(RColorBrewer)
library(NMF)
library(data.table)

Args <- commandArgs()
count<-readRDS(Args[6])
# filter out low-gene cells (often empty wells):at least 1000 genes.
counts<-count[, colSums(count>0)>1000]
cell.labels <-substr(colnames(counts),3,5)
grp <- factor(cell.labels,levels=c("mes","npc"))
names(grp)<-colnames(counts)
cdr <- scale(colMeans(counts > 0))
dge <- DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
cpms <- cpm(dge)
sca <- FromMatrix(exprsArray = log2(cpms + 1), cData = data.frame(wellKey = names(grp),grp = grp, cdr = cdr))
zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
summaryCond <-summary(zlmdata,logFC = TRUE,doLRT=TRUE) 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='grpnpc' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='grpnpc' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
FCTHRESHOLD <- log2(2)
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD],as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
#########################G vs G#####
counts.g<-counts[,c(1:4,112:140)]
cell.labels <-substr(colnames(counts.g),1,5)
grp<- factor(cell.labels,levels=c("G_mes","G_npc"))
names(grp)<-colnames(counts.g)
cdr <- scale(colMeans(counts.g > 0))
dge <- DGEList(counts = counts.g)
dge <- edgeR::calcNormFactors(dge)
cpms <- cpm(dge)
sca <- FromMatrix(exprsArray = log2(cpms + 1), cData = data.frame(wellKey = names(grp),grp = grp, cdr = cdr))

zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
summaryCond <-summary(zlmdata,logFC = TRUE,doLRT=TRUE) 
summaryDt <- summaryCond$datatable
fcHurdle.g <- merge(summaryDt[contrast=='grpG_npc' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='grpG_npc' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle.g[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig.g <- merge(fcHurdle.g[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig.g, fdr)
#######S vs S
counts.s<-counts[,c(5:103,141:146)]
cell.labels <-substr(colnames(counts.s),1,5)
grp<- factor(cell.labels,levels=c("S_mes","S_npc"))
names(grp)<-colnames(counts.s)
cdr <- scale(colMeans(counts.s > 0))
dge <- DGEList(counts = counts.s)
dge <- edgeR::calcNormFactors(dge)
cpms <- cpm(dge)
sca <- FromMatrix(exprsArray = log2(cpms + 1), cData = data.frame(wellKey = names(grp),grp = grp, cdr = cdr))
zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
summaryCond <-summary(zlmdata,logFC = TRUE,doLRT=TRUE) 
summaryDt <- summaryCond$datatable
fcHurdle.s <- merge(summaryDt[contrast=='grpS_npc' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='grpS_npc' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle.s[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
FCTHRESHOLD <- log2(2)
fcHurdleSig.s <- merge(fcHurdle.s[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig.s, fdr)
###########M vs M
counts.m<-counts[,c(104:111,147:153)]
cell.labels <-substr(colnames(counts.m),1,5)
grp<- factor(cell.labels,levels=c("M_mes","M_npc"))
names(grp)<-colnames(counts.m)
cdr <- scale(colMeans(counts.m > 0))
dge <- DGEList(counts = counts.m)
dge <- edgeR::calcNormFactors(dge)
cpms <- cpm(dge)
sca <- FromMatrix(exprsArray = log2(cpms + 1), cData = data.frame(wellKey = names(grp),grp = grp, cdr = cdr))
zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
summaryCond <-summary(zlmdata,logFC = TRUE,doLRT=TRUE) 
summaryDt <- summaryCond$datatable
fcHurdle.m <- merge(summaryDt[contrast=='grpM_npc' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='grpM_npc' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle.m[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
FCTHRESHOLD <- log2(2)
fcHurdleSig.m <- merge(fcHurdle.m[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig.m, fdr)
#save(fcHurdle,file="fcHurdle.RData")
#save(fcHurdle.g,file="fcHurdle.g.RData")
#save(fcHurdle.s,file="fcHurdle.s.RData")
#save(fcHurdle.m,file="fcHurdle.m.RData")
whole<-data.frame(fcHurdleSig$coef,fcHurdleSig$fdr)
rownames(whole)<-fcHurdleSig$primerid
G<-data.frame(fcHurdleSig.g$coef,fcHurdleSig.g$fdr)
rownames(G)<-fcHurdleSig.g$primerid
S<-data.frame(fcHurdleSig.s$coef,fcHurdleSig.s$fdr)
rownames(S)<-fcHurdleSig.s$primerid
M<-data.frame(fcHurdleSig.m$coef,fcHurdleSig.m$fdr)
rownames(M)<-fcHurdleSig.m$primerid
library(gplots)
data<-list(G,S,M,whole)
names(data)<-c("G1","S","G2M","whole")
save(data,file="data.RData")