library(SingleCellExperiment)
library(scater)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(scran)
library(Rtsne)


Args <- commandArgs()
mesc_npc<- read.table(Args[6], header = TRUE, skip = 1, row.names = 1)
colnames(mesc_npc) <- gsub(".bam", "", colnames(mesc_npc), fixed = T)
counts<- as.matrix(mesc_npc[ ,c(-1:-5)])
table(substr(colnames(mesc_npc),1,3))#mesc:183,npc:54
genes.length<-mesc_npc$Len
genes<-read.delim("genes_id_names",sep="\t",header=F)
colnames(genes)<-c("ids","names")
geneslist <- list(names=rownames(counts),ids=genes[match(rownames(counts), genes$names),]$ids)
rownames(counts) <- make.unique(as.character(geneslist$ids))
sce <- SingleCellExperiment(list(counts=counts))
rowData(sce)$GeneLength <- genes.length
sce$sg<-factor(substr(colnames(counts),1,3),levels=c("mes","npc"))
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- make.unique(as.character(geneslist$names))
# replace missing symbols with the Ensembl identifier
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
    column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
png(paste(Args[6],"_qc.png",sep=""))
par(mfrow=c(1,2))
plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab="Number of expressed genes",
    ylab="Library size (millions)")
plot(sce$total_features_by_counts, sce$pct_counts_Mt, xlab="Number of expressed genes",
    ylab="Mitochondrial proportion (%)")
dev.off()
#Outliers are defined based on the median absolute deviation (MADs) from the median value of each metric across all cells. 
#We remove cells with log-library sizes that are more than 3 MADs below the median log-library size. 
#A log-transformation improves resolution at small values, especially when the MAD of the raw values is comparable to or greater than the median.
# We also remove cells where the log-transformed number of expressed genes is 3 MADs below the median value.
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce$sg)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce$sg)
keep <- !(libsize.drop | feature.drop )
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=sum(keep))
sce$PassQC <- keep
saveRDS(sce, file=paste(Args[6],"_preQC.rds",sep=""))
sce <- sce[,keep]
dim(sce)
set.seed(100)
counts<-counts(sce)
load("../annotation/mm.merge.pair.RData")
########################classify cell cycle in mesc 
umi <- SingleCellExperiment(assays = list(counts = as.matrix(counts[,1:111])),colData=colnames(counts[,1:111]))
assignments <- cyclone(umi, merge.pair, gene.names=rownames(counts(umi)))
table(assignments$phases) 
png(paste(Args[6],"_rawcycle.png",sep=""))
plot(assignments$score$G1, assignments$score$G2M, 
    xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()
G1 <- counts[,1:111][,assignments$phases=="G1"]
S <- counts[,1:111][,assignments$phases=="S"]
G2M <- counts[,1:111][,assignments$phases=="G2M"]
colnames(G1)<-paste("G",colnames(G1),sep="_")
colnames(S)<-paste("S",colnames(S),sep="_")
colnames(G2M)<-paste("M",colnames(G2M),sep="_")
mesc<-cbind(G1,S,G2M)
dim(mesc)
#plot tsne #log2(TPM+1)
cell <-substr(colnames(mesc),0,1)
sg<- factor(cell,levels=c("G","S","M"))
table(sg)
merge<-c(as.vector(merge.pair$G1$first),as.vector(merge.pair$G1$second),
		as.vector(merge.pair$S$first),as.vector(merge.pair$S$second),
		as.vector(merge.pair$G2M$first),as.vector(merge.pair$G2M$second))
merge <- merge[!duplicated(merge)]
length(merge) #890
mesc <- mesc[rowSums(mesc)>0,]
# remove genes that are not seen in a sufficient number of cells
mesc<-mesc[rowSums(mesc>0)>=5, ] #0.05
merge.o<-intersect(merge,rownames(mesc))
lib_size = colSums(mesc)
norm <- t(t(mesc)/lib_size * median(lib_size))
#log transform
mat <- log2(as.matrix(norm)+1)
mat<-mat[merge.o,]
cell.labels <-substr(colnames(mesc),0,1)
sg <- factor(cell.labels,levels=c("G","S","M"))
set.seed(2222)
png(paste(Args[6],"_rawcycle.mesc.png",sep=""))
joint_tsne <- Rtsne(t(mat),initial_dims=2, theta=0.5,perplexity=9,
	check_duplicates=FALSE, max_iter=200, stop_lying_iter=50, mom_switch_iter=50,dims=2)
sgcol=c("#FF3366","#00CC66","#3399FF")[sg]
plot(joint_tsne$Y,col=sgcol, pch=16, main='tSNE')
legend("topright", legend=c("G1","S","G2M"),col=c("#FF3366","#00CC66","#3399FF"), pch=16, cex=0.8)
dev.off()
#########################classify cell cycle in mnpc
umi <- SingleCellExperiment(assays = list(counts = as.matrix(counts[,112:153])),colData=colnames(counts[,112:153]))
assignments <- cyclone(umi, merge.pair, gene.names=rownames(counts(umi)))
G1 <- counts[,112:153][,assignments$phases=="G1"]
S <- counts[,112:153][,assignments$phases=="S"]
G2M <- counts[,112:153][,assignments$phases=="G2M"]
colnames(G1)<-paste("G",colnames(G1),sep="_")
colnames(S)<-paste("S",colnames(S),sep="_")
colnames(G2M)<-paste("M",colnames(G2M),sep="_")
npc<-cbind(G1,S,G2M)
npc <- npc[rowSums(npc)>0,]
# remove genes that are not seen in a sufficient number of cells
npc<-npc[rowSums(npc>0)>=3, ] #0.05
merge.n<-intersect(merge,rownames(npc))
lib_size = colSums(npc)
norm <- t(t(npc)/lib_size * median(lib_size))
#log transform
mat.n <- log2(as.matrix(norm)+1)
mat.n<-mat.n[merge.n,]
cell.labels <-substr(colnames(npc),0,1)
sg.n <- factor(cell.labels,levels=c("G","S","M"))
set.seed(222)
png(paste(Args[6],"_rawcycle.mnpc.png",sep=""))
joint_tsne.n <- Rtsne(t(mat.n),initial_dims=2, theta=0.5,perplexity=8,
	check_duplicates=FALSE, max_iter=200, stop_lying_iter=50, mom_switch_iter=50,dims=2)
sgcol.n=c("#FF3366","#00CC66","#3399FF")[sg.n]
plot(joint_tsne.n$Y,col=sgcol.n, pch=16, main='tSNE')
legend("topleft", legend=c("G1","S","G2M"),col=c("#FF3366","#00CC66","#3399FF"), pch=16, cex=0.8)
dev.off()
#############obtain the raw count matrix 
counts<- as.matrix(mesc_npc[ ,c(-1:-5)])
G.mesc<-colnames(mesc[,1:4])
G.mesc<-substr(G.mesc,3,nchar(G.mesc)+2)
G.mesc<-counts[,G.mesc]
S.mesc<-colnames(mesc[,5:103])
S.mesc<-substr(S.mesc,3,nchar(S.mesc)+2)
S.mesc<-counts[,S.mesc]
M.mesc<-colnames(mesc[,104:111])
M.mesc<-substr(M.mesc,3,nchar(M.mesc)+2)
M.mesc<-counts[,M.mesc]
colnames(G.mesc)<-paste("G",colnames(G.mesc),sep="_")
colnames(S.mesc)<-paste("S",colnames(S.mesc),sep="_")
colnames(M.mesc)<-paste("M",colnames(M.mesc),sep="_")
#################
G.npc<-colnames(npc[,1:29])
G.npc<-substr(G.npc,3,nchar(G.npc)+2)
G.npc<-counts[,G.npc]
S.npc<-colnames(npc[,30:35])
S.npc<-substr(S.npc,3,nchar(S.npc)+2)
S.npc<-counts[,S.npc]
M.npc<-colnames(npc[,36:42])
M.npc<-substr(M.npc,3,nchar(M.npc)+2)
M.npc<-counts[,M.npc]
colnames(G.npc)<-paste("G",colnames(G.npc),sep="_")
colnames(S.npc)<-paste("S",colnames(S.npc),sep="_")
colnames(M.npc)<-paste("M",colnames(M.npc),sep="_")
count<-cbind(G.mesc,S.mesc,M.mesc,G.npc,S.npc,M.npc)
saveRDS(count,file=paste(Args[6],"_cycle.rds",sep=""))



