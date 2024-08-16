library(rtracklayer)
library(DiffBind)
library(ChIPQC)
library(soGGi)
library(ChIPseeker)
library(DT)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(Rsubread)


peaks <- dir("./", pattern = "*filter.bed", full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
names(myPeaks) <- c("HCT116_G1", "HCT116_G2", "HCT116_mixture", "HCT116_S")
blklist <- import.bed("/data/nier/ref/hg19/ref/hg19_blacklist.JDB.bed")
myGRangesList<-GRangesList(myPeaks)   
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
consensusToCount <- reducedConsensus 
consensusToCount <- consensusToCount[!consensusToCount %over% blklist & !seqnames(consensusToCount) %in% "chrM"]
#head(as.data.frame(elementMetadata(consensusToCount)))
myCol <- brewer.pal(4, "Set2")
#plot overlap in samples
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(starts_with("HCT116")) %>% 
  vennDiagram(main = "Overlap for HCT116 open regions")
#plot PCA
myPlot <- as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(-consensusIDs) %>% 
  as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = gsub("_\\d", "", Samples)) %>% ggplot(aes(x = PC1, y = PC2, 
                                                           colour = Group)) + geom_point(size = 5)
myPlot
#Counting for diff
write.table(as.data.frame(consensusToCount),file="merge_peak.txt",quote=F,sep="\t")
#plot hc clust
df<- as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(-consensusIDs) %>% as.matrix 
hc <- hclust(dist(t(df)))
plot(hc)
df1<-df[,c(1,2,4)]
hc <- hclust(dist(t(df1)))
plot(hc)
#upstream 2kb
gene<-read.table("/data/nier/ref/hg19/gencode.v19.chr_patch_hapl_scaff.genes_with_promoters.sorted.bed",sep="\t",stringsAsFactors = F,header=F)
genes<-gene[,c(4,1,2,3,6)]
colnames(genes)<-c("GeneID","Chr","Start","end","Strand")
setwd("/data/nier/ATAC/HCT116/bam")
library(Rsamtools)
bamsToCount <- dir("./", full.names = TRUE, pattern = "*.\\.bam$")
indexBam(bamsToCount)
fcResults <- featureCounts(bamsToCount, annot.ext = genes, isPairedEnd = TRUE, 
                           countMultiMappingReads = FALSE, maxFragLength = 100)
myCounts <- fcResults$counts
colnames(myCounts)<-c("G1","G2M","mixture","S")
write.table(myCounts,file="merge_peak_gene_count.txt",quote=F,sep="\t")

# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
countdata <- read.table("/data/nier/ATAC/HCT116/RNA/bam/hct116_RNA_gencode_count.txt", header = TRUE, skip = 1, row.names = 1)
# Remove .bam + '..' from column identifiers
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("..", "", colnames(countdata), fixed = T)
dim(countdata)
# Remove length/char columns
countdata <- countdata[ ,c(5,6,7)]
# Make sure ID's are correct
#head(countdata)
data<- merge(myCounts,countdata,by="row.names",all=T)
dat<-data[,c(2:5,7,8)]
rownames(dat)<-data$Row.names
colnames(dat)<-c(colnames(dat)[1:4],"RNA_rep1","RNA_rep2")
library(DESeq2)
metaData <- data.frame("sample"=colnames(dat), 
                       row.names = colnames(dat),"group"=c(rep("ATAC",4),rep("RNA",2)),
                       "phase"=c("G1","G2M","mixture","S","bulk","bulk"))
dat[is.na(dat)]<-0
#deseq normalize
# DESeq2 can not handle them.
dat[dat>1e9] <- 1e9
#remove batch
library("sva")
batch <- c(rep(1, 4), rep(2, 2))
adjusted <- ComBat_seq(dat, batch=batch, group=NULL)
coldata <- matrix(NA, nrow=ncol(adjusted), ncol=1)
dds <- DESeqDataSetFromMatrix(countData = adjusted,
                              colData = coldata,
                              design= ~ 1)
dds <- estimateSizeFactors(dds) # run this or DESeq() first
# return the normalized tag counts
dat1<-counts(dds, normalized=T)
myPlot <- dat1[,1:4] %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = Samples) %>% ggplot(aes(x = PC1, y = PC2, colour = Group)) + geom_point(size = 5)

myPlot
dat2<-stats::cor(dat1)

##DEG genes
atacDDS <- DESeqDataSetFromMatrix(adjusted, metaData, ~group)
atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)
plotPCA(atac_Rlog, intgroup = "group", ntop = nrow(atac_Rlog))
ATACvsRNA<-results(atacDDS, c("group", "ATAC", "RNA"))
ATACvsRNA<-ATACvsRNA[order(abs(ATACvsRNA$log2FoldChange),decreasing = T),]
head(ATACvsRNA)
summary(abs(ATACvsRNA$log2FoldChange))
ATACvsRNA1<- subset(ATACvsRNA,log2FoldChange<0.5 & log2FoldChange>-0.5)
dim(ATACvsRNA1)

filter_genes<-rownames(ATACvsRNA1)
stats::cor(dat1[filter_genes,])
library(psych)
corPlot(dat1[filter_genes,], cex = 1.2)
##enrichment annalysis
library(clusterProfiler)
library(ChIPseeker)
filter_genes1<-as.data.frame(filter_genes)
filter_genes1$entrez <- mapIds(x = org.Hs.eg.db,keys = as.character(filter_genes1$filter_genes),
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals = "first")
filter_genes_entrez<- subset(filter_genes1, is.na(entrez) == FALSE)
dim(filter_genes_entrez)
go1 <- enrichGO(filter_genes_entrez$entrez, OrgDb = "org.Hs.eg.db", ont = "BP", maxGSSize = 5000)

filter_genes2<-as.data.frame(setdiff(rownames(ATACvsRNA),filter_genes))
colnames(filter_genes2)<-"filter_genes"
filter_genes2$entrez <- mapIds(x = org.Hs.eg.db,keys = as.character(filter_genes2$filter_genes),
                               column = "ENTREZID",
                               keytype = "SYMBOL",
                               multiVals = "first")
filter_genes_entrez1<- subset(filter_genes2, is.na(entrez) == FALSE)
dim(filter_genes_entrez1)
go2 <- enrichGO(filter_genes_entrez1$entrez, OrgDb = "org.Hs.eg.db", ont = "BP", maxGSSize = 5000)
go1_all<- go1[,c("ID", "Description", "p.adjust")]
go2_all<- go2[,c("ID", "Description", "p.adjust")]

head(go1, 10) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle1")
head(go2, 10) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% datatable(elementId = "goEle1")

gcSample=list("high"= as.character(filter_genes_entrez$entrez),"low"= as.character(filter_genes_entrez1$entrez))
str(gcSample)
ck <- compareCluster(geneCluster = gcSample, fun = "enrichGO",OrgDb = "org.Hs.eg.db", ont = "BP")
pdf("compare_high_low.pdf",width=12.5,height = 4.5)
dotplot(ck)
dev.off()
head(ck@compareClusterResult$Description)