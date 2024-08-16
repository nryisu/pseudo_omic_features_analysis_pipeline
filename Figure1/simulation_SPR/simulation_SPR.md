1. calculated the different SPR reads number in raw data
   The mESC has total reads:  188502781, and the SPR supposed to be 55%.
   The MEF has total reads:  380641607, and the SPR supposed to be 20%.
   according to the equation,if SPR is 52%:
   380641607/188502781=2.019289;
   x*0.55 * 188502781 +y*0.2 *380641607 =0.52 * 188502781 (1);
   x*188502781+y*380641607=188502781 (2); (y<0.2101715,x<0.4243969)
   (1) simplifying：0.55x+ 0.4038578y=0.52
   (2) simplifying：x =1-2.019289y; 0.55-0.7067512y=0.52
   y=0.04244775,x= 0.9142857
   the SPR is 50%: y=(0.55-0.50)/0.7067512= 0.07074625; x=1-2.019289 *0.07074625=0.8571429
   the SPR is 48%: y=(0.55-0.48)/0.7067512= 0.0990447; x=1-2.019289 *0.0990447=0.8000001
   the SPR is 46%: y=(0.55-0.46)/0.7067512= 0.1273433; x=1-2.019289 *0.1273433=0.7428571
   the SPR is 44%: y=(0.55-0.44)/0.7067512= 0.1556418; x=1-2.019289 *0.1556418=0.6857142
   the SPR is 42%: y=(0.55-0.42)/0.7067512= 0.1839403; x=1-2.019289 *0.1839403=0.6285714
   the SPR is 40%: y=(0.55-0.40)/0.7067512=  0.2122388; x=1-2.019289 * 0.2122388=0.5714285
   the SPR is 30%: y=(0.55-0.30)/0.7067512= 0.3537313; x=1-2.019289 *0.3537313=0.2857143
   the SPR is 20%: y= 0.4952238;x=0
2. split reads to 1G size for running fast, then randomly sample reads from data. The random function is referred to
   https://www.biostars.org/p/6544/  and https://bioinformatics.stackexchange.com/questions/807/random-access-on-a-fastq-file 

​	For example: if the SPR is 50%:

```shell
for X in mESCs:
python sample_reads.py 0.8571429 ESC_1.fastq.gz ESC_2.fastq.gz  ESC_1_random.fq  ESC_2_random.fq 
for y in MEFs:
python sample_reads.py 0.07074625 MEF_1.fastq.gz  MEF_2.fastq.gz  MEF_1_random.fq  MEF_2_random.fq 
```

3.  After randomly sample reads, then merge the sample ESCs and MEFs reads, process them according to the WGS method for calling CNV segments, and calculated the enrich score for visualization（Figure 1d）.

```R
#For total number CNV segment with replication timing
library("ggpubr")
SPR<-c("55","52","50","48","46","44","42","40","30")
number<-c(437,356,307,293,289,264,232,240,190)
df<-data.frame(SPR,number)
df$SPR<-as.numeric(df$SPR)
par(mar=c(0.5,0.5,0.5,0.5))
ggscatter(df, x = "SPR", y = "number",add = "reg.line",
   add.params = list(color = "blue", fill = "lightgray"),
   conf.int = TRUE) + stat_cor(method = "pearson") + ylim(150,450) +xlim(15,65)
#For Gain ES
library("ggpubr")
SPR<-c("55","52","50","48","46","44","42","40","30","20")
Gain_ES<-c(6.706,6.481,6.415,3.717,4.46,5.0623,2.728,1.899,1.736,0.597)
df<-data.frame(SPR,Gain_ES)
df$SPR<-as.numeric(df$SPR)
par(mar=c(0.5,0.5,0.5,0.5))
ggscatter(df, x = "SPR", y = "Gain_ES",add = "reg.line",
   add.params = list(color = "blue", fill = "lightgray"),
   conf.int = TRUE) + stat_cor(method = "pearson") + ylim(0,8) +xlim(15,65)
#For Loss ES
library("ggpubr")
SPR<-c("55","52","50","48","46","44","42","40","30","20")
Loss_ES<-c(3.765,3.03,2.7,2.304,2.289,2.117,1.739,1.555,1.334,0.819)
df<-data.frame(SPR,Loss_ES)
df$SPR<-as.numeric(df$SPR)
par(mar=c(0.5,0.5,0.5,0.5))
ggscatter(df, x = "SPR", y = "Loss_ES",add = "reg.line",
   add.params = list(color = "blue", fill = "lightgray"),
   conf.int = TRUE) + stat_cor(method = "pearson") + ylim(0,4) +xlim(15,65)
```

