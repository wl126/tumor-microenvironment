#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("beeswarm")

library(limma)
library(beeswarm)

inputFile=""       
gene=""                 
yMin=0                       
yMax=165                     
setwd("C:\\Users\\lexb4\\Desktop\\TMEimmune\\24.geneDiff")      


rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
tumor=data[,group==0]
tumorNum=ncol(tumor)
normalNum=ncol(data)-tumorNum


Type=c(rep(1,normalNum),rep(2,tumorNum))
rt=rbind(expression=data[gene,],Type=Type)
rt=as.matrix(t(rt))


wilcoxTest=wilcox.test(expression ~ Type, data=rt)
pvalue=wilcoxTest$p.value
if(pvalue<0.001){
     pvalue="<0.001"
}else{
     pvalue=paste0("=",sprintf("%.03f",pvalue))
}


pdf(file="diff.pdf",width=6,height=5)
par(mar = c(4,6,3,3))
labels=c("Normal","Tumor")
boxplot(expression ~ Type, data = rt,names=labels,xlab="",
     ylab = paste(gene," expression",sep=""),
     cex.main=1.5, cex.lab=1.3, cex.axis=1.2,ylim=c(yMin,yMax),outline = FALSE)
beeswarm(expression ~ Type, data = rt, col = c("blue","red"),lwd=0.1,
     pch = 16, add = TRUE, corral="wrap")
ySeg=yMax*0.94
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.96);segments(2,ySeg, 2,ySeg*0.96)
text(1.5,ySeg*1.05,labels=paste("p",pvalue,sep=""),cex=1.2)
dev.off()

