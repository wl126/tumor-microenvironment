#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library("limma")
library("ggpubr")

inputFile=""        
gene=""                  
setwd("")    


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
normal=data[,group!=0]
tumor=data[,group==0]


normal=rbind(normal,gene=normal[gene,])
normal=as.matrix(t(normal[c("gene",gene),]))
rownames(normal)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(normal))
normal=avereps(normal)


tumor=rbind(tumor,gene=tumor[gene,])
tumor=as.matrix(t(tumor[c("gene",gene),]))
rownames(tumor)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(tumor))
tumor=avereps(tumor)


samSample=intersect(row.names(normal),row.names(tumor))
normal=normal[samSample,]
tumor=tumor[samSample,]
data=cbind(normal,tumor)
data=as.data.frame(data[,c(1,3)])
colnames(data)=c("Normal","Tumor")


pdf(file="pairDiff.pdf",width=5.5,height=5)
ggpaired(data, cond1 = "Normal", cond2 = "Tumor",fill = "condition", palette = "jco",
    xlab="",ylab = paste0(gene," expression"))+
    stat_compare_means(paired = TRUE, label = "p.format", label.x = 1.35)
dev.off()


