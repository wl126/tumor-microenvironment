#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")

library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

immuneFile=""         
expFile=""                     
gene=""                              
pFilter=0.05                              
setwd("")     


immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])


rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]


data=rbind(data,gene=data[gene,])
exp=t(data[c("gene",gene),])
sameSample=intersect(row.names(immune),row.names(exp))
immune1=immune[sameSample,]
exp1=exp[sameSample,]


outTab=data.frame()
x=as.numeric(exp1[,1])

for(j in colnames(immune1)[1:22]){
	y=as.numeric(immune1[,j])
	if(sd(y)>0.001){
		df1=as.data.frame(cbind(x,y))
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pValue=corT$p.value
		p1=ggplot(df1, aes(x, y)) + 
			ylab(j)+xlab(gene)+
			geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
		if(pValue<pFilter){
			pdf(file=paste0(j,".pdf"),width=5,height=5)
			print(p1)
			dev.off()
			outTab=rbind(outTab,cbind(Cell=j,pValue))
		}
	}
}
write.table(outTab,file="immuneCor.result.txt",sep="\t",row.names=F,quote=F)

