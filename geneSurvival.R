#install.packages("survival")
#install.packages("survminer")


library(survival)
library(survminer)

inputFile=""     
gene=""                 
setwd("")    


rt=read.table(inputFile,header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365        

a=ifelse(rt[,gene]<=median(rt[,gene]),"Low","High")

diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)

fit=survfit(Surv(futime, fustat) ~ a, data = rt)

if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
titleName=gene
surPlot=ggsurvplot(fit, 
				data=rt,
				conf.int=TRUE,
				pval=pValue,
				pval.size=6,
				risk.table=T,
				#ncensor.plot = TRUE,
				legend.labs=c("high","low"),
				legend.title=titleName,
				xlab="Time(years)",
				break.time.by = 1,
				risk.table.title="",
				palette=c("red", "blue"),
				risk.table.height=.25)          
pdf(file=paste0("sur.",gene,".pdf"), width = 6.5, height = 5.5,onefile = FALSE)
print(surPlot)
dev.off()

