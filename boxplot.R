

#install.packages("ggplot2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(ggplot2)
library(ggpubr)

pFilter=0.001                      #pvalue????????
riskFile="risk.all.txt"            #?????ļ?
drugFile="DrugPredictions.csv"     #ҩ?????????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\32.oncoPredict")     #???ù???Ŀ¼

#??ȡ?????????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#??ȡҩ?????????ļ?
senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))

#???ݺϲ?
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
rt=cbind(risk, senstivity)

#???ñȽ???
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#??ҩ??????ѭ??, ????????ͼ
for(drug in colnames(rt)[2:ncol(rt)]){
	rt1=rt[,c(drug, "risk")]
	colnames(rt1)=c("Drug", "Risk")
	rt1=na.omit(rt1)
	rt1$Drug=log2(rt1$Drug+1)
	#????????
	test=wilcox.test(Drug ~ Risk, data=rt1)
	diffPvalue=test$p.value
	#????????????ҩ??????????ͼ
	if(diffPvalue<pFilter){
		boxplot=ggboxplot(rt1, x="Risk", y="Drug", fill="Risk",
					      xlab="Risk",
					      ylab=paste0(drug, " senstivity"),
					      legend.title="Risk",
					      palette=c("#00A087B2",  "#E64B35B2")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		#????ͼ??
		pdf(file=paste0("drugSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}



