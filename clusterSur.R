

#install.packages("survival")
#install.packages("survminer")


#???ð?
library(survival)
library(survminer)

clusterFile="GRGcluster.txt"     #???͵Ľ????ļ?
cliFile="time.txt"               #?????????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\14.clusterSur")      #???ù???Ŀ¼

#??ȡ???͵Ľ????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
#??ȡ?????????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#???ݺϲ?
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#????????????,?۲첻ͬ????֮?䲡?˵??????Ƿ????в???
length=length(levels(factor(rt$GRGcluster)))
diff=survdiff(Surv(futime, fustat) ~ GRGcluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ GRGcluster, data = rt)
#print(surv_median(fit))

#????????????
bioCol=c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
         "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="GRGcluster",
		           legend.labs=levels(factor(rt[,"GRGcluster"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           break.time.by = 2,
		           palette = bioCol,
		           surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.25)

#????ͼ??
pdf(file="survival.pdf", width=6.5, height=5.25, onefile=FALSE)
print(surPlot)
dev.off()




