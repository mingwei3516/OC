

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(reshape2)
library(ggpubr)

expFile="uniSigExp.txt"      #?????????ļ?
cluFile="GRGcluster.txt"     #???͵Ľ????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\16.clusterdiff")       #???ù???Ŀ¼

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#??ȡ???͵Ľ????ļ?
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
sameSample=intersect(row.names(data), row.names(cluster))
expClu=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])

#??ȡ?????????Ļ???
sigGene=c()
for(i in colnames(expClu)[1:(ncol(expClu)-1)]){
	if(sd(expClu[,i])<0.001){next}
	if(length(levels(factor(expClu[,"GRGcluster"])))>2){
		test=kruskal.test(expClu[,i] ~ expClu[,"GRGcluster"])
	}else{
		test=wilcox.test(expClu[,i] ~ expClu[,"GRGcluster"])
	}
	pvalue=test$p.value
	if(pvalue<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "GRGcluster")
expClu=expClu[,sigGene]

#??????ת????ggplot2?????ļ?
data=melt(expClu, id.vars=c("GRGcluster"))
colnames(data)=c("GRGcluster", "Gene", "Expression")

#????ͼ????ɫ
bioCol=c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
         "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
bioCol=bioCol[1:length(levels(factor(data[,"ARGcluster"])))]

#????????ͼ
p=ggboxplot(data, x="Gene", y="Expression", color="GRGcluster",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="GRGcluster",
	     palette = bioCol,
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=GRGcluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#????????ͼ
pdf(file="boxplot.pdf", width=9, height=6)
print(p1)
dev.off()




