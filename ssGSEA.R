

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")


#???ð?
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="merge.txt"               #?????????ļ?
clusterFile="GRGcluster.txt"      #???͵Ľ????ļ?
gmtFile="immune.gmt"              #???߻??????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\18.ssGSEA")     #???ù???Ŀ¼

#??ȡ?????????ļ?,?????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#??ȡ???????ļ?
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA????
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#??ssGSEA???ֽ??н???
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#????ssGSEA???ֽ???
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#??ȡ???͵Ľ????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#???ݺϲ?(?????͵Ľ?????????ϸ???Ĵ??ֽ??кϲ?)
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#??????ת????ggplot2?????ļ?
data=melt(scoreCluster, id.vars=c("GRGcluster"))
colnames(data)=c("GRGcluster", "Immune", "Fraction")

#????????ͼ
bioCol=c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
         "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
bioCol=bioCol[1:length(levels(factor(data[,"GRGcluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="GRGcluster",
     xlab="",
     ylab="Immune infiltration",
     legend.title="GRGcluster",
     palette=bioCol)
p=p+rotate_x_text(50)

#????ͼ??
pdf(file="boxplot.pdf", width=8, height=6)
p+stat_compare_means(aes(group=ARGcluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()



