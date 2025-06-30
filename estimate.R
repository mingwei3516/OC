

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#???ð?
library(limma)
library(estimate)
inputFile="merge.txt"      #?????????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\30.estimate")      #???ù???Ŀ¼

#??ȡ?ļ?,?????????ļ?????????
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#???????????ı????????ļ?
out=rbind(ID=colnames(data), data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#????estimate??
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

#??????΢?????Ĵ??ֽ???????, ????ÿ????Ʒ?Ĵ???
scores=read.table("estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)
write.table(out, file="TMEscores.txt", sep="\t", quote=F, col.names=F)



