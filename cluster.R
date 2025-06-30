

#install.packages("survival")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


#???ð?
library(limma)
library(survival)
library(ConsensusClusterPlus)

expFile="uniSigExp.txt"      #?????????ļ?
workDir="C:\\Users\\17860\\Desktop\\OV\\13.cluster"     #????Ŀ¼
setwd(workDir)       #???ù???Ŀ¼

#??ȡ?????ļ?
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#????Ʒ???з???
maxK=9      #??????kֵ(???????Խ???Ʒ?ֳɼ?????)
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123,
              plot="png")


#???????ͽ???
clusterNum=3      #?ֳɼ???????
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("ARGcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$ARGcluster))
cluster$ARGcluster=letter[match(cluster$ARGcluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="ARGcluster.txt", sep="\t", quote=F, col.names=F)



