

#install.packages("pheatmap")


library(pheatmap)       #???ð?
expFile="uniSigExp.txt"          #?????????ļ?
clusterFile="GRGcluster.txt"     #???͵Ľ????ļ?
cliFile="clinical.txt"           #?ٴ??????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\17.heatmap")     #???ù???Ŀ¼

#??ȡ?????????ļ?
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)
#??ȡ???͵Ľ????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ??????ͷ???????
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
rownames(expCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expCluster))
expCluster=cbind(expCluster, Project)

#?ϲ??ٴ?????
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

#??ȡ??ͼ??????
data=data[order(data$GRGcluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
data=t(data[,1:ncol(exp)])

#??????ͼע?͵???ɫ
bioCol=c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
         "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
ann_colors=list()
prgCluCol=bioCol[1:length(levels(factor(Type$GRGcluster)))]
names(prgCluCol)=levels(factor(Type$GRGcluster))
ann_colors[["GRGcluster"]]=prgCluCol

#??ͼ???ӻ?
pdf("heatmap.pdf", width=7.5, height=5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         show_colnames=F,
         scale="row",
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()




