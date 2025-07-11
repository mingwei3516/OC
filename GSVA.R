


#???ð?
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

expFile="merge.txt"                 #?????????ļ?
clusterFile="GRGcluster.txt"        #???͵Ľ????ļ?
gmtFile="c2.cp.kegg.symbols.gmt"    #???????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\19.GSVA")      #???ù???Ŀ¼

#??ȡ?????????ļ?,?????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#GSVA????
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#??ȡ???͵Ľ???
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#???ݺϲ?(?????ͺ?GSVA?Ľ??????кϲ?)
gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
gsvaCluster=cbind(gsvaCluster, Project)

#???ñȽ???
adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$GRGcluster)
comp=combn(levels(factor(allType)), 2)

#?ԱȽ???????ѭ??, ?۲???Щͨ·?ڲ?ͬ?ķ???֮???Ǿ??в?????
for(i in 1:ncol(comp)){
	#??Ʒ????
	treat=gsvaCluster[gsvaCluster$GRGcluster==comp[2,i],]
	con=gsvaCluster[gsvaCluster$GRGcluster==comp[1,i],]
	data=rbind(con, treat)
	#??ͨ·???в???????
	Type=as.vector(data$GRGcluster)
	ann=data[,c(ncol(data), (ncol(data)-1))]
	data=t(data[,-c((ncol(data)-1), ncol(data))])
	design=model.matrix(~0+factor(Type))
	colnames(design)=levels(factor(Type))
	fit=lmFit(data, design)
	contrast=paste0(comp[2,i], "-", comp[1,i])
	cont.matrix=makeContrasts(contrast, levels=design)
	fit2=contrasts.fit(fit, cont.matrix)
	fit2=eBayes(fit2)
	
	#????????ͨ·?Ĳ???????
	allDiff=topTable(fit2,adjust='fdr',number=200000)
	allDiffOut=rbind(id=colnames(allDiff),allDiff)
	write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
	
	#??????????????ͨ·
	diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
	diffSigOut=rbind(id=colnames(diffSig),diffSig)
	write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
	
	#??????ͼע?͵???ɫ
	bioCol=c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
	         "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
	ann_colors=list()
	m6aCluCol=bioCol[1:length(levels(factor(allType)))]
	names(m6aCluCol)=levels(factor(allType))
	ann_colors[["GRGcluster"]]=m6aCluCol[c(comp[1,i], comp[2,i])]

	#???Ʋ???ͨ·??ͼ
	termNum=20     #??????ʾͨ·????Ŀ
	diffTermName=as.vector(rownames(diffSig))
	diffLength=length(diffTermName)
	if(diffLength<termNum){termNum=diffLength}
	hmGene=diffTermName[1:termNum]
	hmExp=data[hmGene,]
	pdf(file=paste0(contrast,".heatmap.pdf"), width=10, height=6)
	pheatmap(hmExp, 
	         annotation=ann,
	         annotation_colors = ann_colors,
	         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
	         cluster_cols =F,
	         show_colnames = F,
	         gaps_col=as.vector(cumsum(table(Type))),
	         scale="row",
	         fontsize = 8,
	         fontsize_row=6,
	         fontsize_col=8)
	dev.off()
}



