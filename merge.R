

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sva")


#???ð?
library(limma)
library(sva)

setwd("C:\\Users\\17860\\Desktop\\OV\\06.merge")     #???ù???Ŀ¼
files=c("symbol.txt", "GSE26193.txt")      #?????ļ?????

#???????ļ?????ѭ??,??ȡ??????????????
geneList=list()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile, header=T, sep="\t",check.names=F)
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    header[1]=gsub("symbol", "TCGA", header[1])
    geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

#???ݺϲ?
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    header[1]=gsub("symbol", "TCGA", header[1])
    #??ȡ?????ļ????????????ļ?????????
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))
    #??????TCGA??????,ɾ????????Ʒ
    if(header[1] == "TCGA"){
		group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
		group=sapply(strsplit(group,""), "[", 1)
		rt=rt[,group==0]
		rt=t(rt)
		row.names(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(rt))
		rt=avereps(rt)
		rt=t(rt)
    }
    #?ж??????Ƿ?ȡ??log2,????û??ȡlog2,????ֵ?Զ?ȡlog2
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    if(header[1] != "TCGA"){
    	rt=normalizeBetweenArrays(rt)
    }
    #???ݺϲ?
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#?????ݽ??????ν????????????????ı???????
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)




