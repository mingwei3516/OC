

library(limma)
library(survival)
library(survminer)

expFile="diffGRGexp.txt"      #?????????ļ?
cliFile="time.txt"            #?????????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\8.GRGsur")     #???ù???Ŀ¼

#??ȡ?????ļ????????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data2=data
data=t(data)
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#??ȡ?????????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

#???ݺϲ?
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#?Ի???????ѭ????????Ԥ?????صĻ???
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	#cox????
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<0.05){
		outTab=rbind(outTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				        )
	}
}

#?????????صĽ???
write.table(outTab, file="uniCox.txt", sep="\t", row.names=F, quote=F)

#???????????????????ı???��
sigGenes=as.vector(outTab[,1])
uniSigExp=data2[sigGenes,]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

#?????????????????????????ݺ????????ݺϲ????ļ?
uniSigExpTime=rt[,c("futime","fustat",sigGenes)]
uniSigExpTime=cbind(id=row.names(uniSigExpTime),uniSigExpTime)
write.table(uniSigExpTime, file="uniSigExpTime.txt", sep="\t", row.names=F, quote=F)

############????ɭ??ͼ????############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
	#??ȡ?????ļ?
	rt <- read.table(coxFile,header=T,sep="\t",check.names=F,row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#????ͼ??
	height=nrow(rt)/13+5
	pdf(file=forestFile, width = 7,height = height)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#????ɭ??ͼ???ߵĻ?????Ϣ
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		
	#????ɭ??ͼ
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
	axis(1)
	dev.off()
}

#???ú???, ?Ե????صĽ??????п??ӻ?, ????ɭ??ͼ
bioForest(coxFile="uniCox.txt", forestFile="forest.pdf", forestCol=c("red","blue"))


####