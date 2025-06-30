

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("pheatmap")
#install.packages("vioplot")
#install.packages("corrplot")


#???ð?
library(limma)
library(pheatmap)
library(ggpubr)
library(vioplot)
library(corrplot)

immFile="CIBERSORT-Results.txt"     #????ϸ???????Ľ????ļ?
riskFile="risk.all.txt"             #?????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\28.immunePlot")     #???ù???Ŀ¼

#??ȡ????ϸ???????Ľ????ļ??????????ݽ???????
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=rownames(risk)[risk[,"risk"]=="low"]
highSample=rownames(risk)[risk[,"risk"]=="high"]

#?ߵͷ???????????ϸ????��
lowSameSample=intersect(row.names(data), lowSample)
highSameSample=intersect(row.names(data), highSample)
data=t(data[c(lowSameSample,highSameSample),])
conNum=length(lowSameSample)
treatNum=length(highSameSample)

##########??????״ͼ##########
pdf("barplot.pdf", width=25, height=12)
col=rainbow(nrow(data),s=0.7,v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,ylab="Relative Percent",xaxt="n",yaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()

##########??????????ͼ??##########
pdf(file="corHeatmap.pdf", width=11, height=11)
par(oma=c(0.5,1,1,1.2))
corData=t(data)
corData=corData[,colMeans(corData)>0]
M=cor(corData)
corrplot(M,
         order="hclust",          #????ϸ??????????ʽ
         method = "color",        #ͼ??չʾ????ʽ
         diag = TRUE,             #?Ƿ?չʾ?Խ???
         tl.col="black",          #??????ɫ
         addCoef.col = "black",   #????ϵ??????????ɫ
         number.cex=0.75,         #????ϵ????????С
         col=colorRampPalette(c("blue", "white", "red"))(50))    #ͼ?ε???ɫ,?????غ?ɫ????????��ɫ
dev.off()

##########????С????ͼ##########
rt=t(data)
pdf("vioplot.pdf", height=8, width=12)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63), ylim=c(min(rt), max(rt)+0.02),
     main="", xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#??ÿ??????ϸ??ѭ????????С????ͼ???ͷ???????ɫ??ʾ???߷????ú?ɫ??ʾ
for(i in 1:ncol(rt)){
	  if(sd(rt[1:conNum,i])==0){
	    rt[1,i]=0.00001
	  }
	  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
	    rt[(conNum+1),i]=0.00001
	  }
	  lowData=rt[1:conNum,i]
	  highData=rt[(conNum+1):(conNum+treatNum),i]
	  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = 'green')
	  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(lowData, highData)
	  p=wilcoxTest$p.value
	  mx=max(c(lowData,highData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topleft", 
       c("Low risk", "High risk"),
       lwd=3.5, bty="n", cex=1.2,
       col=c("green","red"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)
dev.off()


######V