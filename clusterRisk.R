
#install.packages("ggplot2")
#install.packages("ggalluvial")
#install.packages("dplyr")
#install.packages("ggpubr")


#???ð?
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(ggpubr)

cluFile="GRGcluster.txt"      #???͵Ľ????ļ?
riskFile="risk.all.txt"       #?????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\25.clusterRisk")     #???ù???Ŀ¼

#??ȡ?????ļ?
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
sameSample=intersect(row.names(cluster), row.names(risk))
data=cbind(risk[sameSample,,drop=F], cluster[sameSample,,drop=F])

#######????ͼ########
#???ñȽ???
data$ARGcluster=factor(data$GRGcluster, levels=levels(factor(data$GRGcluster)))
group=levels(factor(data$GRGcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#????ͼ????ɫ
bioCol=c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
         "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
bioCol=bioCol[1:length(levels(factor(data$GRGcluster)))]
	
#????????ͼ
boxplot=ggboxplot(data, x="GRGcluster", y="riskScore", color="GRGcluster",
			      xlab="GRGcluster",
			      ylab="Risk score",
			      legend.title="GRGcluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)

#????ͼ??
pdf(file="clusterRisk.pdf", width=5.5, height=4.5)
print(boxplot)
dev.off()
#######????ͼ########


#׼??ɣ??ͼ?????ļ?
rt=data[,c("GRGcluster", "risk", "fustat")]
colnames(rt)=c("GRGcluster", "Risk", "Fustat")
rt[,"Fustat"]=ifelse(rt[,"Fustat"]==0, "Alive", "Dead")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#????ɣ??ͼ
pdf(file="ggalluvial.pdf", width=6, height=5.5)
mycol=rep(c("#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
            "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #??????????ɫ??forward˵????????ɫ??ǰ??????״ͼһ?£?backward˵????????ɫ??????????״ͼһ?¡?
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size=3??????????С
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()



