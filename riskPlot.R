

#install.packages("pheatmap")


library(pheatmap)           #???ð?
riskFile="risk.all.txt"     #?????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\24.riskHeatmap")      #???ù???Ŀ¼

#??ȡ?????ļ?
rt=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$riskScore),]      #???ݲ??˵ķ??յ÷ֶ???Ʒ????????

#??????ͼע?͵???ɫ
ann_colors=list()
bioCol=c("#0088FF", "#FF5555")
names(bioCol)=c("low", "high")
ann_colors[["risk"]]=bioCol

#???Ʒ?????ͼ
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(risk=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="riskHeatmap.pdf", width=7, height=4)
pheatmap(rt1, 
		 annotation=annotation,
		 annotation_colors = ann_colors, 
		 cluster_cols = F,
		 show_colnames = F,
		 cluster_rows = T,
		 color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
		 scale="row",
		 fontsize_col=8,
		 fontsize=7,
		 fontsize_row=8)
dev.off()



