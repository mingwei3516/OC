

#install.packages("RCircos")


library("RCircos")       #???ð?
setwd("C:\\Users\\17860\\Desktop\\OV\\12.RCircos")    #???ù???Ŀ¼

#??ʼ??Ȧͼ
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t")
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#????Ȧͼ????
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.6
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

#?????ļ?
pdf(file="RCircos.pdf", width=8, height=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#??ȡ?????????ļ???????ɢ??ͼ
RCircos.Scatter.Data=read.table("Rcircos.scatter.txt", header=T, sep="\t", check.names=F)
data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)

#??ȡ????ע???ļ?????ע??????????
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()




