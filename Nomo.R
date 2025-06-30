

#install.packages("survival")
#install.packages("survminer")
#install.packages("regplot")
#install.packages("rms")


#???ð?
library(survival)
library(regplot)
library(rms)
library(survminer)

riskFile="risk.all.txt"     #?????ļ?
cliFile="clinical.txt"      #?ٴ??????ļ?
setwd("C:\\Users\\17860\\Desktop\\OV\\26.Nomo")     #???ù???Ŀ¼

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#??ȡ?ٴ??????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#?ϲ?????
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)

#????????ͼ
# 移除 'Gender' 变量，宫颈癌全是女性所以要移除
rt <- rt[, !(names(rt) %in% c("Gender"))]

# 运行 Cox 回归模型
res.cox <- coxph(Surv(futime, fustat) ~ . , data = rt)

res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=TRUE,
              droplines=TRUE,
              observation=rt[2,],
              rank="sd",
              failtime = c(1,3,5),
              prfail = F)
dev.copy2pdf(file="Nomo.pdf", width=8, height=6, out.type="pdf")

#????????ͼ?ķ??յ÷?
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#У׼????
pdf(file="calibration.pdf", width=5, height=5)
#1??У׼????
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
	 xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3??У׼????
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5??У׼????
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
	   col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

#?ۼƷ???????
nomoRisk=ifelse(rt$Nomogram>median(rt$Nomogram), "High", "Low")
fit=survfit(Surv(futime, fustat) ~ nomoRisk, data=rt)
gg=ggsurvplot(fit,
           conf.int = T,
           risk.table.col="strata",
           ggtheme = theme_bw(),
           #palette = "lancet",
           fun = 'cumhaz')
pdf(file="cumulative.pdf", width=5, height=4.8, onefile=F)
print(gg)
dev.off()



