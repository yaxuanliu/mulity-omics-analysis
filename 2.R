
#引用包
library(survival)
library(survminer)
library(timeROC)

riskFile="risk.tcga.txt"     #风险文件
riskFiletest="risk.cgga.txt"     #风险文件
riskFiletest2="risk.REMBRANDTE.txt"
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA\\13.ROC")     #修改工作目录

#读取风险输入文件
rt=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 3 5年的ROC曲线######
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.tcga.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()

#读取风险输入文件
rt=read.table(riskFiletest, header=T, sep="\t", check.names=F, row.names=1)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 3 5年的ROC曲线######
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.cgga.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()

#读取风险输入文件
rt=read.table(riskFiletest2, header=T, sep="\t", check.names=F, row.names=1)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 3 5年的ROC曲线######
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.rem.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()

