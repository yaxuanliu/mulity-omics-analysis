#install.packages("survival")
#install.packages("survminer")

library(survival)
library("survminer")
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA\\16.Risk survival")              #设置工作目录

bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                   #读取输入文件
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=paste0("p=",pValue),
                     pval.size=4,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     surv.median.line = "hv",
                     palette=c("#996699", "#336699"),
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 7,height =6)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="risk.cgga.txt",outFile="risk.cgga.pdf")
bioSurvival(inputFile="risk.REMBRANDTE.txt",outFile="risk.REMBRANDTE.pdf")
bioSurvival(inputFile="risk.tcga.txt",outFile="risk.tcga.pdf")
