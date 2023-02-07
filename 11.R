
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA2\\23.gsesurvival")    #工作目录（需修改）

library("survminer")
library(survival)
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365                                        #如果以月为单位，除以30；以年为单位，除以365
outTab=data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,3)
  #pValue=format(pValue, scientific = TRUE)
  
  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)
  
    ggsurvplot(fit, 
             data=rt,
             conf.int=F,
             pval=paste0("p=",pValue),
             pval.size=4,
             risk.table=F,
             legend=c(0.8, 0.9),
             legend.labs=c("High score", "Low score"),
             legend.title="",
             xlab="Time(years)",
             break.time.by = 3,
             surv.median.line = "hv",
             palette=c("#996699", "#336699"),
  )
  ggsave(filename=paste(gene,".survival.pdf",sep=""),
         width = 4.5,            
         height =4            
  )
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)
