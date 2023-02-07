
#install.packages("survival")
#install.packages("caret")
#install.packages("glmnet")
#install.packages("survminer")
#install.packages("timeROC")


#引用包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
library(randomForestSRC)

coxPfilter=0.01        #单因素cox方法显著性的过滤标准
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA\\12.model")      #设置工作目录

#读取输入文件
tcga.exp=read.table("TCGAexpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
cgga.exp <- read.table("CGGAexpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
REMBRANDTE.exp <- read.table("REMBRANDTexpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
tcga.exp$futime=tcga.exp$futime/365
cgga.exp$futime=cgga.exp$futime/365
REMBRANDTE.exp$futime=REMBRANDTE.exp$futime/365

############定义森林图函数############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width=7, height=6)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  LOGindex = 10 
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol[1],forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,10^a1)
  dev.off()
}
############绘制森林图函数############

############对数据进行分组#############
  train=tcga.exp
  #单因素cox分析
  outUniTab=data.frame()
  sigGenes=c("futime","fustat")
  for(i in colnames(train[,3:ncol(train)])){
    #cox分析
    cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    
    #保留显著性基因
    if(coxP<coxPfilter){
      sigGenes=c(sigGenes,i)
      outUniTab=rbind(outUniTab,
                      cbind(id=i,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }
  uniSigExp=train[,sigGenes]
  uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
  if(ncol(uniSigExp)<10){next}
  
  #随机森林
  pbc.obj <- rfsrc(Surv(futime, fustat) ~ ., data = uniSigExp,
                   ntree = 1000,
                   mtry = length(colnames(uniSigExp))/3,
                   nodesize=15,
                   samptype="swor",
                   sampsize=237,
                   block.size = 10,
                   importance = "random",
                   splitrule="logrank",
                   nsplit = 10
  )
  vs.pbc <- var.select(object = pbc.obj)
  topvars <- vs.pbc$topvars
  RandomExp <- uniSigExp[,c("futime", "fustat", topvars)]
  RandomExpOut <- cbind(id=row.names(RandomExp), RandomExp)

  #############构建COX模型#############
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = RandomExp)
  multiCox=step(multiCox, direction = "both")
  multiCoxSum=summary(multiCox)
  
  #输出模型的公式
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  outMultiTab=outMultiTab[,1:2]
  
  #输出TCGA风险文件
  TCGA.riskScore=predict(multiCox,type="risk",newdata=tcga.exp)      #利用train得到模型预测train样品风险
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  medianTrainRisk=median(TCGA.riskScore)
  risk=as.vector(ifelse(TCGA.riskScore>medianTrainRisk,"high","low"))
  TCGARiskOut=cbind(id=rownames(cbind(tcga.exp[,outCol],TCGA.riskScore,risk)),cbind(tcga.exp[,outCol],riskScore=TCGA.riskScore,risk))
  
  #输出CGGA组风险文件
  CGGA.riskScore=predict(multiCox,type="risk",newdata=cgga.exp)     #利用train得到模型预测test样品风险
  riskTest=as.vector(ifelse(CGGA.riskScore>medianTrainRisk,"high","low"))
  CGGARiskOut=cbind(id=rownames(cbind(cgga.exp[,outCol],CGGA.riskScore,riskTest)),cbind(cgga.exp[,outCol],riskScore=CGGA.riskScore,risk=riskTest))
  
  #输出REMBRANDTE组风险文件
  REMBRANDTE.riskScore=predict(multiCox,type="risk",newdata=REMBRANDTE.exp)     #利用train得到模型预测test样品风险
  REMBRANDTE.riskTest=as.vector(ifelse(REMBRANDTE.riskScore>medianTrainRisk,"high","low"))
  REMBRANDTEOut=cbind(id=rownames(cbind(REMBRANDTE.exp[,outCol],REMBRANDTE.riskScore,REMBRANDTE.riskTest)),cbind(REMBRANDTE.exp[,outCol],riskScore=REMBRANDTE.riskScore,risk=REMBRANDTE.riskTest))
  
  #比较高低风险组的生存差异，得到差异的pvalue	
  diff=survdiff(Surv(futime, fustat) ~risk,data = TCGARiskOut)
  pValue=1-pchisq(diff$chisq, df=1)
  CGGA.diff=survdiff(Surv(futime, fustat) ~risk,data = CGGARiskOut)
  CGGA.pValue=1-pchisq(CGGA.diff$chisq, df=1)
  REMBRANDTE.diff=survdiff(Surv(futime, fustat) ~risk,data = REMBRANDTEOut)
  REMBRANDTE.pValue=1-pchisq(REMBRANDTE.diff$chisq, df=1)
  
  
  #ROC曲线下面积
  predictTime=1    #预测时间
  roc=timeROC(T=TCGARiskOut$futime, delta=TCGARiskOut$fustat,
              marker=TCGA.riskScore, cause=1,
              times=c(predictTime), ROC=TRUE)
  CGGA.roc=timeROC(T=CGGARiskOut$futime, delta=CGGARiskOut$fustat,
                   marker=CGGARiskOut$riskScore, cause=1,
                   times=c(predictTime), ROC=TRUE)	
  REMBRANDTE.roc=timeROC(T=REMBRANDTEOut$futime, delta=REMBRANDTEOut$fustat,
                         marker=REMBRANDTEOut$riskScore, cause=1,
                         times=c(predictTime), ROC=TRUE)
    #输出单因素结果
    write.table(outUniTab,file="uni.trainCox.txt",sep="\t",row.names=F,quote=F)
    write.table(uniSigExpOut,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)
    bioForest(coxFile="uni.trainCox.txt",forestFile="uni.foreast.pdf",forestCol=c("red","green"))
    #lasso结果
    write.table(RandomExpOut,file="Random.SigExp.txt",sep="\t",row.names=F,quote=F)
    pdf("Random.pdf", width=11.5, height=6)
    plot(pbc.obj)
    dev.off()
    #输出多因素结果
    write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
    write.table(TCGARiskOut,file="risk.tcga.txt",sep="\t",quote=F,row.names=F)
    write.table(CGGARiskOut,file="risk.cgga.txt",sep="\t",quote=F,row.names=F)
    write.table(REMBRANDTEOut,file="risk.REMBRANDTE.txt",sep="\t",quote=F,row.names=F)
    #所有样品的风险值
    allRiskOut=rbind(TCGARiskOut, CGGARiskOut, REMBRANDTEOut)
    write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)