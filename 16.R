#引用包
library(limma)
library(ggpubr)

tciaFile="TCIA.txt"        #免疫打分文件
expFile="risk.tcga.txt"      #表达数据文件
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA2\\36.TCIA")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(rt)
rownames(data)=substr(rownames(data),start = 1,stop = 12)

#读取TCIA的打分文件
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(ips), row.names(data))
ips=ips[sameSample, , drop=F]
data=data[sameSample, "risk", drop=F]
data=cbind(ips, data)

#设置比较组
data$risk=as.factor(data$risk)
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对TCIA打分进行循环,分别绘制小提琴图
for(i in colnames(data)[1:(ncol(data)-1)]){
  rt=data[,c(i, "risk")]
  gg1=ggviolin(data, x="risk", y=i, fill = "risk", 
               xlab="", ylab=i,
               legend.title="risk",
               palette = c("#336699", "#996699"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0(i, ".pdf"), width=4.8, height=4.25)
  print(gg1)
  dev.off()
}
dataOUT=cbind(id=row.names(data),data)
write.table(dataOUT,file = "TCIA.result.txt",sep = "\t",row.names = F,quote = F)
