#引用包
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)
expFile="symbol.all.txt"      #表达输入文件
riskFile="risk.all.txt"       #风险输入文件
geneFile="gene.txt"       #免疫检查点的基因文件
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA2\\34.checkpoint")     #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因文件
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
# data=log2(data+1)

#删除正常样品
# group=sapply(strsplit(row.names(data),"\\-"),"[",4)
# group=sapply(strsplit(group,""),"[",1)
# group=gsub("2","1",group)
# data=data[group==0,]
# row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
# data=avereps(data)

#合并数据
risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c("riskScore",sameGene)]
rtOUT=cbind(id=row.names(rt1),rt1)
write.table(rtOUT,file = "Corresult.txt",row.names = F,quote = F,sep = '\t')

#生成相关性热图
pdf("correlation.pdf",height=8,width=8)         #保存图片的文件名称
corrplot(corr=cor(rt1),
         method = "color",
         #hclust.method="ward.D2",
         order = "original",
         addCoef.col = "black",
         number.cex = 0.33,
         diag = TRUE,
         type = "full",
         tl.col="black",
         tl.cex = 0.8,
         col=colorRampPalette(c("#336699", "white", "#996699"))(50),
)

dev.off()
