library(limma)

expFile="gsvaOut.txt"                                          #表达数据文件
MarFile="Hallmark.txt"
cliFile="time.all.txt"                                                 #临床数据
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA2\\23.1.presurvival")       #工作目录（需修改）

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=t(data)

#读取关键通路
Hallmark=read.table(MarFile,sep="\t",check.names=F,header=F)     #读取临床文件

#通路取交集
data=data[,Hallmark[,1]]

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)