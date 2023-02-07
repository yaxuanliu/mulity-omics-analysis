
library(maftools)       #引用包
setwd("C:\\Users\\86188\\Desktop\\GBM-scRNA2\\25.maftools")      #设置工作目录

#读取风险文件,得到注释文件
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#读取基因突变的文件，选取突变频率最高的15个基因进行可视化
geneNum=20      #设置基因的数目
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#注释的颜色
ann_colors=list()
col=c("#336699", "#996699")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

#定义颜色
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

#绘制高风险组瀑布图
maf_high=read.maf(maf="high.maf")    #读取高风险组的突变文件
pdf(file="high.pdf", width=9, height=6)
oncoplot(maf=maf_high, top = 20, colors = vc_cols, keepGeneOrder=T, showTumorSampleBarcodes = T)
# oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#绘制低风险组瀑布图
maf_low=read.maf(maf="low.maf")    #读取低风险组的突变文件
pdf(file="low.pdf", width=9, height=6)
oncoplot(maf=maf_low, top = 20, colors = vc_cols, keepGeneOrder=T, showTumorSampleBarcodes = T)
dev.off()

#高低风险比较
HLrisk <- mafCompare(m1 = maf_high, m2 = maf_low, m1Name = 'High risk', m2Name = 'Low risk', minMut = 5)
pdf(file="forest.pdf", width=9, height=6)
forestPlot(mafCompareRes = HLrisk, pVal = 0.05,lineWidth = 2, color = c('#336699', '#996699'), geneFontSize = 0.8)
dev.off()

#体细胞突变互动
pdf(file="somaticInteract_high.pdf", width=7.5, height=6)
somaticInteractions(maf = maf_high, top = 25,showSum = F,pvalue = c(0.05, 0.01))
dev.off()
pdf(file="somaticInteract_low.pdf", width=7.5, height=6)
somaticInteractions(maf = maf_low, top = 25,showSum = F,pvalue = c(0.05, 0.01))
dev.off()

#预测癌症驱动基因
luad.sig <- oncodrive(maf=maf_high, minMut=5, AACol="HGVSp_Short", pvalMethod="zscore")
pdf(file="Oncodrive_high.pdf", width=7.5, height=6)
plotOncodrive(res = luad.sig, fdrCutOff = 0.1, useFraction = TRUE)
dev.off()
luad.sig_low <- oncodrive(maf=maf_low, minMut=5, AACol="HGVSp_Short", pvalMethod="zscore")
plotOncodrive(res = luad.sig_low, fdrCutOff = 0.1, useFraction = TRUE)

#Oncogenic 信号通路
Hpath=OncogenicPathways(maf = maf_high)
write.table(Hpath, file="Hpath.txt", sep="\t", quote=F, row.names=F)
Lpath=OncogenicPathways(maf = maf_low)
write.table(Lpath, file="Lpath.txt", sep="\t", quote=F, row.names=F)
