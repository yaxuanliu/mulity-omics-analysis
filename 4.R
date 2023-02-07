
library(dplyr)
library(ggplot2)
setwd("C:\\Users\\86188\\Desktop\\17.3.clinical")    #设置工作目录
dat=read.table("clinical.txt",header=T,sep="\t",check.names=F,row.names = 1) 
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
dat[,"Age"]=ifelse(dat[,"Age"]=="unknow", "unknow", ifelse(dat[,"Age"]>60,">60","<=60"))

#样本合并
sameSample=intersect(row.names(dat),row.names(risk))
risk=risk[sameSample,]
dat=dat[sameSample,]
dat=cbind(dat,risk=risk[,"risk"])

# 按Risk分成High和Low，计算各列数值。
gname <- "risk"
vname <- setdiff(colnames(dat), gname)
pie.high <- pie.low <- list()
fisher.p <- c()

for (i in vname) {
  tmp <- table(dat[,gname], dat[,i])
  p <- format(fisher.test(tmp)$p.value,digits = 2)
  names(p) <- i
  fisher.p <- c(fisher.p, p)
  
  pie.dat <- 
    tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
  
  # 表格内的两行对应Risk的两类：Risk high和Risk low
  pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "high"),]
  pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "low"),]
}
# 
# rthigh=rt[rt$risk=="high",]
# rthigh=rthigh[,1:ncol(rthigh)-1]
# rtlow=rt[rt$risk=="low",]
# rtlow=rtlow[,1:ncol(rtlow)-1]

#定义颜色
black  <- "#1E1E1B"
blue   <- "#3C4E98"
yellow <- "#E4DB36"
orange <- "#E19143"
green  <- "#57A12B"
cherry <- "#8D3A86"
red <- "#ff3333"

# 创建颜色
Age.col <- c("grey80",black)
Gender.col <- alpha(cherry, c(0.5, 1))
PRS_type.col <- alpha(blue, c(0.5, 1))
Radio_therapy.col <- alpha(green, c(0.5, 1))
Chemo_therapy.col <- c(yellow, orange)
IDH_mutation.col <- alpha(red, c(0.5, 1))
etion.col <- alpha("#339999",c(0.5, 1))
MGMTp_methylation.col <- alpha("#33cccc", c(0.5, 1))

# 硬核base plot一块一块画，当然也可以把其中的pie chart提取出来后期AI或者PPT拼接也是比较方便的
pdf("pieTable.pdf",width = 10, height = 5)
showLayout <- F # 默认不在最终pdf的首页显示layout结构，不过建议初次绘制的时候改为TRUE看一下，方便理解

# 设置画面布局，相同数字代表同一区块，数字越多代表该区块所占面积越大（一共25个区域）
layout(matrix(c( 1, 1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,7, 7, 7,8,8,8,9,9,9,
                 10,10,10, 11,11,11, 12,12,12, 13,13,13,14,14,14,15,15, 15,16,16,16,17,17,17, 18,18,18,
                 10,10,10, 11,11,11, 12,12,12, 13,13,13,14,14,14,15,15, 15,16,16,16,17,17,17, 18,18,18,
                 19,19,19, 20,20,20, 21,21,21,22,22,22, 23,23,23, 24,24,24,25,25,25, 26,26,26,27,27,27,
                 19,19,19, 20,20,20, 21,21,21,22,22,22, 23,23,23, 24,24,24,25,25,25, 26,26,26,27,27,27,
                 28,28,28,29,29,29, 30,30,30, 31,31,31, 32,32,32, 33,33,33,34,34,34, 35,35,35,36,36,36,
                 37,37,37,37,37,37,37,37,37,37,37,37, 37,37,37, 37,37,37,37,37,37, 37,37,37,37,37,37),
              byrow = T,nrow = 7))

if(showLayout) {
  layout.show(n = 37) # 直观展示画布分布
}

#-------------------------#
# 画布区域1-6：绘制图抬头 #
#-------------------------#

par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) # 基础参数，各边界距离为0
plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "GBM",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Age",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Gender",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "PRS_type",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Radio_therapy",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Chemo_therapy",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "IDH_mutation",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "1p19q_1etion",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "MGMTp_methylation",cex = 1.5, col = "white") # 显示图标题

#--------------------------------------#
# 画布区域7-12：绘制High组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "High\n(n = 355)",cex = 1.5, col = "white") # 显示图标题

# High group
pie(pie.high$Age$Pct, 
    col = Age.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Gender$Pct, 
    col = Gender.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$PRS_type$Pct, 
    col = PRS_type.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Radio_therapy$Pct, 
    col = Radio_therapy.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Chemo_therapy$Pct, 
    col = Chemo_therapy.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$IDH_mutation$Pct, 
    col = IDH_mutation.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$etion$Pct, 
    col = etion.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$MGMTp_methylation$Pct, 
    col = MGMTp_methylation.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------------#
# 画布区域13-18：绘制Low组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Low\n(n = 191)",cex = 1.5, col = "white") # 显示图标题

# Low group
pie(pie.low$Age$Pct, 
    col = Age.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Gender$Pct, 
    col = Gender.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$PRS_type$Pct, 
    col = PRS_type.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Radio_therapy$Pct, 
    col = Radio_therapy.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Chemo_therapy$Pct, 
    col = Chemo_therapy.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$IDH_mutation$Pct, 
    col = IDH_mutation.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$etion$Pct, 
    col = etion.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$MGMTp_methylation$Pct, 
    col = MGMTp_methylation.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------#
# 画布区域19-24：绘制空抬头和p值 #
#--------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Age"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Gender"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n",# 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["PRS_type"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Radio_therapy"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Chemo_therapy"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["IDH_mutation"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["etion"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["MGMTp_methylation"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#----------------------#
# 画布区域25：绘制图例 #
#----------------------#

plot(0,0,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
legend("topleft",
       legend = c("<=65",">65",
                  "Female","Male",
                  "False","True",
                  "False","True",
                  "False","True",
                  "False","True",
                  "False","True",
                  "False","True"),
       fill = c(Age.col,
                Gender.col,
                PRS_type.col,
                Radio_therapy.col,
                Chemo_therapy.col,
                IDH_mutation.col,
                etion.col,
                MGMTp_methylation.col),
       border = NA, # 图例颜色没有边框
       bty = "n", # 图例没有边框
       cex = 1.5,
       #box.lwd = 3,
       x.intersp = 0.05,
       y.intersp = 1,
       text.width = 0.063, # 图例的间隔
       horiz = T) # 图例水平放置
abline(h = par("usr")[1], col = "black") # 底部封上黑线
abline(h = par("usr")[3], col = "black") # 底部封上黑线
abline(v = par("usr")[2], col = "black") # 右侧封上黑线
abline(h = par("usr")[4], col = "black") # 底部封上黑线
# 关闭图像句柄
invisible(dev.off())



