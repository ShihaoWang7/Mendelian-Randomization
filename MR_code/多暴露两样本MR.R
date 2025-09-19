#设置工作目录
setwd("C:\\Users\\wsh20\\Desktop")
#读入TwoSampleMR包
library(TwoSampleMR)
# 自动一次提取多个变量的工具变量，只需要把id填入
Multiple_exp_dat <- extract_instruments(
  outcomes = c("ukb-e-24006_EAS", "ukb-e-24005_EAS", "ukb-e-24004_EAS","ukb-e-24003_EAS","ukb-e-24024_EAS"),p1 = 5e-06,
  clump = TRUE, r2 = 0.001,
  kb = 10000, access_token = NULL
)
#本地数据导入
pm2.5.exp_dat <- read.table("pm2.5.exp_dat.txt",sep = "\t",header = T)
pm10.exp_dat <- read.table("pm10.exp_dat.txt",sep = "\t",header = T)
NO2.exp_dat <- read.table("NO2.exp_dat.txt",sep = "\t",header = T)
NOx.exp_dat <- read.table("NOx.exp_dat.txt",sep = "\t",header = T)
noise.exp_dat <- read.table("noise.txt",sep = "\t",header = T)
#合并本地数据
#两个数据框
Multiple_exp_dat <- rbind(NO2.exp_dat,NOx.exp_dat)
#多个数据框
library(tidyverse)
df1 <- pm2.5.exp_dat
df2 <- pm10.exp_dat
df3 <- NO2.exp_dat
df4 <- NOx.exp_dat
df5 <- noise.exp_dat
l <- data.frame()
for(i in 1:5){
  df<-get(paste0("df",i))
  l<-rbind(l,df)
}
Multiple_exp_dat <- l 
write.table(file = "pollution_snps.txt",Multiple_exp_dat,sep = "\t",quote = F)
# 提取多个变量合并SNP的结局
# 当然这边也可以填入多个研究结局 outcomes=c("ieu-a-2","ieu-a-7")
Multiple_exp_dat <- read.table("pollution_snps.txt",sep = "\t",header = T)
outcome <- extract_outcome_data(snps = Multiple_exp_dat$SNP,
                                outcomes = "ieu-a-1189")
rm(list = ls())
#####提取结局snp(本地数据)
data <- fread("pollution_snps.txt",sep = "\t",header = T)
data1 <- fread("ASD.txt",sep = "\t",header = T,fill=TRUE)
data2 <- merge(data,data1,by = "SNP")

write.table(file = "outcome.txt",data2,sep = "\t",quote = F)

outcome_dat <- read_outcome_data(snps = data$SNP,
                                 filename ="outcome.txt",
                                 sep = "\t",
                                 eaf_col = "Freq",
                                 snp_col = "SNP", 
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 pval_col = "P",
                                 samplesize_col = "N"
)
outcome_dat$outcome <- "TS"
outcome_dat$samplesize.outcome <- 72517

# harmonise # 【以这个函数,查看原理】
Multiple_har <- harmonise_data(data,outcome_dat,action = 2)
outcome_dat$samplesize.outcome <- 500199
Multiple_har <- Multiple_har[-c(43,44,45),]
write.table(file = "Multiple_har.txt",Multiple_har,sep = "\t",quote = F)
Multiple_har <- fread("Multiple_har.txt",sep = "\t",header = T)
Multiple_har <- Multiple_har[-c(32,33,34),]
# 可能存在双向问题再用，进行过滤
Multiple_har=steiger_filtering(Multiple_har)

# 进行mr分析
Multiple_res <- mr(Multiple_har)
OR <-generate_odds_ratios(Multiple_res)
write.table(file = "TS_MR.txt",OR,sep = "\t",quote = F)
# 进行MRPRESSO
Multiple_har <- as.data.frame(Multiple_har)
res_MRPRESSO <- run_mr_presso(dat = Multiple_har)
#将列表作文本文件输出
capture.output(res_MRPRESSO, file = "TS_res_MRPRESSO.txt")
# 进行异质性
heterogeneity <- mr_heterogeneity(Multiple_har)
write.table(heterogeneity,file = "TS_异质性.txt",sep = "\t",quote = F)
# 进行水平多效性
pleiotropy <- mr_pleiotropy_test(Multiple_har)
write.table(pleiotropy,file = "TS_多效性.txt",sep = "\t",quote = F)

#分割数据画图
x <- Multiple_har[grep(pattern = "ukb-b-10817",Multiple_har$id.exposure),]
y <- Multiple_har[grep(pattern = "ukb-b-18469",Multiple_har$id.exposure),]
z <- Multiple_har[grep(pattern = "ukb-b-12417",Multiple_har$id.exposure),]
a <- Multiple_har[grep(pattern = "ukb-b-9942",Multiple_har$id.exposure),]
b <- Multiple_har[grep(pattern = "ukb-b-19490",Multiple_har$id.exposure),]
x1 <- x[-3,]
a1 <- a[-5,]
y1 <- y[-14,]
z1 <- z[-3,]
b1 <- b[-5,] 
#逐个剔除检验
single <- mr_leaveoneout(x)
pdf(file="PM2.5_归一法.pdf")
mr_leaveoneout_plot(single)
dev.off()

single <- mr_leaveoneout(y)
pdf(file="PM10_归一法.pdf")
mr_leaveoneout_plot(single)
dev.off()

single <- mr_leaveoneout(z)
pdf(file="NOx_归一法.pdf")
mr_leaveoneout_plot(single)
dev.off()

single <- mr_leaveoneout(a)
pdf(file="NO2_归一法.pdf")
mr_leaveoneout_plot(single)
dev.off()

single <- mr_leaveoneout(b)
pdf(file="noise_归一法.pdf")
mr_leaveoneout_plot(single)
dev.off()
#散点图
pdf(file="散点图.pdf")
mr_scatter_plot(Multiple_res,Multiple_har)
dev.off()
#森林图 漏斗图

res_single <- mr_singlesnp(x)
pdf(file="PM2.5_森林图.pdf")
mr_forest_plot(res_single)
dev.off()
pdf(file="PM2.5_漏斗图.pdf")
mr_funnel_plot(res_single)
dev.off()

res_single <- mr_singlesnp(y)
pdf(file="PM10_森林图.pdf")
mr_forest_plot(res_single)
dev.off()
pdf(file="PM10_漏斗图.pdf")
mr_funnel_plot(res_single)
dev.off()

res_single <- mr_singlesnp(z)
pdf(file="NOx_森林图.pdf")
mr_forest_plot(res_single)
dev.off()
pdf(file="NOx_漏斗图.pdf")
mr_funnel_plot(res_single)
dev.off()

res_single <- mr_singlesnp(a)
pdf(file="NO2_森林图.pdf")
mr_forest_plot(res_single)
dev.off()
pdf(file="NO2_漏斗图.pdf")
mr_funnel_plot(res_single)
dev.off()

res_single <- mr_singlesnp(b)
pdf(file="noises_森林图.pdf")
mr_forest_plot(res_single)
dev.off()
pdf(file="noises_漏斗图.pdf")
mr_funnel_plot(res_single)
dev.off()

#清除空间
rm(list = ls())