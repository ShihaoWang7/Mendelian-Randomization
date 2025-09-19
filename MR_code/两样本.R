#打开R包
library(TwoSampleMR)

#设置工作环境
setwd("D:/5-孟德尔随机化GWAS数据/BMI/SNP_gwas_mc_merge_nogc.tbl.uniq")

#读取bmi数据
a<-read.table("SNP_gwas_mc_merge_nogc.tbl.txt",header = T)
View(a)

#相关性设置，并将文件放到TwoSampleMR包所在文件夹
b<-subset(a,p<5e-08)
write.csv(b, file="exposure.csv")

#读取exposure数据
bmi<-system.file("3- Teaching cases/exposure.csv",package = "TwoSampleMR")
bmi_exp_dat_clumped<-read_exposure_data(filename = bmi,sep = ",",snp_col = "SNP",beta_col = "beta",se_col = "se",effect_allele_col = "effect_allele",other_allele_col = "other_allele",eaf_col = "eaf",clump = TRUE)

#设置工作环境
setwd("D:/5-孟德尔随机化GWAS数据/Alzheimer IGAP/IGAP_summary_statistics")

#读取Alzheimer数据
c<-read.table("IGAP_stage_1.txt",header = T)

#将暴露SNP从结局中提取出来
d<-merge(bmi_exp_dat_clumped,c,by.x = "SNP",by.y = "MarkerName")

#读取outcome数据
outcome_dat<-read_outcome_data(snps = bmi_exp_dat_clumped$SNP,filename = "outcome.csv",sep = ",",snp_col = "SNP",beta_col = "beta",se_col = "se",effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "p")
dat<-harmonise_data(exposure_dat = bmi_exp_dat_clumped,outcome_dat = outcome_dat)

#Perform MR
mr(dat)
generate_odds_ratios(mr_res = mr(dat))
mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
mr_heterogeneity(dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
mr_pleiotropy_test(dat)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))


setwd("F://")
a <- read.table("bbj-a-10.txt",header=T)
View(a)
b <- strsplit(a$X0.0002103.0.006598.0.0111736.0.1517.rs143225517,":") 
View(b)
a$p <- c()
 for (i in 1:nrow(a)) {a$p[i] <- b[[i]][4]
 }

