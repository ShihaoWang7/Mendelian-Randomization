#设置工作目录
setwd("C:\\Users\\wsh20\\Desktop")
rm(list = ls())
#读入数据包
library(TwoSampleMR)
library(data.table)

data1 <- fread("PTSD.txt",sep = "\t",quote = F)

data1 <- extract_outcome_data(snps = data$SNP, outcomes = 'ukb-b-10817')
data1$outcome <- "pm2.5"
data1 <- data1[,-c(13:19)]
               
data2 <- extract_outcome_data(snps = data$SNP, outcomes = 'ukb-b-18469')
data2$outcome <- "pm10"
data2 <- data2[,-c(13:19)]


data3 <- extract_outcome_data(snps = data$SNP, outcomes = 'ukb-b-12417')
data3$outcome <- "NOx"
data3 <- data3[,-c(13:19)]

data4 <- extract_outcome_data(snps = data$SNP, outcomes = 'ukb-b-9942')
data4$outcome <- "NO2"
data4 <- data4[,-c(13:19)]

data5 <- extract_outcome_data(snps = data$SNP, outcomes = 'ukb-b-19490')
data5$outcome <- "noise"
data5 <- data5 [,-c(13:19)]

df1 <- data1
df2 <- data2
df3 <- data3
df4 <- data4
df5 <- data5
l <- data.frame()
for(i in 1:5){
  df<-get(paste0("df",i))
  l<-rbind(l,df)
}
exp_dat <- l 
exp_dat <- data.frame(l$SNP,exp_dat[,c(11,12,9,10,8,4,5,7)])
colnames(exp_dat) <- c("SNP","exposure","id.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure")

outcome_dat <- read_outcome_data(snps = data$SNP,
                           filename ="outcome.txt",
                           sep = "\t",
                           eaf_col = "Freq",
                           snp_col = "SNP", 
                           beta_col = "BETA",
                           se_col = "SE",
                           effect_allele_col = "A1",
                           other_allele_col = "A2",
                           pval_col = "P"
)

outcome_dat$outcome <- "TS"
outcome_dat$samplesize.outcome <- 500199



mvdat <- mv_harmonise_data(exp_dat, outcome_dat)
res <- mv_multiple(mvdat)
res_OR<-generate_odds_ratios(res$result)
res_OR
write.table(res_OR, file="MV孟德尔随机化结果_two sample.xls",sep="\t",quote=F)

##################MVMR所用的SNP是根据什么规则来的######################
exp1 <- extract_instruments(outcomes="UKB-B-10817")  
exp2 <- extract_instruments(outcomes="ukb-b-18469")  
exp3 <- extract_instruments(outcomes="ieu-a-302")

sum=dim(exp1)[1]+dim(exp2)[1]+dim(exp3)[1]

sum

exp<-rbind(exp_pm2.5,exp_pm10,exp_NO,exp_NO2)  #按列合并
exp<-exp["SNP"]  #提出素有SNP
exp<-unique(exp)  #去非重复SNP
dim(exp)  #计算最后SNP

####################LASSO去除高度共线MR###############
mv_lasso_feature_selection(mvdat)   #LASSO

mv_lasso<-mv_subset(
  mvdat,
  features = mv_lasso_feature_selection(mvdat),   #LASSO后进行MR
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
mv_lasso
mv_lasso_OR<-generate_odds_ratios(mv_lasso$result)
mv_lasso_OR
####################单个SNP的MVMR################
mv_residual(
  mvdat,
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
#######################PRESSO####################
# mr_presso可以完成多变量
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")

SummaryStats<-cbind(mvdat[["exposure_beta"]][,1],
                    mvdat[["exposure_beta"]][,2],
                    mvdat[["exposure_beta"]][,3],
                    mvdat[["exposure_beta"]][,4],
                    mvdat[["exposure_beta"]][,5],
                    mvdat[["exposure_se"]][,1],
                    mvdat[["exposure_se"]][,2],
                    mvdat[["exposure_se"]][,3],
                    mvdat[["exposure_se"]][,4],
                    mvdat[["exposure_se"]][,5],
                    mvdat[["outcome_beta"]],
                    mvdat[["outcome_se"]])
SummaryStats<-data.frame(SummaryStats)
colnames(SummaryStats) <- c("PM2.5_beta","NOx_beta","PM10_beta","noise_beta","NO2_beta","PM2.5_se","NOx_se","PM10_se","noise_se","NO2_se","outcome_beta","outcome_se")
write.table(SummaryStats, file="mvdat_snp.xls",sep="\t",quote=F)
library(MRPRESSO)
mr_presso(BetaOutcome = "outcome_beta",
          BetaExposure = c("PM2.5_beta","NOx_beta","PM10_beta","noise_beta","NO2_beta"), 
          SdOutcome = "outcome_se", 
          SdExposure = c("PM2.5_se","NOx_se","PM10_se","noise_se","NO2_se"),
          OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, 
          data = SummaryStats,
          NbDistribution = 1000, 
          SignifThreshold = 0.05)
#################MendelianRandomization包的几种方法############

#if (!requireNamespace("MendelianRandomization"))
#  install.packages("MendelianRandomization")
library(MendelianRandomization)
#读取实例数据
#准备数据
MRMVInputObject <- mr_mvinput(bx = cbind(ldlc, hdlc, trig),
                              bxse = cbind(ldlcse, hdlcse, trigse),
                              by = chdlodds, 
                              byse = chdloddsse)

MRMVInputObject

MRMVInputObject_1<- mr_mvinput(bx = cbind(SummaryStats$PM2.5_beta,SummaryStats$NOx_beta,SummaryStats$PM10_beta,SummaryStats$noise_beta,SummaryStats$NO2_beta),
                               bxse = cbind(SummaryStats$PM2.5_se,SummaryStats$NOx_se,SummaryStats$PM10_se,SummaryStats$noise_se,SummaryStats$NO2_se),
                               by = SummaryStats$outcome_beta, 
                               byse = SummaryStats$outcome_se)
MRMVInputObject_1
#IVW方法
MRMVObject <- mr_mvivw(MRMVInputObject, 
                       model = "default",
                       correl = FALSE,
                       distribution = "normal",
                       alpha = 0.05)

MRMVObject

MRMVObject <- mr_mvivw(MRMVInputObject_1, 
                       model = "random",
                       correl = FALSE,
                       distribution = "normal",
                       alpha = 0.05)
MRMVObject
#egger方法
MRMVObject<-mr_mvegger(
  MRMVInputObject,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05)
MRMVObject

MRMVObject<-mr_mvegger(
  MRMVInputObject_1,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05)
MRMVObject

#LASSO
MRMVObject<-mr_mvlasso(
  MRMVInputObject,
  orientate = 1,
  distribution = "normal",
  alpha = 0.05,
  lambda = numeric(0)
)
MRMVObject

MRMVObject<-mr_mvlasso(
  MRMVInputObject_1,
  orientate = 1,
  distribution = "normal",
  alpha = 0.05,
  lambda = numeric(0)
)
MRMVObject

#median
MRMVObject<-mr_mvmedian(
  MRMVInputObject,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject

MRMVObject<-mr_mvmedian(
  MRMVInputObject_1,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject


####################RMVMR###########################
#remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
#rawdat_mvmr<-rawdat_rmvmr
head(rawdat_mvmr)

F.data <- format_mvmr(BXGs = SummaryStats[,c(1,2,3,4)],
                      BYG = SummaryStats[,9],
                      seBXGs = SummaryStats[,c(5,6,7,8)],
                      seBYG = SummaryStats[,10],
                      RSID = colnames(SummaryStats))
head(F.data)

sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

res <- ivw_mvmr(r_input = F.data)
write.table(sres, file="F statistic.xls",sep="\t",quote=F)





####################RMVMR###########################
#remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
#rawdat_mvmr<-rawdat_rmvmr
head(rawdat_mvmr)

F.data <- format_mvmr(BXGs = SummaryStats[,c(2,3,4,5)],
                      BYG = SummaryStats[,1],
                      seBXGs = SummaryStats[,c(6,7,8,9)],
                      seBYG = SummaryStats[,10],
                      RSID = rownames(SummaryStats))
head(F.data)

sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

res <- ivw_mvmr(r_input = F.data)












###################本地暴露读取#################
exp_l<-mv_extract_exposures_local(
  c("exp_a.txt","exp_b.txt"),
  sep = " ",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  id_col = "id.exposure"
)
a <- read.table("PTSD.txt",header = T,sep = "\t")
outcome_dat1<-format_data(
  a,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
)



