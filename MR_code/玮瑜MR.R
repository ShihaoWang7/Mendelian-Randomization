#TSMR分析

#对于暴露和结局的数据都来自Open GWAS,即存在暴露ID号
#以BMI和CHD为例：BMI的ID：ieu-a-2  CHD的ID：ieu-a-7
#首先在R里面安装TSMR包
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")

#运行TSMR
library(TwoSampleMR)
#提取暴露BMI的工具变量
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')

#解释extract_instruments函数
#extract_instruments(outcomes,p1 = 5e-08,clump = TRUE,
p2 = 5e-08,r2 = 0.001,kb = 10000,access_token = ieugwasr::check_access_token(),
force_server = FALSE)

#如果想要调整clump参数
bmi <- extract_instruments(outcomes = 'ieu-a-2',
                           clump = TRUE, r2 = 0.01,
                           kb = 5000, access_token = NULL)

#如果想要调整P值
bmi_1 <- extract_instruments(outcomes = "",
                             p1 = 5e-06,
                             clump = TRUE, r2 = 0.001,
                             kb = 10000, access_token = NULL)


#提取工具变量在结局中的信息
chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, outcomes = 'ieu-a-7')

#将暴露和结局的数据进行合并，产生用于进行MR分析的数据
#第一种代码：
dat <- harmonise_data(bmi_exp_dat, chd_out_dat)
#第二种代码：
dat <- harmonise_data(
  exposure_dat=bmi_exp_dat,
  outcome_dat=chd_out_dat,
  action= 2
)

#MR分析的主要结果:默认用5种方法进行MR分析
res <- mr(dat)
res

#换算成OR值
OR <-generate_odds_ratios(res)
OR

#如果MR分析中限定方法，如只用mr_egger和mr_ivw
mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))

#使用随机效应模型
RE <-mr(dat,method_list=c('mr_ivw_mre'))
REOR <-generate_odds_ratios(RE)

#固定效应模型
FE <-mr(dat,method_list=c('mr_ivw_fe'))
FEOR <-generate_odds_ratios(FE)

#离群值检验
#安装MRPRESSO包
devtools::install_github("rondolab/MR-PRESSO",force = TRUE)

#安装完运行
library(MRPRESSO)
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
          SignifThreshold = 0.05)

#单个snp分析-（三种方法）
#默认Wald比值
res_single <- mr_singlesnp(dat)
ORR <-generate_odds_ratios(res_single)

#数据敏感性分析
#异质性检验
het <- mr_heterogeneity(dat)
het
#多效性检验
pleio <- mr_pleiotropy_test(dat)
pleio
#逐个剔除检验
single <- mr_leaveoneout(dat)
mr_leaveoneout_plot(single)
#散点图
mr_scatter_plot(res,dat)
#森林图
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
#漏斗图
mr_funnel_plot(res_single)


#保存数据
#安装xlsx
install.packages('xlsx')
#运行xlsx
library(xlsx)
write.xlsx(ORR, "D:xx.xls")

#本地为暴露文件，结局为 CHD的ID：ieu-a-7
#运行TSMR
library(TwoSampleMR)
#读取暴露本地数据
exp_dat <- read_exposure_data(
  filename = 'Blood selenium.csv',
  sep= ",",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col ="EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "P"
)
exp_dat$exposure <- "Blood selenium"
#读取工具变量在结局当中的信息
CHD_out <- extract_outcome_data(
  snps=exp_dat$SNP,
  outcomes='ieu-a-7',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL
)
mydata <- harmonise_data(
  exposure_dat=exp_dat,
  outcome_dat=CHD_out,
  action= 2
)
res <- mr(mydata)
res
OR <-generate_odds_ratios(res)
OR
#异质性
het <- mr_heterogeneity(mydata)
het
#多效性
pleio <- mr_pleiotropy_test(mydata)
pleio
#逐个剔除检验
single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
#散点图
mr_scatter_plot(res,mydata)
#森林图
res_single <- mr_singlesnp(mydata)
mr_forest_plot(res_single)
#漏斗图
mr_funnel_plot(res_single)

#结局为本地文件，怎么读取
outcome_dat<-read_outcome_data(snps = as_exp_dat$SNP,
                               filename = "LBD.GWAS.Summary.Stats.EBI.submission.formatted.tsv",
                               sep = "\t",
                               eaf_col = "effect_allele_frequency",
                               snp_col = "variant_id", 
                               beta_col = "beta",
                               se_col = "standard_error",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = "p_value"
)

outcome_dat$outcome <- "LBD"
