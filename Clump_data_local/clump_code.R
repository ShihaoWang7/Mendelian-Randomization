
#install packages
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

#install.packages("ggplot2")
rm(list = ls())
#library
library(TwoSampleMR)
library(ggplot2)
library(parallel)
library(doParallel)
library(yulab.utils)
library(ieugwasr)
# 查看CPU核心数量
no_cores <- detectCores()
# 指定要使用的CPU核心数量
cl <- makeCluster(10)
# 创建一个注册的CPU集群
registerDoParallel(cl)

setwd("C:\\Users\\wangs\\Desktop\\HCC\\pqtl")

  expo_rt<- read_exposure_data(
    filename = "HIF1A_GRCh37.txt",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FRQ",
    pval_col = "P",
    samplesize_col = "N"
  )
  expo_rt1 <- expo_rt[expo_rt$pval.exposure < 5e-05,] 
  exp <- expo_rt1[,c("SNP","pval.exposure")]
  colnames(exp) <- c("rsid","pval")
  clump_dat <- ld_clump_local(dat = exp, 
                              clump_kb = 5000,
                              clump_r2 = 0.01, 
                              clump_p = 1, 
                              bfile = "C:\\Users\\wangs\\Desktop\\MR\\data_maf0.01_rs_ref\\data_maf0.01_rs_ref", 
                              plink_bin = "C:\\Users\\wangs\\Desktop\\MR\\plink\\plink")
  rownames(expo_rt) <- expo_rt$SNP
  expo_rt <- expo_rt[clump_dat$rsid,]
  write.table(expo_rt,file = "HIF1A_5E-06_10000_0.001.txt",sep = "\t",quote = F,row.names = F)
  




