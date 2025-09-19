
#install packages
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")

#install.packages("ggplot2")
rm(list = ls())
#library
library(TwoSampleMR)
library(ggplot2)
library(foreach)
library(ggpubr)
setwd("C:\\Users\\wsh20\\Desktop")
iddf=read.table("id.txt",header =T,sep = "\t")


dxw=as.vector(iddf$id)

  expo_rt<- read_exposure_data(
    filename = "FAM83D.txt",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "Freq",
    pval_col = "p"
  )
 expo_rt$exposure <- "Idiopathic urticaria" 
foreach(i=dxw, .errorhandling = "pass") %do%{
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = paste0("危险/",i,".tsv"),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    samplesize_col = "n"
  )
  
  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf<- mean( harm_rt$f)
  harm_rt<-harm_rt[harm_rt$f>10,]
  harm_rt=steiger_filtering(harm_rt)
  harm_rt<-harm_rt[harm_rt$steiger_dir == TRUE,]
  harm_rt <- harm_rt[harm_rt$mr_keep == TRUE,]

  harm_rt$outcome <- i
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result) 
  dir.create("CXCL1") 
  filename="CXCL1"
  write.table(harm_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
  write.table(result_or[,3:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
  pleiotropy=mr_pleiotropy_test(harm_rt)
  write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
  heterogeneity=mr_heterogeneity(harm_rt)
  write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
  presso=run_mr_presso(harm_rt)
  capture.output(presso,file = paste0(filename,"/presso.txt"))
  
  p1 <- mr_scatter_plot(mr_result, harm_rt)
  ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
  
  
  
  singlesnp_res<- mr_singlesnp(harm_rt)
  singlesnpOR=generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
  sen_res<- mr_leaveoneout(harm_rt)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
  res_single <- mr_singlesnp(harm_rt)
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
  p5 <- ggarrange(p3[[1]], p1[[1]], p2[[1]], p4[[1]], ncol = 2, nrow = 2,
            labels = c("A","B","C","D"), # 添加标签
            font.label = list(size = 14, face = "bold")) # 设置标签字体样式
  ggsave(p5, file=paste0(filename,"/mergeplot.pdf"), width=10, height=10)
}                    
  
  
  rm(list = ls())