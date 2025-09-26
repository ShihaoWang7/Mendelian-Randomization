####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com

setwd("/Users/shengxinkuangren/Desktop/BWMR/No1/ebi-a-GCST90001943")

library(ggplot2)
source("BWMR_updated.R") 

harm_rt=read.table("harmonise.txt",header = T,sep = "\t")

myBWMR <- BWMR(gammahat = harm_rt$beta.exposure,
                 Gammahat = harm_rt$beta.outcome,
                 sigmaX = harm_rt$se.exposure,
                 sigmaY = harm_rt$se.outcome) 
####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com
beta=myBWMR[["beta"]]
lci95 <- myBWMR[["beta"]]-1.96*myBWMR[["se_beta"]]
uci95 <- myBWMR[["beta"]]+1.96*myBWMR[["se_beta"]]
or <- exp(myBWMR[["beta"]])
or_lci95 <- exp(lci95)
or_uci95 <- exp(uci95)
pval<-myBWMR[["P_value"]]

myres=data.frame(method="BWMR",beta,lci95,uci95,or,or_lci95,or_uci95,pval)
myres$se=abs(log(or)/qnorm(pval/2))
write.table(myres,"BWMR_OR.txt",row.names = F,sep = "\t",quote = F)

####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com
