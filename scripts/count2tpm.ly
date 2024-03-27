#!/usr/bin/env Rscript

#******************************************************************************#

# 日期:2024/03/25
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(argparser))
argv <- list()
# 参数设置
p <- arg_parser("calculate TPM")
p <- add_argument(p, "-i", help="input file", type="character",default = "featureCounts.txt")
#p <- add_argument(p, "-o", help="output file", type="character",default = "count2tpm.txt")

# 参数解析
argv <- parse_args(p)
ip <- argv$i

df <- read.table(ip,header = T,check.names=FALSE)
# 自定义函数：来源https://zhuanlan.zhihu.com/p/513391213
countToTpm <- function(counts, effLen)
{
rate <- log(counts) - log(effLen)
denom <- log(sum(exp(rate)))
exp(rate - denom + log(1e6))
}

# 计算
a <- ncol(df)
for (i in 7:a){
  df[,i] <- countToTpm(df[,i],df[,6])  
}
# 输出
write.table(df,"count2norm.txt",row.names = F,quote = FALSE,sep = '\t')
