#!/usr/bin/env Rscript
#******************************************************************************#

# 日期:2024/03/25
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors = FALSE) 
argv <- list()
options(warn = -1)

# 参数修改,往后不用修改
#argv$matrix <- "featureCounts_clean.txt" ##表达矩阵
#argv$sample <- "sample.txt" ##样品分组
#argv$s <- '1' ##输入1或-1选择前/后-对照组/实验组,默认1,即前为实验组

# 参数设置
p <- arg_parser("calculate DEG")
p <- add_argument(p, "-i", help="input file", type="character",default = "featureCounts_clean.txt")
p <- add_argument(p, "-s", help="sample file", type="character",default = "sample.txt")
p <- add_argument(p, "-F", help="log2FoldChange", type="character",default = "1")
p <- add_argument(p, "-P", help="padj", type="character",default = "0.05")

#p <- add_argument(p, "-t", help="sample file", type="character",default = "1")
# 参数解析
argv <- parse_args(p)
argv$matrix <- argv$i
argv$sample <- argv$s
argv$s <- "1"
# 数据读取
matrix <- read.csv(argv$matrix,sep = "\t",
                   row.names = 1,check.names = F) 
sample <- read.table(argv$sample,header = T) %>% 
  `colnames<-`(c('sample', 'treatment'))
data <- select(matrix,sample[,1])

# 构建差异表达表

## Step1：构建DEseq2对象(dds)
dds <- DESeqDataSetFromMatrix(data, 
                              colData=sample, design= ~ treatment)
## Step2：标准化数据
dds <- DESeq(dds)##计算并输出

## Step3：从DESeq矩阵中提取差异表达表
nn <- unique(sample$treatment)%>%as.data.frame()
if(argv$s==1)
{res <- results(dds, contrast = c("treatment",nn[1,1],nn[2,1]))
suppressMessages(sig<-print(paste(nn[1,1],"vs",nn[2,1])))
res$signal<-sig
}else
{res <- results(dds, contrast = c("treatment",nn[2,1],nn[1,1]))
suppressMessages(sig<-print(paste(nn[2,1],"vs",nn[1,1])))
res$signal<-sig
}

## Step4：将差异表达表转化数据框,并去除缺失值
DEG_DESeq2 <- as.data.frame(res) %>%
  na.omit()

# 输出差异表达表
DEG_DESeq2 <- rownames_to_column(DEG_DESeq2,var="id") 
write.csv(DEG_DESeq2,"result.csv",row.names = F,quote = T)

# 加载 R 包
library(tidyverse)
options(stringsAsFactors = FALSE)

# 参数修改,往后不用修改
argv$input_file <- 'result.csv'##差异基因文件名,第一列为id
argv$result_DET <- 'DEG.csv'##输出文件名


# 加载数据
Dat<-read.csv(argv$input_file)

## Step1 添加上下调信号
FF <- argv$F%>%as.numeric()
PP <- argv$P%>%as.numeric()
Dat$threshold = factor(ifelse(Dat$padj < PP & abs(Dat$log2FoldChange) >= FF, 
                       ifelse(Dat$log2FoldChange>= FF ,'Up','Down'),'NoSignifi'),
		       levels=c('Up','Down','NoSignifi'))
## Step2 统计频数
det.stat <- list()
det.stat$'---------cutoff-------' <- '' 
det.stat$log2FC<-FF
det.stat$padj<-PP
B <- Dat$threshold%>% table()%>%as.data.frame()
det.stat$'---------------------' <- '' 
det.stat$Up<- B[1,2]
det.stat$Down<- B[2,2]
det.stat$NoSignifi <- B[3,2]
det.stat$NoSignifi <- B[3,2]
det.stat <- t(t(det.stat))


gene <- filter(Dat,
               abs(log2FoldChange)>FF&padj<PP)%>%
  pull(id)%>%as.data.frame()
colnames(gene)<-'id'
Dat2 <- dplyr::select(Dat,id,log2FoldChange,pvalue,padj,threshold,signal)


## 输出结果
write.csv(gene,file='DEG.csv',row.names = F)
write.table(det.stat,file='2.基因表达结果统计.txt',row.names = T,col.names = F)
write.csv(Dat2,file='1.基因表达结果.csv',row.names = F)
