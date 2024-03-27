#******************************************************************************#

# 日期:2024/03/25
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
options(stringsAsFactors = FALSE) 
options(warn = -1)
argv <- list()


# 参数设置
p <- arg_parser("calculate Volcano")
p <- add_argument(p, "-i", help="input file", type="character",default = "1.基因表达结果.csv")
p <- add_argument(p, "-F", help="log2FoldChange", type="character",default = "1")
p <- add_argument(p, "-P", help="padj", type="character",default = "0.05")
argv <- parse_args(p)

# 加载数据
Dat<-read.csv(argv$i)

# 参数修饰
FF <- argv$F%>%as.numeric()
PP <- argv$P%>%as.numeric()
FF1 = ifelse( argv$F>1,-2,-1)

# 基础绘图
pdf("vol.pdf", width = 8)
ggplot(Dat,aes(x=log2FoldChange,y=-log10(padj)))+
  geom_point(size = 2,
             aes(color = threshold),
             show.legend = F) +
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名x称
  geom_hline(yintercept = -log10(0.05), 
             linetype = 'dotdash', color = 'grey30') +
  geom_vline(xintercept = c(-1, 1), 
             linetype = 'dotdash', color = 'grey30') +
  scale_color_manual(values = c('#1500FF', '#A9A9A9', '#FF0102')) +
  labs(x = 'Log2(fold change)', y = '-log10(p-value)') +
  theme_half_open()
dev.off()


