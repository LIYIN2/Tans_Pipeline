#!/usr/bin/env Rscript
#******************************************************************************#

# 日期:2024/03/25
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
argv <- list()
options(warn = -1)
# 参数设置
p <- arg_parser("calculate correlation")
p <- add_argument(p, "-i", help="input file", type="character",default = "count2tpm.txt")
p <- add_argument(p, "-s", help="sample file", type="character",default = "sample.txt")
argv <- parse_args(p)

group <- read.table(argv$s,header = T,row.names = 1)
sample <- row.names(group)
count2tpm <- read_delim(argv$i) %>% dplyr::select(sample)

COR <- cor(count2tpm)
print("1")

## 绘制热图
pdf("cor.pdf", width = 8)
pheatmap(COR,display_numbers = T,
         annotation_col = group,
         fontsize = 10,cellheight = 20,
         cellwidth = 20,cluster_rows = T,
         cluster_cols = T)
dev.off()
