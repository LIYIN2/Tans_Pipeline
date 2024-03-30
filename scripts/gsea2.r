#!/usr/bin/env Rscript
#******************************************************************************#

# 日期:2024/03/30
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse))
argv <- list()
options(warn = -1)

p <- arg_parser("calculate GSEA")
p <- add_argument(p, "-i", help="input file", type="character",default = "GSEA.RData")
p <- add_argument(p, "-n", help="geneset name", type="character",default = "molecular transducer activity")
argv <- parse_args(p)

print("111")
gset="molecular transducer activity"
outname <- paste0(gset,'.','pdf')
load(argv$i)
gseaplot2(GSEA_enrichment,gset,
          title = gset,
          color="red",
          base_size = 15,
          pvalue_table = 0
)

ggsave(outname, width = 10)
