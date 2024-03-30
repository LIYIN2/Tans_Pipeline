#!/usr/bin/env Rscript
#******************************************************************************#

# 日期:2024/03/30
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(clusterProfiler)) 
suppressPackageStartupMessages(library(GSEABase)) 
suppressPackageStartupMessages(library(enrichplot))
argv <- list()
options(warn = -1)


# 参数设置
p <- arg_parser("calculate GSEA")
p <- add_argument(p, "-i", help="input file", type="character",default = "geneSet.csv")
p <- add_argument(p, "-e", help="expression file", type="character",default = "1.基因表达结果.csv")
argv <- parse_args(p)

exp <- read.csv(argv$e) %>% 
  dplyr::select(c("id","log2FoldChange"))
geneSet <- read.csv(argv$i)


geneList <- exp$log2FoldChange #获取GeneList
names(geneList) <- exp$id
geneList <- sort(geneList, decreasing = T) #从高到低排序


GSEA_enrichment <- GSEA(geneList, 
                        TERM2GENE = geneSet, 
                        pvalueCutoff = 0.05, 
                        minGSSize = 10, 
                        maxGSSize = 500, 
                        eps = 0, 
                        pAdjustMethod = "BH") 
result <- data.frame(GSEA_enrichment)
dir.create("GSEA富集结果")

##显示最显著的15个通路
dotplot(GSEA_enrichment,showCategory=15,
        color="p.adjust")
ggsave("./GSEA富集结果/GSEA_top15.pdf", width = 6)
write.csv(result,"./GSEA富集结果/GSEA_result.csv", row.names = FALSE,quote = TRUE)
##按照上调和下调展示通路
dotplot(GSEA_enrichment,split = ".sign")+facet_grid(~.sign)+
  theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
        axis.title = element_text(size = 10,color ="black"), 
        axis.text = element_text(size= 10,color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        #legend.position = "top",##标签位置
        legend.text = element_text(size= 10),
        legend.title= element_text(size= 10))
ggsave("./GSEA富集结果/GSEA_up-down.pdf", width = 6)
save(GSEA_enrichment, file = "./GSEA富集结果/GSEA.RData",compress="xz")

