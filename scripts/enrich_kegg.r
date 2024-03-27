#!/bin/bash

installed_packages <- rownames(installed.packages())
org_db_packages <- installed_packages[grep("^org.*db$", installed_packages)]
sapply(org_db_packages, function(pkg) library(pkg, character.only = TRUE))


suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
argv <- list()
options(show.info = FALSE)
options(stringsAsFactors = F)
options(warn = -1)

# 参数设置
p <- arg_parser("calculate KEGG_rich")
p <- add_argument(p, "-i", help="input file", type="character",default = "DEG.csv")
p <- add_argument(p, "-f", help="input file", type="character",default = "DEG.csv")
argv <- parse_args(p)

# 提取目录部分
data_dir <- dirname(argv$i)
# 设置工作目录
setwd(data_dir)
gene <- read.csv(argv$i)
colnames(gene) <- 'id'

#
dir.create("KEGG富集结果")
org_db_data <- eval(parse(text = org_db_packages[1]))
pathway2gene <- AnnotationDbi::select(org_db_data, keys = keys(org_db_data), columns = c("Pathway", "Ko")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1]%<>% gsub("map","ko",.)
  colnames(keggpathid2name.df) <- c("Pathway","Name")
  return(keggpathid2name.df)
}
pathway2name <- get_path2name()


gene2des <-left_join(pathway2gene,pathway2name) %>% unique() %>%na.omit()


rich <- enricher(gene = gene$id,
                 TERM2GENE = gene2des[c("Pathway","GID")],
                 TERM2NAME = gene2des[c("Pathway","Name")],
                 pvalueCutoff = 1,###原值为0.05
                 pAdjustMethod = 'BH',
                 qvalueCutoff = 1,###原值为0.2
                 maxGSSize = 500)


###取前二十的富集结果，可手动计算
rich_results<-as.data.frame(rich)
rich_results1 <- rich_results[rich_results$p.adjust<0.05 &rich_results$qvalue<0.2,]


rich_20 <- rich_results[1:20,]
richa <- rich_20[c("BgRatio","GeneRatio")]
richa$BgRatio <- sub("/.*","",richa[,1])
richa$GeneRatio <- sub("/.*","",richa[,2])
a1 <- as.numeric(richa$BgRatio)
b1 <- as.numeric(richa$GeneRatio)
rich_20$Rich.factor <- b1/a1


##输出kegg结果


write.csv(rich_results , "./KEGG富集结果/all_result.csv", row.names = FALSE,quote = TRUE)
write.csv(rich_results1, "./KEGG富集结果/rich_result.csv",row.names = FALSE,quote = TRUE)
write.csv(rich_20,"./KEGG富集结果/rich_20.csv",row.names = FALSE,quote = TRUE)


##画图


rich2.csv <- read.csv("./KEGG富集结果/rich_20.csv",encoding = "UTF-8")
labels=(levels(factor(rich2.csv$Description))[as.factor(rich2.csv$Description)])


rich2.csv$number <- factor(rev(1:nrow(rich2.csv)))
names(labels) = rev(1:20)
rich2.csv$shortname<-labels


p <- ggplot(data=rich2.csv, aes(x=number, y=Rich.factor)) +
  geom_point(mapping = aes(size=Count,colour=-log10(qvalue)))+
  coord_flip() + theme_test() +
  scale_color_gradient(low = "darkgreen",high = "red")+
  scale_x_discrete(labels=labels) +
  labs(title = "Top20 of KEGG enrichment",x=" ",y="Rich factor",
       colour="-log10(qvalue)",size="Gene number")+theme_bw()
ggsave("./KEGG富集结果/kegg.pdf", width = 6)


