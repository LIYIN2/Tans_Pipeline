#******************************************************************************#

# 日期:2024/03/25
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#

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
argv <- list()
options(show.info = FALSE)
options(stringsAsFactors = F)
options(warn = -1)

# 参数设置
p <- arg_parser("calculate GO_rich")
p <- add_argument(p, "-i", help="input file", type="character",default = "DEG.csv")
argv <- parse_args(p)

# 提取目录部分
data_dir <- dirname(argv$i)
# 设置工作目录
setwd(data_dir)
gene <- read.csv(argv$i)
colnames(gene) <- 'id'

####go富集
dir.create("GO富集结果")
rich <- enrichGO(gene = gene$id,                              
                 keyType = "GID",                                   
                 OrgDb = org_db_packages,                          
                 ont = "all",                                        
                 pvalueCutoff = 1,                             
                 qvalueCutoff = 1,                                 
                 pAdjustMethod = "BH",
                 readable = FALSE)   
###取前二十的富集结果，求generatio这一步调用do包，也可手动计算
rich_results<-as.data.frame(rich)
rich_results <- rich_results[order(rich_results$qvalue),]
rich_results1 <- rich_results[rich_results$p.adjust<0.05 &rich_results$qvalue<0.2,]
rich_20 <- rich_results[1:20,]


richa <- rich_20[c("BgRatio","GeneRatio")] 
richa$BgRatio <- sub("/.*","",richa[,1])
richa$GeneRatio <- sub("/.*","",richa[,2])


a1 <- as.numeric(richa$BgRatio)
b1 <- as.numeric(richa$GeneRatio)
rich_20$Rich.factor <- b1/a1






##输出go结果


write.csv(rich_results , "./GO富集结果/all_result.csv", row.names = FALSE,quote = TRUE)
write.csv(rich_results1, "./GO富集结果/rich_result.csv",row.names = FALSE,quote = TRUE)
write.csv(rich_20,"./GO富集结果/rich_20.csv",row.names = FALSE,quote = TRUE)


##画图


rich2.csv <- read.csv("./GO富集结果/rich_20.csv",encoding = "UTF-8")
labels=(levels(factor(rich2.csv$Description))[as.factor(rich2.csv$Description)])


rich2.csv$number <- factor(rev(1:nrow(rich2.csv)))
names(labels) = rev(1:20)
rich2.csv$shortname<-labels




p <- ggplot(data=rich2.csv, aes(x=number, y=Rich.factor)) +
  geom_point(mapping = aes(size=Count,colour=-log10(qvalue),shape=ONTOLOGY))+
  coord_flip() + theme_test() +
  scale_color_gradient(low = "darkgreen",high = "red")+
  scale_x_discrete(labels=labels) +
  labs(title = "Top20 of GO enrichment",x=" ",y="Rich factor",
       colour="-log10(qvalue)",size="Gene number")+theme_bw()
ggsave("./GO富集结果/go.pdf", width = 6)

