#******************************************************************************#

# 日期:2024/03/25
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(argparser))
options(stringsAsFactors = F)
argv <- list()
options(warn = -1)

#参数设置
p <- arg_parser("calculate org.db")
p <- add_argument(p, "-i", help="emapper file", type="character",default = "out.emapper.annotations.txt")
##p <- add_argument(p, "-s", help="genus_id", type="character",default = "NULL")
##p <- add_argument(p, "-t", help="species_id", type="character",default = "NULL")

##参数解析
argv <- parse_args(p)

#数据处理
emapper <- read_delim(argv$i,
                      "\t", escape_double = FALSE, trim_ws = TRUE)
emapper[emapper==""]<-NA
gene_info <- emapper %>% dplyr::select(GID = query, GENENAME = Preferred_name) %>% na.omit()
##提取GO信息
gene2go <- emapper %>% dplyr::select(GID = query, GO = GOs) %>%
  separate_rows(GO, sep = ",")
gene2go$GO[gene2go$GO=="-"]<-NA
gene2go<-na.omit(gene2go)
gene2go$EVIDENCE <- "IEA"

##提取KEGG信息
gene2ko <- emapper %>% dplyr::select(GID = query, Ko = KEGG_ko) %>% na.omit()
gene2pathway <- emapper %>% dplyr::select(GID = query, Pathway = BRITE) %>%
  separate_rows(Pathway, sep = ",")
gene2pathway$Pathway[gene2pathway$Pathway=="-"]<-NA
gene2pathway<-na.omit(gene2pathway)

#get_path2name <- function(){
#  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
#  keggpathid2name.df[,1]%<>% gsub("map","ko",.)
#  colnames(keggpathid2name.df) <- c("Pathway","Name")
#  return(keggpathid2name.df)
#}
#pathway2name <- get_path2name()
#gene2pathway <- gene2Pathway %>% left_join(pathway2name, by = "Pathway") %>%
#  dplyr::select(GID, Pathway) %>% na.omit()

##org 库构建
makeOrgPackage(gene_info=gene_info,
                 go=gene2go,
                 ko=gene2ko,
                 pathway=gene2pathway,
                 version="0.0.1",
                 maintainer = "my <email@example.com>",
                 author = "my",
                 outputDir = ".",  #输出文件位置
                 tax_id=1,
                 genus=NULL,
                 species=NULL,
                 goTable="go")
#安装org.db
files <- list.files(pattern = "^org")
lapply(files, function(pkg) install.packages(pkg, repos = NULL, type = "sources"))
cat("恭喜构建成功！\n")
