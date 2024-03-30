#!/usr/bin/env Rscript
#******************************************************************************#

# 日期:2024/03/30
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

argv <- list()
options(warn = -1)

p <- arg_parser("calculate heatmap")
p <- add_argument(p, "-i", help="heatmap_id file", type="character",default = "heatmap_id.txt")
p <- add_argument(p, "-n", help="norm file", type="character",default = "count2tpm.txt")
p <- add_argument(p, "-s", help="group file", type="character",default = "group.txt")

argv <- parse_args(p)


id <- read.csv(argv$i) %>%
  `colnames<-`("Geneid")
matrix <- read_delim(argv$n)

heatmap <- left_join(id,matrix) %>% 
  column_to_rownames(var = "Geneid")

sample_info <- read_delim(argv$s)%>% 
  `colnames<-`(c("id","treatment"))
pdf("heat.pdf", width = 6)
Heatmap(heatmap, name = "Expression",
        #col = col_fun, 
        row_km = 2, 
        column_km = 3,
        column_gap = unit(0.01, 'npc'),
        rect_gp = gpar(col = 'white', lwd = 1),
        column_title = "Gene Expression Heatmap",
        row_title = " ",
        row_names_gp = gpar(fontsize = 10, fontface = 'italic'),
        show_column_names = T,
        top_annotation = HeatmapAnnotation(
          group = sample_info$treatment,col = list(
            group = c(c = '#fc8d59', s = '#99d594'))),
        column_title_gp = gpar(fontsize = 20, fontface = "bold"))
dev.off()

