#!/usr/bin/env Rscript
#******************************************************************************#

# 日期:2024/03/30
# 作者:liy<liyin59375@gmail.com>

#******************************************************************************#
suppressPackageStartupMessages(library(cowplot))
#suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(PCAtools))
argv <- list()
options(warn = -1)

p <- arg_parser("calculate PCA")
p <- add_argument(p, "-i", help="norm file", type="character",default = "count2tpm.txt")
p <- add_argument(p, "-s", help="group file", type="character",default = "group.txt")
argv <- parse_args(p)

#加载数据
gene_exp <- read.table(argv$i, 
                       header = T, row.names = 1,check.names = F)

sample_info <- read.table(argv$s,header = T,row.names = 1)
pca <- pca(gene_exp, metadata = sample_info)
pca_prcp_contrib <- pca$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:2] %>% signif(digits = 3)

##加入分组信息
rotated <- pca$rotated
a <-  rownames_to_column(rotated,var="sample")
b <- rownames_to_column(sample_info,var="sample")
pca_rotated_plus <- left_join(a,b,by="sample")

ggplot(pca_rotated_plus, aes(x = PC1, y = PC2)) +
  geom_point(size = 8, aes(fill = sample, shape = group)) +
  stat_ellipse(aes(color = sample)) +
  labs(x=paste0('PC1','(',pca_prcp_contrib[1],'%',')'),
       y=paste0('PC2','(',pca_prcp_contrib[2],'%',')'))+
  scale_shape_manual(values = 21:24) +
  #scale_fill_brewer(palette = 'Set2') + 
  #scale_color_brewer(palette = 'Set2') + 
  theme_classic()+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  theme(
    #legend.position = 'top',  # 调整图例位置为顶部
    #legend.justification = c(0,1),  # 图例位置左上角
    axis.title = element_text(face='bold',size=13)  # x,y轴标题加粗，字体大小为13
  )

ggsave("PCA.pdf", width = 10)

