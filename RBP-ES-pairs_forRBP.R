library(tidyverse)
library(tidygraph)
library(ggraph)
library(cowplot)
library(conflicted)  
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

args=commandArgs(T)
f = args[1]
geneid = args[2]

inpath = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFiguresRBP/files/'
outpath = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFiguresRBP/figures/'

inf = str_c(inpath,f,sep='/')
outf = str_c(outpath,geneid,'_RBP-ES_forRBP.svg',sep='')

d <- read.table(inf,header=TRUE,sep='\t')
n_can <- length(unique(d$Cancer))
if (n_can < 6){n_row = 1}else{n_row = ceiling(n_can/6)}
cancers <- unique(d$Cancer)
my_Plot <- function(d_graph){
  ggraph(d_graph, layout = "stress") +
    geom_edge_diagonal(aes(colour = factor(Type)),
        arrow = arrow(length = unit(3, 'mm')),start_cap = circle(3, "mm"),end_cap = circle(3, "mm"),show.legend = T)+
    scale_edge_color_manual(values = c('#DB5289', '#6980BD'))+
    geom_node_point(size=5,alpha = 0.3,shape=24,color='black',fill='black')+
    geom_node_text(aes(label = name),repel = TRUE,size=2.8) +
    theme(axis.text.x =element_blank(),legend.position = "none")+
  facet_wrap(~Cancer)
}

my_graph <-  function(cancer){ as_tbl_graph(subset(d[which(d$Cancer == cancer), ][1:20,], !is.na(from)), directed = FALSE) }
g_list <- lapply(cancers, my_graph)
p_list <- lapply(g_list, my_Plot)
svg(outf,ceiling(n_can/n_row)*3.6,n_row*3.6)
plot_grid(plotlist = p_list, align = "h", nrow = n_row)
dev.off()


