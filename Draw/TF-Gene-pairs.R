#setwd('/Users/yangmengyuan/Documents/All_Project/SexAnnoDB/network_with_threshold/')
library(tidyverse)

args=commandArgs(T)
outpath = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFigures/figures/'
f = args[1]
geneid = args[2]
inpath = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFigures/files'
#f = 'ENSG00000163220_TF-Gene.txt'
inf = str_c(inpath,f,sep='/')  #'/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFigures/files/ENSG00000163220_TF-Gene.txt'
#geneid = 'ENSG00000163220'

d <- read.table(inf,header=TRUE,sep='\t')
n <- length((unique(d$cancer)))
m = length((unique(d$TF)))
p <- ggplot(d, aes(x = Score, y = TF,color=NetName)) +
  geom_segment(aes(x = 0, y = TF, xend = Score, yend =TF )) +
  geom_point( aes(size = Threshold)) + scale_colour_manual(values = c('#DB5289', '#6980BD'))
p <- p + theme_bw() 
p <- p +  facet_grid ( ~cancer)
p
outf = str_c(outpath,geneid,'_TF-Gene.svg',sep='')
svg(outf,n*1.5+2,m*0.2+1)
print(p)
dev.off()

