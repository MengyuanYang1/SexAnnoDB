#setwd('/Users/yangmengyuan/Documents/All_Project/SexAnnoDB/network_with_threshold/')
library(tidyverse)

args=commandArgs(T)
outpath = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFigures/figures/'
f = args[1]
geneid = args[2]
inpath = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFigures/files'
#f = 'ENSG00000163220_TF-Gene.txt'
inf = str_c(inpath,f,sep='/')  #'/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFigures/files/ENSG00000163220_TF-Gene.txt'
#geneid = 'ENSG00000163220'

d <- read.table(inf,header=TRUE,sep='\t')
n <- length((unique(d$Cancer)))
m = length((unique(d$RBP)))
l = length((unique(d$ESID)))

p <- ggplot(d, aes(x = Score, y = RBP,color=NetName,fill=NetName)) +
  geom_segment(aes(x = 0, y = RBP, xend = Score, yend =RBP )) +
  geom_point( aes(size = Threshold),shape=24) + scale_colour_manual(values = c('#DB5289', '#6980BD'))+ scale_fill_manual(values = c('#DB5289', '#6980BD'))
p <- p + theme_bw() 
p <- p +  facet_grid ( ESID~Cancer)
p
outf = str_c(outpath,geneid,'_RBP-ES.svg',sep='')
svg(outf,n*1.5+2,m*0.6*l+1)
print(p)
dev.off()

