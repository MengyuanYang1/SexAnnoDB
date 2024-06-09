library(oce)
library(plyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggdendro)
library(ComplexHeatmap)

args=commandArgs(T)
f1 = args[1];outsvg = args[2]
#f1 = '/Users/yangmengyuan/Documents/All_Project/SexDiff/SexDiffinCancer/Scripts/ENSG00000278318.txt'
#outsvg = '/Users/yangmengyuan/Documents/All_Project/SexDiff/SexDiffinCancer/Scripts/ENSG00000278318_heatmap.svg'

d <- read.table(f1 ,header =TRUE, stringsAsFactors = FALSE,row.names= 1,sep=',');
d.m <-melt(d)
d.m <- ddply(d.m, .(variable), transform,rescale = rescale(value))
id_list = c(sort(rownames(d),decreasing=T))
d.m$ID <-  factor(rownames(d), id_list )
p <- ggplot(d.m, aes(variable,ID)) + geom_tile(aes(fill = value),colour = "white")
p <- p + scale_fill_gradient(low = "white",high = "#3399FF")  ###"#3399FF" blue, "#9933FF" purple
p <- p + xlab('')+ ylab('')+labs(title =''  )## +xlim(0,1)
p <- p +theme(axis.text.x = element_text(size = 12,color = 'black', vjust = 0.2, hjust =0, angle =90),
      axis.text.y = element_text(size = 12,color = 'black',vjust=0.2,hjust=0.95,angle = 0),
      axis.title.x =  element_text(size = 12,color = 'black'),
      axis.title.y =  element_text(size = 12,color = 'black'),
      title =  element_text(size = 12,color = 'black'),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major =element_blank(),
      panel.grid.minor = element_blank(),
)
raws = dim(d)[1]
if(raws <= 2) {w = 2.3} else {w = (raws - 2)*0.4 +3}
svg(outsvg,16,w)
p
dev.off()


