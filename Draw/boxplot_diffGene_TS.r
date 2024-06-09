library(ggplot2)
library(ggpubr)
library(ggsignif)
library("gridExtra")

args=commandArgs(T)
f = args[1]
svgf = args[2]

#f = 'ENSG00000004939.txt'
#svgf = 'ENSG00000004939.svg'

boxplot_for_isoform <- function(f,svgf){
	d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
	d <- data.frame(d)
	y_position = max(d$FPKM)*0.9
	d_tmp = d[,c('cancer_type','sign')]
	d_tmp<-d_tmp[!duplicated(d_tmp, fromLast=TRUE),] 
	d_tmp$sign_x <- seq(from=0.775, to=dim(d_tmp)[1], by=1)
	d_tmp$sign_xend <- seq(from=1.125, to=dim(d_tmp)[1]+1, by=1)
	d_tmp$y_position <- rep(y_position, times=dim(d_tmp)[1])
	cancer_types = levels(factor(d$cancer_type))
	d$cancer_type <- factor(d$cancer_type,levels=cancer_types)
	d_tmp$cancer_type <- factor(d_tmp$cancer_type,levels=cancer_types)
	n = length(cancer_types)
	x = 'cancer_type';y = 'FPKM';group = 'group';x_lab = '';
	p1 <- drawfigure(d,d_tmp,x_lab,x,y,group)
	if (n<2){wid = 3}else{wid =3+(n-2)*0.3}
	svg(svgf,width=wid,height=3)
	print(p1)
	dev.off()
}
drawfigure <- function(d,d_sign,enst,x,y,group){
  p1 <- ggboxplot(d,x,y,color =group) 
  p1 <- p1 +labs(x=enst,y="FPKM")
  p1 <- p1 + scale_colour_manual(values = c('#DB5289','#6980BD'))
  #p1 <- p1 + scale_colour_manual(values = c('#DB5289','black'))
  p1 <- p1 + geom_signif(xmin=d_sign$sign_x,xmax=d_sign$sign_xend, y_position=d_sign$y_position, annotation=d_sign$sign)
  size = 8
  p1 <- p1 +theme_bw() + theme(# axis.ticks.x=element_blank(),
     axis.text.x = element_text(size = size,color = 'black', vjust =0.2, hjust = 0.5, angle =90),
     axis.text.y = element_text(size = size,color = 'black', vjust = 0, hjust =1, angle = 0),
     axis.title.x =  element_text(size = size,color = 'black'),
     axis.title.y =  element_text(size = size,color = 'black'),
     title =  element_text(size = size,color = 'black'),
     plot.title = element_text(hjust = 0.5),
		 panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
		 axis.line = element_line(size=0.4, colour = "black"),
		 panel.border = element_rect(colour ="black",size=0.5)#legend.position = "none"
		 )
return (p1)
  }

boxplot_for_isoform(f,svgf)
