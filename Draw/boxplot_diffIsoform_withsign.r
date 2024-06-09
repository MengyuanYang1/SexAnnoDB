library(ggplot2)
library(ggpubr)
library("gridExtra")


args=commandArgs(T)
f = args[1]
svgf = args[2]
SE_id = args[3]

boxplot_for_isoform <- function(f,svgf){
	d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
	d <- data.frame(d)
	d$label <- paste(d$isoformid,d$tissue,sep='\n')
	d$label <- as.factor(d$label)
	y_position = max(d$TPM)*0.9
	d_tmp = d[,c('label','sign')]
	d_tmp<-d_tmp[!duplicated(d_tmp, fromLast=TRUE),] 
	d_tmp$sign_x <- seq(from=0.775, to=dim(d_tmp)[1], by=1)
	d_tmp$sign_xend <- seq(from=1.125, to=dim(d_tmp)[1]+1, by=1)
	d_tmp$y_position <- rep(y_position, times=dim(d_tmp)[1])
	labels = levels(factor(d$label))
	d$label <- factor(d$label,,levels=labels)
	d_tmp$label <- factor(d_tmp$label,levels=labels)
	print(d_tmp)
	n = length(labels)
	x = 'label'
	y = 'TPM'
	group = 'diagnosis'
	x_lab = ''
	p1 <- drawfigure(d,d_tmp,x_lab,x,y,group)
	if (n<2){wid = 3}
	else{wid =3+(n-2)*0.3}
	svg(svgf,width=wid,height=3)
	print(p1)
	dev.off()
}


boxplot_for_SE <- function(f,svgf){
  d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
  ensts = levels(factor(d$isoformid))
  n = length(ensts)
  d$tissue = sub("-","\n", d$tissue)
  plot = list()
  for (i in 1:length(ensts)){
	enst = ensts[i]
	print(enst)
	d_sub <- d[which(d$isoformid == enst),]
	print(dim(d_sub))
	plot[[i]]  <- drawfigure(d_sub,enst)
  }
  print(length(plot))
  glist <- lapply(plot, ggplotGrob)
  pngf <- sub("svg","png",svgf)
  ggsave(pngf,width =4, height = 3.5*length(plot),marrangeGrob(glist, nrow = length(plot), ncol = 1,list(top=NULL)))
}

drawfigure <- function(d,d_sign,enst,x,y,group){
	#x = label,y = TPM, group = diagnosis
  #symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01,1), symbols = c("****", "***", "**", "*"))
  p1 <- ggboxplot(d,x,y,color=group,palette = "npg") # add.params = list(fill = "diagnosis"))
  #p1 <- ggplot(d,mapping = aes(x =tissue_p, y =psi)) + geom_violin(color="white") + geom_boxplot(width=0.2,position=position_dodge(0.9))
  #p1 <- p+stat_compare_means(label.x = 1.5,hide.ns = TRUE, method = "t.test",aes(group =diagnosis),symnum.args = symnum.args,label = "p.signif",label.x.npc='center',label.y.npc= 'top') ##label = "p.signif",
  p1 <- p1 +labs(x=enst,y="TPM") #,tittle = "")
  #p1 <- p1 +  facet_grid(isoformid ~ .)
  print(d_sign$sign_x)
  print(d_sign$sign_xend)
  p1 <- p1 + geom_signif(xmin=d_sign$sign_x,xmax=d_sign$sign_xend, y_position=d_sign$y_position, annotation=d_sign$sign)
  size = 8
  p1 <- p1 +theme_bw() + theme(# axis.ticks.x=element_blank(),
                 axis.text.x = element_text(size = size,color = 'black', vjust =0.2, hjust = 0.5, angle =90),
                 axis.text.y = element_text(size = size,color = 'black', vjust = 0, hjust =1, angle = 0),
                 axis.title.x =  element_text(size = size,color = 'black'),
                 axis.title.y =  element_text(size = size,color = 'black'),
                 title =  element_text(size = size,color = 'black'),
                 plot.title = element_text(hjust = 0.5),
		 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
		 axis.line = element_line(size=0.4, colour = "black"),
		 panel.border = element_rect(colour ="black",size=0.5)
		 #legend.position = "none"
		)
return (p1)
  }


boxplot_for_isoform(f,svgf)
