library(ggplot2)
library(ggpubr)
library("gridExtra")
library(ggsignif)
#setwd('Downloads/down_from36/')
args=commandArgs(T)
f = args[1]
svgf = args[2]
sign_f = args[3]
sign_m = args[4]
cancer_type =  args[5]
#f = 'ENSG00000170323-cg19879471_ENSG00000170323-HNSC-eQTM.txt'
#svgf = 'ENSG00000170323-cg19879471_ENSG00000170323-HNSC-eQTM.svg'
#sign_f='***'
#sign_m='ns'
size=10
d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
d$y <- log10(d$y)
p1 <- ggplot(d, aes(x=x, y=y,color=Group)) +geom_point() + facet_grid(.~Group)
FPKM_max = max(d$y)*1.05
anno <- data.frame(xstar = c(0.5, 0.5), ystar = c(FPKM_max, FPKM_max),lab = c(sign_f, sign_m),Group = c("Female", "Male"))
p1 <- p1 +labs(x='Beta value',title = cancer_type,y="log10(FPKM)")
p1 <- p1 + xlim(0,1.05)
p1 <- p1 + scale_colour_manual(values = c('#DB5289','#6980BD'))
p1 <- p1 +theme_bw() + theme(# axis.ticks.x=element_blank(),
  axis.text.x = element_text(size = size,color = 'black', vjust =0.2, hjust = 0.5, angle =0),
  axis.text.y = element_text(size = size,color = 'black', vjust = 0, hjust =1, angle = 0),
  axis.title.x =  element_text(size = size,color = 'black'),
  axis.title.y =  element_text(size = size,color = 'black'),
  title =  element_text(size = size,color = 'black'),
  plot.title = element_text(hjust = 0.5),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
  axis.line = element_line(size=0.4, colour = "black"),
  panel.border = element_rect(colour ="black",size=0.5),
  legend.position = "none",
  #strip.background = (fill='#DB5289')#c('#DB5289','#6980BD')
)
p1 <- p1 + geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab),color='black') 
svg(svgf,width=4,height=2.5)
print(p1)
dev.off()


