library(ggplot2)
library(ggpubr)
library("gridExtra")
library(ggsignif)
library(getopt)

#spec <- matrix(
#  c("inf",  "i", 2, "character", "This is inf",
#   "svgf", "o", 2, "character",  "This is outf!",
#   "sign_f",  "f", 2, "character",  "This is sign_f!",
#   "sign_m",   "m", 2, "character",  "This is sign_m!",
#   "cancer_type","c",2, "character", "This is cancer type!"),
#  byrow=TRUE, ncol=5)

#print(spec)
#opt <- getopt(spec=spec)

#f = opt$inf
#svgf = opt$svgf
#sign_f = opt$sign_f
#sign_m = opt$sign_m
#cancer_type =  opt$cancer_type
args=commandArgs(T)
f = args[1]
svgf = args[2]
sign_f = args[3]
sign_m = args[4]
cancer_type =  args[5]

#f = 'ENSG00000131018-rs10499233_exon_skip_462972-BRCA-sQTL.txt'
#svgf = 'ENSG00000131018-rs10499233_exon_skip_462972-BRCA-sQTL.svg'
#sign_f='***'
#sign_m='ns'
size=10
anno <- data.frame(xstar = c(2, 2), ystar = c(1.05, 1.05),lab = c(sign_f, sign_m),Group = c("Female", "Male"))
d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
p1 <- ggplot(d, aes(x=x, y=y,color=Group)) +geom_boxplot() + facet_grid(.~Group)
p1 <- p1 +labs(x='',title = cancer_type,y="PSI")
#p1 <- p1 + ylim(0,1.1)
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
#p1
svg(svgf,width=4,height=2.5)
print(p1)
dev.off()









