library(ggplot2)

args=commandArgs(T)
f = args[1]
svgf = args[2]
#f = '/data8/myang9/ExonSkipAD/output/diff_analysis/boxplot_diffSE/file/CB/SE_ZNF75A_chr16_+_3313048_3313175_3311619_3311948_3316911_3316924.txt'
#svgf = '/data8/myang9/ExonSkipAD/output/diff_analysis/boxplot_diffSE/figure/CB/SE_ZNF75A_chr16_+_3313048_3313175_3311619_3311948_3316911_3316924.svg'
SE_id = args[3]

boxplot_for_SE <- function(f,svgf){
  d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
  #print(head(d))
  d$tissue_p = sub(",","\n", d$tissue_p)
  print(head(d))
  p1 <- ggplot(d,mapping = aes(x =tissue_p, y =psi,fill= diagnosis)) + geom_violin(color="white") + geom_boxplot(width=0.2,position=position_dodge(0.9))
  p1 <- p1 +labs(x="",y="PSI",title = SE_id)
  p1 <- p1 + ylim(0,1) 
  p1 <- p1 +theme_bw() + theme(# axis.ticks.x=element_blank(),
                 axis.text.x = element_text(size = 12,color = 'black', vjust =0.5, hjust = 0, angle = 90),
                 axis.text.y = element_text(size = 12,color = 'black', vjust = 0, hjust =1, angle = 0),
                 axis.title.x =  element_text(size = 12,color = 'black'),
                 axis.title.y =  element_text(size = 12,color = 'black'),
                 title =  element_text(size = 12,color = 'black'),
                 plot.title = element_text(hjust = 0.5),
		 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
		 axis.line = element_line(size=0.5, colour = "black"),
		 panel.border = element_rect(colour ="black",size=1)
		)
  svg(svgf,width=8,height=4.5)
  print(p1)
  dev.off()
  }

boxplot_for_SE(f,svgf)
