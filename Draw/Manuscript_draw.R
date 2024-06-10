setwd('/Users/yangmengyuan/Documents/All_Project/SexAnnoDB/SexAnnoDB_20240510/Figures/')
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(ggsignif)
library(stringr)
library(reshape2)
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(forcats)
library(tidyverse)  #包含ggplot2包的重要R包集合
library(network)
library(gg.gap)
library(ggbreak)
library(patchwork)
require(scales)#数据缩放
library(ggtree)#聚类


#### '#DB5289' female
#### '#6980BD' male

#extrafont::loadfonts()
size = 8 
mytheme <- theme_bw() + theme(
  #text=element_text(family="Times New Roman"),
  axis.text.x = element_text(size = size,color = 'black', angle = 0),
  axis.text.y = element_text(size = size,color = 'black', angle = 0),
  axis.title.x =  element_text(size = size,color = 'black'),
  axis.title.y =  element_text(size = size,color = 'black'),
  title =  element_text(size = size,color = 'black'),
  plot.title = element_text(hjust = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  #axis.line = element_line(size=0.5, colour = "black"),
  #panel.border = element_rect(colour ="black",size=0.5)
  #legend.position = 'none'
)



###-----Datasummary----
f = 'Muti-omic_dataSummary.txt'
d <- read.table(f,header=TRUE,sep='\t')
#d$sex_muti <- paste(d$Cancer,d$DataType,d$Type)
d$color <- ifelse(d$Number <= 5, '<5', '>=5')
p <- ggplot(d, aes(x=Type, y=DataType, size=Number,color=color)) +  geom_point(alpha = 0.6) + facet_wrap(~Cancer) +  theme_bw()+labs(x=NULL,y=NULL)
p <- p + theme(axis.text.x = element_text(angle =90,))

p
pdf(outf,width=8,height=6)
print(p)
dev.off()

f = 'QTLSummary_dataSummary.txt'
d <- read.table(f,header=TRUE,sep='\t')
p <- ggplot(d,aes(x=Cancer,y=Number,fill=Group))+geom_bar(stat = 'identity', position = 'dodge', color='black')+labs(x=NULL,y='Number of patients')
p <- p + mytheme + geom_hline(yintercept = 50,linetype = "dashed", size = 0.3)  
p <- p + facet_grid(QTLtype~ .)
p 
outf <- "QTLSummary_dataSummary.pdf"
pdf(outf,width=8,height=3)
print(p)
dev.off()

###-----Resultsummary----
f = 'SignatureResultSummary.txt'
d1 <- d[which(d$Group == 'Female_specific'), ]
d2 <- d[which(d$Group == 'Male_specific'), ]
d3 <- d[which(d$Group == 'Tumor_specific_sex_biased'), ]
p <- ggplot(d1, aes(cancer,Nmuber)) + geom_bar(stat="identity",fill = '#DB5289',color = 'white', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='#Number of genes',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p1 <- p +   geom_text(aes(label=Nmuber), vjust=.5,hjust=0.5,colour="black",size=4)
p1
outf = "SignatureResultSummary_p1.pdf"
pdf(outf,width=5,height=3)
p1 
dev.off()
p <- ggplot(d2, aes(cancer,Nmuber)) + geom_bar(stat="identity",fill = '#6980BD',color = 'white', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='Number of genes',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p2 <- p +   geom_text(aes(label=Nmuber), vjust=.5,hjust=0.5, colour="black",size=4)
p2
outf = "SignatureResultSummary_p2.pdf"
pdf(outf,width=5,height=3)
p2 
dev.off()
p <- ggplot(d3, aes(cancer,Nmuber)) + geom_bar(stat="identity",fill = 'grey',color = 'white', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='Number of genes',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p <- p +   geom_text(aes(label=Nmuber), vjust=.5,hjust=0.5, colour="black",size=4)
p
outf = "SignatureResultSummary_sex-biased.pdf"
pdf(outf,width=7,height=3)
p
dev.off()

f = 'SignatureResultSummary_Tspecifc.txt'
d <- read.table(f,header=TRUE,sep='\t')
d$Cancer_Number <- as.character(d$Cancer_Number)
d <- d[which(d$Signature == 'CodingGene'), ]
d$type <- paste(d$Group,d$Signature,sep='_')
d1 <- d[which(d$Group == 'Female_specific'), ]
d2 <- d[which(d$Group == 'Male_specific'), ]
d3 <- d[which(d$Group == 'Tumor_specific_sex_biased'), ]
p <- ggplot(d1, aes(Cancer_Number,Signature_N)) + geom_bar(stat="identity",fill = '#DB5289',color = 'white', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='Number of cancers',y='Number of genes',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <-p +  scale_y_break(c(10, 25),1.5) + scale_y_break(c(180, 600),1.5)
p3 <- p +   geom_text(aes(label=Signature_N), vjust=.5, colour="black",size=4)
p3 
outf = "SignatureResultSummary_p3.pdf"
pdf(outf,width=3.5,height=3)
p3 
dev.off()
p <- ggplot(d2, aes(Cancer_Number,Signature_N)) + geom_bar(stat="identity",fill = '#6980BD',color = 'white', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='Number of cancers',y='Number of genes',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <-p +  scale_y_break(c(15,70),1.2)  +  scale_y_break(c(600, 900),1.5)
p4 <- p +   geom_text(aes(label=Signature_N), vjust=.5,hjust=0.5, colour="black",size=4)
p4 
outf = "SignatureResultSummary_p4.pdf"
pdf(outf,width=3.5,height=3)
p4 
dev.off()

d3$Cancer_Number<- factor(d3$Cancer_Number, levels = c('1','2','3','4','5','6','7','8','9','15','16'),ordered = TRUE)
p <- ggplot(d3, aes(Cancer_Number,Signature_N)) + geom_bar(stat="identity",fill = 'grey',color = 'white', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='Number of cancers',y='Number of genes',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <-p +  scale_y_break(c(20,40),1.2)  +  scale_y_break(c(400, 1200),0.5)
p <- p +   geom_text(aes(label=Signature_N), vjust=.5,hjust=0.5, colour="black",size=4)
p
outf = "SignatureResultSummary_sex-biased2.pdf"
pdf(outf,width=3.5,height=3)
p
dev.off()

###--------EditingSummary------
f = 'EditingSNVsummary.txt'
d <- read.table(f,header=TRUE,sep='\t')
d <- d[which(d$Group == 'Tumor_specific_sex_biased'), ]
p <- ggplot(d, aes(Cancer,Number)) + geom_bar(stat="identity",position="dodge",color = 'white',fill='grey', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
#p <- p +   scale_fill_manual(values = c('#DB5289','#6980BD','grey'))
#p <-p +  scale_y_break(c(200,300),0.5) 

p <- p +   geom_text(aes(label=Number), vjust=.5,hjust=0.5,colour="black",size=3)
p

outf = "Editing_summary.pdf"
pdf(outf,width=6,height=2.5)
p
dev.off()

f = 'EditingSNVsummary.txt'
d <- read.table(f,header=TRUE,sep='\t')
d$Group
d <- d[which(d$Group == 'Female_specific'), ]
p <- ggplot(d, aes(Cancer,Number)) + geom_bar(stat="identity",position="dodge",color = 'white',fill='#DB5289', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p <- p +   geom_text(aes(label=Number), vjust=.5,hjust=0.5,colour="black",size=3)
p
outf = "Editing_summary_F.pdf"
pdf(outf,width=3,height=2)
p
dev.off()
f = 'EditingSNVsummary.txt'
d <- read.table(f,header=TRUE,sep='\t')
d$Group
d <- d[which(d$Group == 'Male_specific'), ]
p <- ggplot(d, aes(Cancer,Number)) + geom_bar(stat="identity",position="dodge",color = 'white',fill='#6980BD', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p <- p +   geom_text(aes(label=Number), vjust=.5,hjust=0.5,colour="black",size=3)
p
outf = "Editing_summary_M.pdf"
pdf(outf,width=3,height=2)
p
dev.off()

f = 'Editing_candidate.txt'
d <- read.table(f,header=TRUE,sep='\t')
p <- ggplot(data=d, aes(x=cancer, y=chr11_61007595_., color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='Frequency')
p <- p + mytheme  +  theme(legend.position = 'none')#+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
outf = 'Editing_candidate.pdf'
pdf(outf,2.2,1.8)
p
dev.off()

f = 'Editing_candidate_SARC_CPM.txt'
d <- read.table(f,header=TRUE,sep='\t')
p <- ggplot(data=d, aes(x=cancer, y=chr12_68851361_., color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='Frequency')
p <- p + mytheme  +  theme(legend.position = 'none')#+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
outf = 'Editing_candidate_SARC_CPM.pdf'
pdf(outf,1.5,1.8)
p
dev.off()


f = 'ENSG00000135678_SARC_expression.txt'
d <- read.table(f,header=TRUE,sep='\t')
d <- d[which(d$type %in% c('F_T','M_T')), ]
p <- ggplot(data=d, aes(x=type, y=log2(exp), color=type)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='log2(FPKM)')#+ylim(0,6)
p <- p + mytheme  +  theme(legend.position = 'none')#+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
outf = 'ENSG00000135678_SARC_expression.pdf'
pdf(outf,1.5,1.8)
p
dev.off()

f = 'ENSG00000163220.txt'
d <- read.table(f,header=TRUE,sep='\t')
d <- d[which(d$group %in% c('F_T','M_T')), ]
d <- d[which(d$cancer_type %in% c('BRCA','DLBC','KIRP','SARC')), ]


p <- ggplot(data=d, aes(x=cancer_type, y=log2(FPKM), color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='S100A9 log2(FPKM)')#+ylim(0,6)
p <- p + mytheme  +  theme(legend.position = 'none')#+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
outf = 'ENSG00000163220.pdf'
pdf(outf,2.5,1.8)
p
dev.off()


d <- d[which(d$Group == 'Male_specific'), ]
p <- ggplot(d, aes(Cancer,Number)) + geom_bar(stat="identity",position="dodge",color = 'white',fill='#6980BD', width = 0.8) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='',title = '')#+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p <- p +   geom_text(aes(label=Number), vjust=.5,hjust=0.5,colour="black",size=3)
p
outf = "Editing_summary_M.pdf"
pdf(outf,width=3,height=2)
p
dev.off()



###----------MutationSummary-------------
f = 'SexAnnoDB_mutation_enrichment.txt'
d <- read.table(f,header=TRUE,sep='\t')
data_m <- melt(d,id.vars=c("Term")) 
p <- ggplot(data_m,aes(x=variable,y=Term,fill=-log10(value))) #热图绘制
p <- p +geom_raster()+scale_fill_gradient2(low="#003366", high="red", na.value="white")
p <- p + theme(axis.text.x = element_text(angle =90))
p

outf <- 'SexAnnoDB_mutation_enrichment.pdf'
pdf(outf,6,4)
p
dev.off()


f = 'ENSG00000196557.txt'
outf = 'ENSG00000196557_CACNA1H.pdf'
d <- read.table(f,header=TRUE,sep='\t')
d <- d[which(d$cancer_type == 'SARC'), ]
d <- d[which(d$group %in% c('F_T','M_T')), ]

p <- ggplot(data=d, aes(x=group, y=FPKM, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
pdf(outf,2.2,1.8)
p
dev.off()


f = 'SexAnnoDB_Sexbiased_mutation2.txt'
d <- read.table(f,header=TRUE,sep='\t')
d_can <-DataFrame(table(d$CancerType))
#p <- ggplot(d_can, aes(x= reorder(Var1,-Freq),Freq)) + geom_bar(stat="identity",fill = 'white',color = 'black', width = 0.7) #+ coord_flip()
p <- ggplot(d_can, aes(x= Var1,Freq)) + geom_bar(stat="identity",fill = 'white',color = 'black', width = 0.7) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='Number of genes',title = '') #+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p
outf <- 'SexAnnoDB_Sexbiased_mutation2_genenumberforeachgene.pdf'
pdf(outf,4.5,2.5)
p
dev.off()

d_gene <-DataFrame(table(d$Hugo_Symbol))
d_gene <- d_gene[which(d_gene$Freq >2 ), ]

p <- ggplot(d_gene, aes(x= reorder(Var1,-Freq),Freq)) + geom_bar(stat="identity",fill = 'white',color = 'black', width = 0.7) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='Number of genes',title = '') #+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))
p
outf <- 'SexAnnoDB_Sexbiased_mutation2_cancernumberforeachgene.pdf'
pdf(outf,4.5,2.5)
p
dev.off()



###----------MethylationSummary-------------
f = 'SexbiasedDpG_in_SexbiasedGene.txt'
d <- read.table(f,header=TRUE,sep='\t')
d_can <-DataFrame(table(d$Cancer))
#p <- ggplot(d_can, aes(x= reorder(Var1,-Freq),Freq)) + geom_bar(stat="identity",fill = 'white',color = 'black', width = 0.7) #+ coord_flip()
p <- ggplot(d_can, aes(x= Var1,Freq)) + geom_bar(stat="identity",fill = 'white',color = 'black', width = 0.7) #+ coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='Number of genes',title = '') #+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + theme(axis.text.x = element_text(angle =90,))+ ylim(0,1300)
p <-p +  scale_y_break(c(300, 1000),0.5) 
p1 <- p +   geom_text(aes(label=Freq), vjust=-1,hjust=0.5,colour="black",size=3)
p1
p
outf <- 'SexbiasedDpG_in_SexbiasedGene.pdf'
pdf(outf,5,3)
p1
dev.off()

f = 'Figure2_GPC3_CpGsites.txt'
d <- read.table(f,header=TRUE,sep='\t')
d_tmp <- d[which(d$CpGids != 'group'), ]
d_tmp <- d[which(d_tmp$CpGids != 'cg00179598'), ]

data <- melt(d_tmp, id="CpGids")
data<-na.omit(data)
data$value <- as.numeric(data$value)

d_group$type  <- as.vector(d[which(d$CpGids == 'group'), ])
d_group$label <- c(rep("group",dim(d_group)[1]))
d_group$label
data <- ggplot(d,aes(.,y=p,fill=group))

p <-ggplot(data,aes(x=variable,y=CpGids,fill=value)) +geom_raster()+scale_fill_gradient2(low="#003366", high="#990033", mid="white")
p <- p + geom_vline(xintercept=c(144),size=.5)
p <- p + scale_y_discrete(position="left")
p <- p + theme(axis.text.x =element_blank(),axis.text.y =element_blank(),legend.position = "none")
#p <- p + insert_top(group, height = .05)
p

outf = 'Figure2_GPC3_CpGsites.pdf'
pdf(outf,3,2.5)
p
dev.off()


f = 'Figure2_GPC3_CpGsites.txt'
d <- read.table(f,header=TRUE,sep='\t', row.names = 1)
d <- d[which(rownames(d) != 'cg00179598'), ]

d <- as.data.frame(t(d))
d <- melt(d, id.vars="group")
d$value <- as.numeric(d$value)
p <- ggplot(data=d, aes(x=variable, y=value, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + mytheme
p <- p + theme(axis.text.x = element_text(angle =90,))
#### '#6980BD' male
p <- p + labs(x='',y='Beta vale')
p <- p + ylim(0,1)
# +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
outf <- 'Figure2_GPC3_CpGsites_supplementary.pdf'
pdf(outf,6,2.5)
p
dev.off()

###--------Sex-biased gene------------------
boxplot_for_isoform <- function(f,svgf){
  d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
  d <- data.frame(d)
  y_position = max(d$FPKM)*0.9
  d_tmp = d[,c('cancer_type','sign')]
  d_tmp = d_tmp[order(d_tmp$cancer_type),]
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
  svg(svgf,width=3.4,height=2)
  print(p1)
  dev.off()
}
drawfigure <- function(d,d_sign,enst,x,y,group){
  p1 <- ggboxplot(d,x,y,color =group) 
  p1 <- p1 +labs(x=enst,y="log2(FPKM)")
  p1 <- p1 + scale_colour_manual(values = c('#DB5289','#6980BD'))
  #p1 <- p1 + scale_colour_manual(values = c('#DB5289','black'))
  p1 <- p1 + geom_signif(xmin=d_sign$sign_x,xmax=d_sign$sign_xend, y_position=d_sign$y_position, annotation=d_sign$sign)
  size = 8
  p1 <- p1 +theme_bw() + theme(# axis.ticks.x=element_blank(),
    axis.text.x = element_text(size = size,color = 'black', vjust =0.2, hjust = 0.5, angle =0),
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
f = 'ENSG00000163220_TS.txt'
svgf = 'Figure1_Sexbiased_gene.svg'
boxplot_for_isoform(f,svgf)

####sup Figure1
boxplot_for_isoform <- function(f,svgf){
  d <- read.table(f ,header =TRUE, stringsAsFactors = FALSE);
  d <- data.frame(d)
  y_position = max(d$FPKM)*0.9
  d_tmp = d[,c('cancer_type','sign')]
  d_tmp = d_tmp[order(d_tmp$cancer_type),]
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
  svg(svgf,width=3.4,height=2)
  print(p1)
  dev.off()
}
drawfigure <- function(d,d_sign,enst,x,y,group){
  groups = c('F_T','F_N') 
  d$group <- factor(d$group,levels = groups)
  p1 <- ggboxplot(d,x,y,color =group) 
  p1 <- p1 +labs(x=enst,y="log2(FPKM)")
  #p1 <- p1 + scale_colour_manual(values = c('#DB5289','#6980BD'))
  p1 <- p1 + scale_colour_manual(values = c('#DB5289','black'))
  p1 <- p1 + geom_signif(xmin=d_sign$sign_x,xmax=d_sign$sign_xend, y_position=d_sign$y_position, annotation=d_sign$sign)
  size = 8
  p1 <- p1 +theme_bw() + theme(# axis.ticks.x=element_blank(),
    axis.text.x = element_text(size = size,color = 'black', vjust =0.2, hjust = 0.5, angle =0),
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
f = 'ENSG00000163220_FS.txt'
svgf = 'SupFigure1_Female_spcific_gene.svg'
boxplot_for_isoform(f,svgf)


f = 'Top20_KEGGPathway_SexGene.txt'
d <- read.table(f,header=TRUE,sep='\t')
d$P <- as.numeric(d$P)
d$y <- -log10(d$P)
d$y <- as.numeric(d$y)
col <-  reorder(d$group,d$y)
#d$Term <- factor(d$Term,levels = terms)
p <- ggplot(d, aes(x= reorder(Term,y),y)) + geom_bar(stat="identity",fill = 'red',alpha=.6,color = 'white', width = 0.5) + coord_flip()
p <- p + mytheme + labs(x='',y='-log10(pvalue)',title = '') + theme(axis.text.x = element_text(colour =col))
p

svg('Top20_KEGGPathway_SexGene.svg',5,4)
p
dev.off()


f = 'ENSG00000147257.txt'
svgf = 'ENSG00000147257_GPC3.svg'
d <- read.table(f,header=TRUE,sep='\t')
d <- d[which(d$cancer_type == 'ACC'), ]
p <- ggplot(data=d, aes(x=group, y=log2(FPKM), color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
#### '#6980BD' male
d$sign
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svg(svgf,2.2,1.8)
p
dev.off()

f = 'GPC3_CpGsites.txt'
svgf = 'GPC3_CpGsites.svg'
d <- read.table(f,header=TRUE,sep='\t')
#d <- d[which(d$cancer_type == 'ACC'), ]
p <- ggplot(data=d, aes(x=CpGsite, y=betavalue, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
#### '#6980BD' male
p <- p + labs(x='',y='Beta vale')
p <- p + ylim(0,1)
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svg(svgf,3.5,1.8)
p
dev.off()


###---TF-Gene network-----------
f = 'S100A9_SARC_TF_cor.txt'
d <- read.table(f,header=TRUE,sep='\t')
d_f <- d[which(d$gander == 'female'), ]
p <- ggplot(data=d_f, aes(x=ENSG00000163220, y=ENSG00000197124,color=gander))+geom_point(alpha=0.8, size=0.8)+geom_smooth(method = lm,linetype='dashed',colour='black',size=0.2)
p <- p + labs(x='S100A9',y='TF')
p <- p  +scale_colour_manual(values = c('#DB5289'))
p <- p + mytheme #+  theme(legend.position = 'none')
#p <- p + xlim(0,15)
p
d_f_cor <- d_f[,c('ENSG00000163220','ENSG00000126003','ENSG00000183647','ENSG00000198538','ENSG00000197124')]
d_cor <- cor(d_f_cor,method="pearson")
d_m <- d[which(d$gander == 'male'), ]
p <- ggplot(data=d_m, aes(x=ENSG00000163220, y=ENSG00000197124,color=gander))+geom_point(alpha=0.8, size=0.8)+geom_smooth(method = lm,linetype='dashed',colour='black',size=0.2)
p <- p + labs(x='S100A9',y='TF')
p <- p  +scale_colour_manual(values = c('#6980BD'))
p <- p + mytheme #+  theme(legend.position = 'none')
#p <- p + xlim(0,15)
p



f = 'SARC_volcanoplot.txt'
d <- read.table(f,header=TRUE,sep='\t')

p <- ggplot(data=d, aes(x=log2FoldChange, y=-log10(padj), color=Change)) +geom_point(alpha=0.8, size=0.75)
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD','black'))
#p <- p + labs(x='dPSI',y='-log10(pvalue)')
p <- p + geom_hline(aes(yintercept=-log10(0.05)), colour="black", linetype="dashed")
p <- p + geom_vline(aes(xintercept=-1), colour="black", linetype="dashed")
p <- p + geom_vline(aes(xintercept=1), colour="black", linetype="dashed")
p <- p + mytheme+  theme(legend.position = 'none')
p <- p +  xlab('log2FoldChange')+ ylab('-log10(pvalue)')
p <-  p +xlim(-10,10) +ylim(0,50)
p
svgf = 'SARC_volcanoplot.svg'
svg(svgf,2,2.5)
p
dev.off()


f = 'SARC_TF_Gene_Summary_sexBiasedTarget.txt'
svgf = 'SARC_TF_Gene_Summary_sexBiasedTarget.svg'
d <- read.table(f,header=TRUE,sep='\t')
tf <- unique(d$from)
to <- unique(d$to)
node <- c(tf,to)
net <- as.network(d , matrix.type='edgelist')
net %v% "shape" = ifelse(node %in% tf, "TF", "Gene")
p <- ggnet2(net, arrow.size = 12,label = TRUE,edge.size = 1,edge.color = d$Rtype,arrow.gap = 0.025,shape = "shape", shape.palette = c("TF" = 19, "Gene" = 15))

svg(svgf,12,4.5)
p
dev.off()

f = 'ENSG00000012817.txt'
d <- read.table(f,header=TRUE,sep='\t')
d1 <- d[which(d$cancer_type == 'SARC'), ]
d1$Gene <- c('KDM5D')
f = 'ENSG00000172156.txt'
d <- read.table(f,header=TRUE,sep='\t')
d2 <- d[which(d$cancer_type == 'SARC'), ]
d2$Gene <- c('CCL11')
f = 'ENSG00000163220.txt'
d <- read.table(f,header=TRUE,sep='\t')
d3 <- d[which(d$cancer_type == 'SARC'), ]
d3$Gene <- c('S100A9')
f = 'ENSG00000196616.txt'
d <- read.table(f,header=TRUE,sep='\t')
d4 <- d[which(d$cancer_type == 'SARC'), ]
d4$Gene <- c('ADH1B')
f = 'ENSG00000248144.txt'
d <- read.table(f,header=TRUE,sep='\t')
d4 <- d[which(d$cancer_type == 'SARC'), ]
d4$Gene <- c('ADH1C')
f = 'ENSG00000181072.txt'
d <- read.table(f,header=TRUE,sep='\t')
d5 <- d[which(d$cancer_type == 'SARC'), ]
d5$Gene <- c('CHRM2')
f = 'ENSG00000188322.txt'
d <- read.table(f,header=TRUE,sep='\t')
d6 <- d[which(d$cancer_type == 'SARC'), ]
d6$Gene <- c('SBK1')
d_out <- rbind(d1,d2,d3,d4,d5,d6)
d_out <- d_out[which(d_out$group %in% c('F_T','M_T')), ]

p <- ggplot(data=d_out, aes(x=Gene, y=FPKM, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svgf= 'SARC_TF_Gene_sex_biased_gene_drugmatch.svg'
svg(svgf,4,2)
p
dev.off()




###------RBP-ES network--------
f = 'LIHC_sexbiasedRBP_ES_net.txt'
svgf = 'LIHC_sexbiasedRBP_ES_net.svg'


f = 'SexBiased_edge_count.txt'
d <- read.table(f,header=TRUE,sep='\t')
d_f <- d[which(d$Type == 'FealeBiased_edges'), ]
p <- ggplot(d_f, aes(x= Cancer_type,Count)) + geom_bar(stat="identity",fill = '#DB5289',color = 'white', width = 0.7) + coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='#Female biased edges',title = '')+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p
svg('SexBiased_edge_count_2.svg',2.5,4.5)
p
dev.off()
d_m <- d[which(d$Type =='MaleBiased_edges'), ]
p <- ggplot(d_m, aes(x= Cancer_type,Count)) + geom_bar(stat="identity",fill = '#6980BD',color = 'white', width = 0.7) + coord_flip() 
p <- p + mytheme+theme_classic() + labs(x='',y='#Male biased edges',title = '') # +scale_x_discrete(position = "top")
p <- p + scale_y_reverse(position = "right")+scale_x_discrete(position = "top")
p
svg('SexBiased_edge_count_1.svg',2.5,4.5)
p
dev.off()

f = 'SexBiased_gene_drugid_count.txt'
d <- read.table(f,header=TRUE,sep='\t')
d <- d[which(d$type !='DrugTarget'), ]
d$type <- factor(d$type,levels = c('Sexbiased','allTarget'))
p <- ggplot(d, aes(x= cancer_type,count,fill=type)) + geom_bar(stat="identity",color = 'black',alpha=.5, width = 0.7) + coord_flip()
p <- p + mytheme+theme_classic() + labs(x='',y='',title = '')+ scale_y_continuous(position = "right") #+scale_y_discrete(position = "right")
p <- p + scale_fill_manual(values = c('red','white'))
p
svg('SexBiased_gene_drugid_count.svg',4.5,4.5)
p
dev.off()

f = 'RBP_ES_targetES_count.txt'
svgf = 'RBP_ES_targetES_count.svg'
d <- read.table(f,header=TRUE,sep='\t')
d <- d[which(d$Type =='TargetedES'), ]
#d$type <- factor(d$type,levels = c('Sexbiased','allTarget'))
p <- ggplot(d, aes(x= CancerType,y = Count,fill=Type)) + geom_bar(stat="identity",color = 'black',fill='grey', width = 0.7)
p <- p + mytheme + labs(x='',y='',title = '') + theme_classic()
p <- p + theme(axis.text.x = element_text(angle =90))
p
svg(svgf,6,2.5)
p
dev.off()

f = 'ENSG00000112183_RBM24_exp.txt'
svgf = 'ENSG00000112183_RBM24_exp.svg'
d <- read.table(f,header=TRUE,sep='\t')
p <- ggplot(data=d, aes(x=group, y=log2(FPKM), color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svg(svgf,2.2,1.8)
p
dev.off()

f = 'LIHC_sexbiasedRBP_ES_net.txt'
svgf = 'LIHC_sexbiasedRBP_ES_net.svg'
d <- read.table(f,header=TRUE,sep='\t')
from <- unique(d$from)
to <- unique(d$to)
node <- c(from,to)
net <- as.network(d , matrix.type='edgelist')
net %v% "shape" = ifelse(node %in% from, "RBP", "ES")
p <- ggnet2(net, arrow.size = 12,label = TRUE,label.size=4,edge.size = 1,edge.color = '#DB5289',arrow.gap = 0.025,shape = "shape", shape.palette = c("RBP" = 19, "ES" = 15))
p 
svg(svgf,10,8)
p
dev.off()

###---------------------ceRNA-------------------
f = 'ceRNA_exp.txt'
d = read.table(f,header = TRUE,sep='\t')
###---------Coding--------
p <- ggplot(data=d, aes(x=group, y=ENSG00000091831, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
#### '#6980BD' male
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svgf = 'FENSG00000091831_ESR1_exp.svg'
svg(svgf,2.2,1.8)
p
dev.off()
###---------lncRNA--------
p <- ggplot(data=d, aes(x=group, y=ENSG00000237125, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
#### '#6980BD' male
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svgf = 'ENSG00000091831_HAND2-AS1_exp.svg'
svg(svgf,2.2,1.8)
p
dev.off()
###---------miRNA--------
p <- ggplot(data=d, aes(x=group, y=hsa.mir.206, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
#### '#6980BD' male
p <- p + labs(x='',y='log2(RPM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svgf = 'hsa.mir.206_exp.svg'
svg(svgf,2.2,1.8)
p
dev.off()

#### '#DB5289' female
#### '#6980BD' male

###----------lncRNA-Coding--------
p <- ggplot(data=d, aes(x=ENSG00000091831, y=ENSG00000237125,color=group))+geom_point(alpha=0.8, size=0.8)+geom_smooth(method = lm,linetype='dashed',colour='black',size=0.2)
p <- p + labs(x='ESR1',y='HAND2-AS1')
p <- p  +scale_colour_manual(values = c('#6980BD','#DB5289'))
p <- p + mytheme #+  theme(legend.position = 'none')
#p <- p + xlim(0,15)
p
svgf = 'ceRNA_lncRNA_mRNA_cor.svg'
svg(svgf,3,1.8)
p
dev.off()


####--------------------QTL analysis-----------------
f = 'QTL_position_pie.txt'
d = read.table(f,header = TRUE,sep='\t',row.names=NULL)
l_type1 = c('eQTL','eQTM')
l_type2 = c('sQTL','sQTM')

d1 <- d[which(d$Ana_type %in% l_type1), ]
p <- ggplot(d1, aes(POsition,Count,fill =Ana_type )) +  geom_bar(stat="identity", position=position_dodge())# + coord_flip()
p <- p + mytheme + labs(x='',y='Count',title = '') 
#install.packages ("plotrix")
p
#p1 <- gg.gap(plot=p,segments=c(100,1500000),tick_width=10,rel_heights=c(0.25,0,0.1),ylim=c(0,1600000))
#p1


d2 <- d[which(d$Ana_type %in% l_type2), ]
Pos <- c('Distant upstream','Proximal upstream','Skipped exon','Proximal downstream','Distant downstream')
d2$POsition <- factor(d2$POsition,levels = Pos)

p <- ggplot(d2, aes(POsition,Count,fill =Ana_type )) +  geom_bar(stat="identity",alpha=.6,position=position_dodge())# + coord_flip()
p <- p + mytheme + labs(x='',y='Count',title = '') 
#install.packages ("plotrix")
p
p1 <- gg.gap(plot=p,segments=c(5000,36000),tick_width=c(100,1000),rel_heights=c(0.25,0,0.1),ylim=c(0,40000))
p1


for (x in l_type){
  print(x)
  d_tmp <-d[which(d$Ana_type == x),]
  print(d)
  d_percent <- prop.table(d_tmp$Count)
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  p<- ggplot(d_tmp, aes("", Count, fill = POsition)) +
    geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
    coord_polar("y")  +
    labs(x = NULL, y = NULL, fill = NULL, 
         title = "") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, color = "#666666"))
  print(p )
  
  
}


f = 'Durg_canDrug_count.txt'
svgf =  'Durg_canDrug_count.svg'
d <- read.table(f,header=TRUE,sep='\t',row.names=NULL)

p <- ggplot(d, aes(Ana,Count,fill =Type )) +  geom_bar(stat="identity",alpha=.8, position=position_dodge())# + coord_flip()
p <- p + mytheme + labs(x='',y='Count',title = '') 
#install.packages ("plotrix")
p

svg(svgf,4,2)
p
dev.off()



f = 'ENSG00000146648.txt'
d <- read.table(f,header=TRUE,sep='\t')
d1 <- d[which(d$cancer_type == 'BRCA'), ]
d1$Gene <- c('EGFR')
f = 'ENSG00000105810.txt'
d <- read.table(f,header=TRUE,sep='\t')
d2 <- d[which(d$cancer_type == 'BRCA'), ]
d2$Gene <- c('CDK6')
f = 'ENSG00000166501.txt'
d <- read.table(f,header=TRUE,sep='\t')
d3 <- d[which(d$cancer_type == 'BRCA'), ]
d3$Gene <- c('PRKCB')
f = 'ENSG00000131747.txt'
d <- read.table(f,header=TRUE,sep='\t')
d4 <- d[which(d$cancer_type == 'BRCA'), ]
d4$Gene <- c('TOP2A')
f = 'ENSG00000188389.txt'
d <- read.table(f,header=TRUE,sep='\t')
d5 <- d[which(d$cancer_type == 'BRCA'), ]
d5$Gene <- c('PDCD1')
f = 'ENSG00000055118.txt'
d <- read.table(f,header=TRUE,sep='\t')
d6 <- d[which(d$cancer_type == 'BRCA'), ]
d6$Gene <- c('KCNH2')
d_out <- rbind(d1,d2,d3,d4,d5,d6)
d_out <- d_out[which(d_out$group %in% c('F_T','M_T')), ]

p <- ggplot(data=d_out, aes(x=Gene, y=FPKM, color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))  ###### '#DB5289' female
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p 
svgf= 'BRCA_sex_biased_gene_drugmatch.svg'
svg(svgf,3.8,1.9)
p
dev.off()

###------Drug Target-------
f = 'BRCA_volcanoplot.txt'
d <- read.table(f,header=TRUE,sep='\t')

p <- ggplot(data=d, aes(x=log2FoldChange, y=-log10(padj), color=Change)) +geom_point(alpha=0.8, size=0.75)
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD','black'))
#p <- p + labs(x='dPSI',y='-log10(pvalue)')
p <- p + geom_hline(aes(yintercept=-log10(0.05)), colour="black", linetype="dashed")
p <- p + geom_vline(aes(xintercept=-1), colour="black", linetype="dashed")
p <- p + geom_vline(aes(xintercept=1), colour="black", linetype="dashed")
p <- p + mytheme+  theme(legend.position = 'none')
p <- p +  xlab('log2FoldChange')+ ylab('-log10(pvalue)')
p <-  p +xlim(-10,10) +ylim(0,50)
p
svgf = 'BRCA_volcanoplot.svg'
svg(svgf,2,2.5)
p
dev.off()

####
f = 'ENSG00000146648.txt'
d <- read.table(f,header=TRUE,sep='\t')
d1 <- d[which(d$cancer_type == 'BRCA'), ]
d1$Gene <- c('EGFR')

d_out <- d1
d_out <- d_out[which(d_out$group %in% c('F_T','M_T')), ]

p <- ggplot(data=d_out, aes(x=group, y=FPKM,color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))      ###### '#DB5289' female
p <- p + labs(x='',y='log2(FPKM)')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p <-  p +ylim(0,15)

p 
svgf= 'BRCA_EGFR.svg'
svg(svgf,2.1,1.9)
p
dev.off()
####----protein----
f = 'GFR_protein.txt'
d <- read.table(f,header=TRUE,sep='\t')
d_out <-d 
p <- ggplot(data=d_out, aes(x=group, y=EGFR,color=group)) +geom_boxplot()
p <- p  +scale_colour_manual(values = c('#DB5289','#6980BD'))      ###### '#DB5289' female
p <- p + labs(x='',y='')
p <- p + mytheme # +  theme(legend.position = 'none')+theme(axis.text.x = element_text(size = size,color = 'black', angle =0))
p <-  p +ylim(0,4)

p 
svgf= 'EGFR_protein.svg'
svg(svgf,2.2,1.9)
p
dev.off()









