#RBG:219 82 137 #DB5289 pink
#RBG:105 128 189 #6980BD blue

#setwd('Downloads/down_from36/')
library(igraph)
library(RColorBrewer)
#f1 = 'ENSG00000000003_TF-Gene_male.txt'
#f2 = 'ENSG00000000003_TF-Gene_female.txt'
#svgf='ENSG00000000003_TF-Gene_network.svg'

args=commandArgs(T)
f1 = args[1]
f2 = args[2]
svgf = args[3]

net1<-read.table(f1,sep='\t',header =TRUE, stringsAsFactors = FALSE)
nodes1 <- data.frame(
  name<-c(unique(net1$from),unique(net1$to)),
  Gtype=c( rep("TF",length(unique(net1$from))),rep("Gene",length(unique(net1$to)))))

col <- data.frame(Gtype = unique(nodes1$Gtype), stringsAsFactors = F)
col$color <- c('grey','#DB5289')
nodes1$color <- col$color[match(nodes1$Gtype, col$Gtype)]
tmp1 <- graph_from_data_frame(net1,vertices = nodes1,directed=TRUE)

net2<-read.table(f2,sep='\t',header =TRUE, stringsAsFactors = FALSE)
nodes2 <- data.frame(
  name<-c(unique(net2$from),unique(net2$to)),
  Gtype=c( rep("TF",length(unique(net2$from))),rep("Gene",length(unique(net2$to)))))
col <- data.frame(Gtype = unique(nodes2$Gtype), stringsAsFactors = F)
col$color <- c('grey','#6980BD')
nodes2$color <- col$color[match(nodes2$Gtype, col$Gtype)]
tmp2 <- graph_from_data_frame(net2,vertices = nodes2,directed=TRUE)
svg(svgf,12,6)
par(mfrow=c(1,2))
plot(tmp1, main='TF-Gene network in male tumor patients',
     layout= layout_nicely,#  layout.kamada.kawai,
     vertex.label.cex = .75,
     edge.arrow.size=.5,
     edge.label.cex = 0.75,
     edge.label =  E(tmp1)$cancer_type,
     edge.width = E(tmp1)$score,
     )

plot(tmp2,  main='TF-Gene network in femla tumor patients',
     layout= layout_nicely,#  layout.kamada.kawai,
     vertex.label.cex = .75,
     edge.arrow.size=.5,
     edge.label.cex = 0.75,
     edge.label =  E(tmp2)$cancer_type,
     edge.width = E(tmp2)$score,
     )

dev.off()



