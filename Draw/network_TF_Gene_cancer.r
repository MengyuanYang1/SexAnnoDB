#RBG:219 82 137 #DB5289 pink
#RBG:105 128 189 #6980BD blue
###E6E6FA purple

#setwd('Downloads/down_from36/')

library(igraph)
library(RColorBrewer)
args=commandArgs(T)
f1 = args[1]
svgf = args[2]
cancer_type = args[3]
#f1 = 'CHOL_RBP_ES_Summary.txt'
#svgf='CHOL_RBP_ES_Summary.svg'
#cancer_type = 'CHOL'

net <- read.table(f1,sep='\t',header =TRUE, stringsAsFactors = FALSE)
net1 <- net[which(net$Rtype==TRUE),]
net1 <- net1[order(-as.numeric(net1$score)),][1:20,]

nodes1 <- data.frame(
  name<-c(unique(net1$from),unique(net1$to)),
  Gtype=c( rep("TF",length(unique(net1$from))),rep("Gene",length(unique(net1$to)))))

col <- data.frame(Gtype = unique(nodes1$Gtype), stringsAsFactors = F)
col$color <- c('#E6E6FA','#DB5289')
nodes1$color <- col$color[match(nodes1$Gtype, col$Gtype)]
tmp1 <- graph_from_data_frame(net1,vertices = nodes1,directed=TRUE)

net2 <- net[which(net$Rtype==FALSE),]
net2 <- net2[order(-as.numeric(net2$score)),][1:20,]
nodes2 <- data.frame(name<-c(unique(net2$from),unique(net2$to)),
Gtype=c( rep("TF",length(unique(net2$from))),rep("Gene",length(unique(net2$to)))))
col <- data.frame(Gtype = unique(nodes2$Gtype), stringsAsFactors = F)
col$color <- c('#E6E6FA','#6980BD')
nodes2$color <- col$color[match(nodes2$Gtype, col$Gtype)]
tmp2 <- graph_from_data_frame(net2,vertices = nodes2,directed=TRUE)
svg(svgf,12,6)
par(mfrow=c(1,2))
plot(tmp1, main=paste0('Top20 male biased TF-Gene regulation in ',cancer_type),
     layout= layout_with_kk,   #  layout.kamada.kawai,
     vertex.label.cex = .75,
     edge.arrow.size=1,
     edge.label.cex = .75,
     #edge.label =  E(tmp1)$cancer_type,
     edge.width = E(tmp1)$score,
)

plot(tmp2,  main=paste0('Top20 female biased TF-Gene regulation in ',cancer_type),
     layout= layout_with_kk,#  layout.kamada.kawai,
     vertex.label.cex = .75,
     edge.arrow.size=1,
     edge.label.cex = 0.75,
     #edge.label =  E(tmp2)$cancer_type,
     edge.width = E(tmp2)$score,
)

dev.off()

#out <- rbind(net1,net2)
net1 <- net1[,1:3]
net2 <- net2[,1:3]
names(net1) <- c('TF','Gene','Score')
names(net2) <- c('TF','Gene','Score')
write.table(net1,paste0('/data2/myang9/SexAnnoDB/output/Network/TF_Gene/',cancer_type,'_TF_Gene_Summary_Maletop20.txt'),sep="\t",row.names=FALSE,quote=FALSE)  
write.table(net2,paste0('/data2/myang9/SexAnnoDB/output/Network/TF_Gene/',cancer_type,'_TF_Gene_Summary_Femaletop20.txt'),sep="\t",row.names=FALSE,quote=FALSE)



