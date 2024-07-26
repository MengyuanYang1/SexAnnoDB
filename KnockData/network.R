library(netZooR)
library(data.table)
library(reticulate)
library(ggplot2)
library(pROC)
library(ROCR)
library(forecast)

setwd('/Users/yangmengyuan/Documents/All_Project/SexAnnoDB/0-SexAnnoDB_20240704/RBP-ES-network/KnockData/')

motif = read.table('PANDA_RBPmotif_ALl.txt',header = T,check.names=FALSE)
ppi =  read.table("my_ppi_rbp_89",header = T,check.names=FALSE)
psi_matrix_knock <- read.table("PANDA_Knock_matrix.txt",header = T,check.names=FALSE,sep="\t",row.names=1)


panda_knock <- panda(motif,psi_matrix_knock,ppi,mode='intersection',progress=TRUE)#
panda_knock_net <- panda_knock@regNet
panda_knock_0 <- reshape2::melt(panda_knock_net)
panda_knock_reshaped <- data.frame("RBP"=as.character(panda_knock_0[,1]),"Gene"=as.character(panda_knock_0[,2]),"Motif"=NA,"Score"=as.numeric(panda_knock_0[,3]),stringsAsFactors = FALSE)

d_diff = read.table('RBP-ES_regulation.txt',header = T,check.names=FALSE)
d_diff$pairs <- paste(d_diff$RBP, d_diff$ES, sep="-")
colnames(panda_knock_reshaped) <-  c('RBP','ES','exp','pre')
d<- panda_knock_reshaped
d$pairs <- paste(d$RBP, d$ES, sep="-")
d$exp = 0
d$exp[d$pairs %in% d_diff$pairs ] <- 1
d$pre1[d$pre >quantile(d$pre,0.9)] <- 1 ; d$pre1[d$pre <quantile(d$pre,0.9)] <- 0
modelroc <- roc(d$exp,d$pre1)
p <- plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
p
outf <- 'ROC.pdf'
pdf(outf,4.5,4.5)
p <- plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
          grid.col=c("green", "red"), max.auc.polygon=TRUE,
          auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()



