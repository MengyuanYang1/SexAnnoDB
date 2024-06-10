library(netZooR)
library(data.table)
library(reticulate)
library(qsmooth)
use_python("/usr/bin/python3")
library(yarn)        
library(RColorBrewer) 
library(tidyverse)


pandaDiffEdges <- function(panda.net1, panda.net2, threshold=0.98, condition_name="cond.1"){
  
  # reshape two PANDA networks
  merged.net <- merge(panda.net1, panda.net2, by=c("TF","Gene"))
  # edge-weight differences
  merged.net$'x-y' <- merged.net$Score.x-merged.net$Score.y
  merged.net$'y-x' <- merged.net$Score.y-merged.net$Score.x
  
  # CDF
  fnx <- ecdf(merged.net$'Score.x')
  fny <- ecdf(merged.net$'Score.y')
  fnxy <- ecdf(merged.net$'x-y')
  fnyx <- ecdf(merged.net$'y-x')
  
  # filter edge weights
  sub.net1 = merged.net
  sub.net1$'net1_threshold' = fnx(merged.net$Score.x)*(fnxy(merged.net$'x-y'))
  #sub.net1 <- merged.net[fnx(merged.net$Score.x)*(fnxy(merged.net$'x-y'))>threshold,]
  sub.net1[,c(condition_name)] <- "turmol_male"
  sub.net1 <- sub.net1[,c("TF","Gene","Score.x",'net1_threshold', condition_name)]
  
  colnames(sub.net1) <- c("TF","Gene","net1_Score",'net1_threshold',"net1_name")
  
  sub.net2 = merged.net
  sub.net2$'net2_threshold' = fny(merged.net$Score.y)*(fnyx(merged.net$'y-x'))
  #sub.net2 <- merged.net[fny(merged.net$Score.y)*(fnyx(merged.net$'y-x'))>threshold,]
  sub.net2[,c(condition_name)] <- "turmol_female"
  sub.net2 <- sub.net2[,c("Score.y", 'net2_threshold',condition_name)]
  colnames(sub.net2) <- c("net2_Score",'net2_threshold',"net2_name")
  
  # merge two panda networks
  merge.sub.net <- cbind(sub.net1,sub.net2)
  return(merge.sub.net)
}






i = "THCA"

load_file_name = paste("saved_",i,".RData",sep = "", collapse = "")
  load(load_file_name)

  diffRes_turmol_male_turmol_female <- pandaDiffEdges(turmol_male_net_reshaped, turmol_female_net_reshaped, condition_name="turmol_male")
  

diffRes_turmol_male_turmol_female_filter = diffRes_turmol_male_turmol_female %>%
	dplyr::filter((net1_threshold > 0.98) | (net2_threshold>0.98))



  diff_edge_name_diffRes_turmol_male_turmol_female = paste(i,"_threshold_diff_edge","_turmol_male_turmol_female",".txt",sep = "", collapse = "")
  diff_edge_name_diffRes_turmol_male_turmol_female_filter = paste(i,"_filter_threshold_diff_edge","_turmol_male_turmol_female",".txt",sep = "", collapse = "")



  write.table(diffRes_turmol_male_turmol_female ,file=diff_edge_name_diffRes_turmol_male_turmol_female,sep = "\t",row.names = F,quote = F)
    write.table(diffRes_turmol_male_turmol_female_filter ,file=diff_edge_name_diffRes_turmol_male_turmol_female_filter,sep = "\t",row.names = F,quote = F)

  
  print("bbbbb")
  
  
