library(data.table)
library(MASS)
options(stringsAsFactors=F)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)


####################################
########## Prepare prior ###########
####################################
prepare <- function(CT, prior, Npair, p_adj_method){
  #select CTs that are in prior data 
  p_CT = paste0('p_',CT)
  CT1 = CT[p_CT %in% colnames(prior)]
  
  for (i in 1:length(CT1)){
    
    prior[,paste0('p_adj_',CT1[i])] = p.adjust(prior[,paste0('p_',CT1[i])],method=p_adj_method, n=Npair)
    prior[prior[,paste0('p_adj_',CT1[i])]>0.9, paste0('p_adj_',CT1[i])] = 1
    prior[,paste0('beta_adj_',CT1[i])] = prior[,paste0('beta_',CT1[i])]*(1-prior[,paste0('p_adj_',CT1[i])])
  }
  
  return(prior)
}

########## read data (3 required input)
Prior.file <- as.character(args[1])
W.file <- as.character(args[2])
output.file <- as.character(args[3])

W = as.data.frame(fread(W.file, header=T))
prior = as.data.frame(fread(Prior.file, header=T))
CT = colnames(W)[-1]

########## 2 default arguments (optional)
Npair <- gsub(x = args[grep(x = args, pattern = "Npair=")], pattern = "Npair=", replacement = "")
p_adj_method <- gsub(x = args[grep(x = args, pattern = "p_adj_method=")], pattern = "p_adj_method=", replacement = "")
Npair = ifelse(length(Npair)==0, nrow(prior), as.numeric(Npair))
p_adj_method = ifelse(length(p_adj_method)==0, 'bonferroni', p_adj_method)
print(Npair)
print(p_adj_method)


########## run the algorithm
prior_adj = prepare(CT, prior, Npair, p_adj_method)
write.table(prior_adj, output.file, quote = F, row.names = F, col.names = T, sep ="\t")

