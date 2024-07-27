library(data.table)
library(MASS)
options(stringsAsFactors=F)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)


####################################
########## HBI algorithm ###########
####################################
EM_IWLS <- function(x, y, s0, num_CT, a, b_vector,mu_vector){
  
  #initial value
  n = dim(x)[1]
  J = dim(x)[2]
  if (abs(rcond(t(x)%*%x)) < 1e-15){
    return(a = 'singular')
  }
  beta_new = solve(t(x)%*%x)%*%t(x)%*%y
  beta <- rep(100,J)
  phi_new <- 1
  beta_all = c(as.vector(beta_new))
  phi_all = c(phi_new)
  while(max(abs(beta_new-beta))>0.0001){
    beta <- beta_new
    phi <- phi_new
    beta_qtl_eff = beta[(J-num_CT+1):J]
    s <- c(rep(s0,J-num_CT), (1+a)/(abs(beta_qtl_eff-mu_vector)+b_vector) )
    W <- diag(c(rep(1,n),phi*s/abs(beta)))
    Z <- c(y, rep(0,J))
    X <- as.matrix(rbind(x, diag(J)))
    if (abs(rcond(t(X)%*%W%*%X)) < 1e-15){
      break
    }
    beta_new <- solve(t(X)%*%W%*%X) %*% t(X)%*%W%*%Z
    phi_new <- as.numeric(1/n * t(Z - X%*%beta_new) %*% W %*% (Z - X%*%beta_new))
    
    beta_all = rbind(beta_all, as.vector(beta_new))
    phi_all = c(phi_all, phi_new)
  }
  
  if (abs(rcond(t(X)%*%W%*%X)) < 1e-15){
    return(a = 'singular')
    
    
  }else{
    var = phi_new*solve(t(X)%*%W%*%X)
    return(list(beta_all=beta_all,phi_all=phi_all,
                beta = beta_new, var = var,
                phi = phi_new,
                s = s))
  }
  
}


####################################
##########  Running HBI  ###########
####################################

run_HBI = function(snp){
  G = as.data.frame(df_geno[, snp])
  colnames(G) = 'genotype'
  rownames(G) = rownames(df_geno)
  G$IID = rownames(G)
  G[G[,1]=='0/0', 1] = 0 
  G[G[,1]=='0/1', 1] = 1 
  G[G[,1]=='1/1', 1] = 2 
  G[G[,1]=='./.', 1] = NA
  G=na.omit(G)
  G[, 1]=as.numeric(G[, 1])
  G[, 1]=scale(G[, 1])
  
  
  big = inner_join(comb, G)
  
  x = as.matrix(data.frame(big[,CT],
                           big[,covariates],
                           big[,CT]*big[,'genotype']))

  y=big[,probe]
  a = 0.5
  this_prior = prior[prior$snp == snp & prior$probe == probe,]
  if (nrow(this_prior) == 0){mu_vector = rep(0, length(b_vector))}
  if (nrow(this_prior) > 0){mu_vector = as.numeric(as.vector(this_prior[,paste0('beta_adj_',CT)]))}
  

  oc = EM_IWLS(x, y ,s0 = 0.005, num_CT=length(b_vector), a, b_vector,mu_vector)
  
  if (mean(oc == "singular") == 1){
    oc = EM_IWLS(x, y ,s0 = 0.005, num_CT=length(b_vector), a, 2*b_vector,mu_vector)
    if (mean(oc == "singular") == 1){
      oc = EM_IWLS(x, y ,s0 = 0.005, num_CT=length(b_vector), a, 5*b_vector,mu_vector)
      if (mean(oc == "singular") == 1){
        oc = EM_IWLS(x, y ,s0 = 0.005, num_CT=length(b_vector), a, 20*b_vector,mu_vector)
        if (mean(oc == "singular") == 1){
          oc = EM_IWLS(x, y ,s0 = 0.005, num_CT=length(b_vector), a, 100*b_vector,mu_vector)
          if (mean(oc == "singular") == 1){
            ans = rep(NA, 4*length(b_vector))
          }
        }
      }
    }
  }
  
  #judge if ans exists
  if (exists('ans') == F){   
    betas = oc$beta
    var = oc$var
    J = dim(betas)[1]
    
    eff = betas[(J-length(b_vector)+1):J]
    eff_var = var[(J-length(b_vector)+1):J, (J-length(b_vector)+1):J]
    se = sqrt(diag(eff_var)) 
    z = eff/se
    p = pnorm(abs(z), lower.tail = F)*2

    ans = c(eff,se, z, p)
    
  }
  return(ans)
}


######################################################################################
##########  determine the number of description lines to skip in VCF file  ###########
######################################################################################
get_skip_lines <- function(file_path) {
  lines <- readLines(file_path)
  skip_lines <- which(grepl("^#CHROM", lines)) - 1
  return(skip_lines)
}



########## read data
phen.file <- as.character(args[1])
probe <- as.character(args[2])
geno.file <- as.character(args[3])
W.file <- as.character(args[4])
Covar.file <- as.character(args[5])
Prior.file <- as.character(args[6])
output.file <- as.character(args[7])


skip_lines <- get_skip_lines(geno.file)
vcf  = as.data.frame(fread(geno.file, skip = skip_lines, header=T))
region = vcf[,c('#CHROM','POS','ID','REF','ALT')]
df=as.data.frame(t(vcf))
colnames(df)=df[3,]
df_geno=df[-c(1:9), ]
df_geno[1:5,1:5]

phen = as.data.frame(fread(phen.file, header=T))
W = as.data.frame(fread(W.file, header=T))
Covar = as.data.frame(fread(Covar.file, header=T))
prior = as.data.frame(fread(Prior.file, header=T))
covariates = colnames(Covar)[-1]
CT = colnames(W)[-1]

########## default arguments (optional)
distance <- gsub(x = args[grep(x = args, pattern = "distance=")], pattern = "distance=", replacement = "")
b_vector <- gsub(x = args[grep(x = args, pattern = "b_vector=")], pattern = "b_vector=", replacement = "")

distance = ifelse(length(distance)==0, 500000, as.numeric(distance))
if (length(b_vector) == 0) {
  b_vector <- rep(0.2, length(CT))
} else {
  b_vector <- as.numeric(unlist(strsplit(b_vector, ",")))
}
#c(0.05,0.05,0.2,0.2,0.2)

print(distance)
print(b_vector)

########## combine phen, W, covar
comb = inner_join(W, Covar, by = 'IID')
ind = comb$IID
m    = phen[, ind]
rownames(m) = phen$probe
m    = t(m)
m    = as.data.frame(m)
m$IID = rownames(m)
dim(m)
probe_position = phen[phen$probe == probe, 'BP']
chr = phen[phen$probe == probe, 'CHR']

comb = inner_join(comb, m[,c('IID',probe)], by = 'IID')


########## flip signs for the prior data 
p_CT = paste0('p_',CT)
CT1 = CT[p_CT %in% colnames(prior)]
CT_noPrior = CT[! p_CT %in% colnames(prior)]

prior = inner_join(prior, region[,c('ID','REF','ALT')], by = c('snp'='ID'))
prior[prior$REF.x == prior$ALT.y & prior$ALT.x == prior$REF.y, paste0('beta_adj_',CT1)] = (-1)*prior[prior$REF.x == prior$ALT.y & prior$ALT.x == prior$REF.y, paste0('beta_adj_',CT1)]
prior[,paste0('beta_adj_',CT_noPrior)]=0
prior = prior[,c('probe','snp',paste0('beta_adj_',CT))]

########## find SNPs within 1Mb
region = region[region$`#CHROM` == chr,]
region = region[abs(region$POS - probe_position) < distance,]
dim(region)

########## run the algorithm
oo = lapply(region$ID, run_HBI)
oc = do.call(rbind, oo)
colnames(oc) = c(paste0('BETA_', CT), paste0('SE_', CT), paste0('STAT_', CT), paste0('P_', CT))
oc = as.data.frame(oc)

snp = region$ID
chr = region$`#CHROM`
bp  = region$POS
REF = region$REF
ALT = region$ALT

df = data.frame(probe = probe, probe_position = probe_position, snp, chr, bp, REF, ALT, oc)
colSums(is.na(df))

write.table(df, output.file, quote = F, row.names = F, col.names = T, sep ="\t")

