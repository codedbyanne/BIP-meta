#!/usr/bin/Rscript

#################################################################################################
# Purpose of script: Residualisation of TOP4. 
# This script was run through slurm om the cluster
# @author: Anne-Kristin Stavrum
################################################################################################

#####################################################################
# Loading Libraries
#####################################################################
library(doParallel)
library(foreach)
library(meta)

#####################################################################
# Setting working directory and paths
#####################################################################

args = commandArgs(trailingOnly=TRUE)
workdir = args[1]
input1 = args[2] #file with betas
covariate_file = args[3] # data_frame_with_covariates
no.clusters = as.integer(args[4]) # number of threads for parallisation
#####################################################################
# Loading Data
#####################################################################
print("Loading data")
betas = get(load(input1))
rm(list=ls()[which(ls()!="betas")])
load(covariate_file)
print(ls())
print(paste("No of clusters: ",no.clusters))

#####################################################################
# Preparing 
#####################################################################
print("Setting up parallell")
cl = makeCluster(no.clusters)
registerDoParallel(cl)

#####################################################################
# Residualisation
#####################################################################
residualise = function(CpG_beta, samplesheet){
  
  y = lm(CpG_beta ~ Sex + Age + smoking + cell_PC1 + cell_PC2 +
           Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + 
           Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
           Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
           anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5, data = samplesheet)

  return(resid(y))
  
}

print("Running the loop")
top4_resid = foreach(i=1:length(betas[,1]), .combine=rbind) %dopar% {
#top4_resid = foreach(i=1:10, .combine=rbind) %dopar% {
  residualise(CpG_beta=betas[i,],samplesheet=variables_df)
}

## create matrix to store results
rownames(top4_resid)<-rownames(betas)
colnames(top4_resid)<-colnames(betas)


# Saving result
print("Saving result")
save(top4_resid,file=file.path(workdir,"TOP4_residualised_betas_SexAgeSmoking_5cellPCs_15ctrlPCs_5ancPCs.RData"))

print("Done")


