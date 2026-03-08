#!/usr/bin/Rscript

#################################################################################################
# Purpose of script: Meta-analysis. This script was run through slurm on the cluster
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
input1 = args[2] #file with summary stats for cohorts to be metaanalysed
outfile = args[3] # name of output file for meta-analysis
nProbes = as.integer(args[4]) # number cohorts that needs to have result for a particular probe
no.clusters = as.integer(args[5]) # number of threads for parallisation
#####################################################################
# Loading Data
#####################################################################
print("Loading data")
load(input1)
#####################################################################
# Preparing 
#####################################################################
print("Setting up parallell")
cl = makeCluster(no.clusters)
registerDoParallel(cl)

#####################################################################
# Meta analysis
#####################################################################
print("Setting up analysis")
n_cohorts = length(res)
probes = unlist(sapply(res,rownames))
probes<-table(probes)
probes<-probes[which(probes >= nProbes)]

# Function for meta-analysing one probe
metaCpG = function(row,probeID){
  
	meanDiff <- sapply(res, function(x) x[[1]][match(probeID,rownames(x))])
  seDiff <- sapply(res, function(x) x[[2]][match(probeID,rownames(x))])
  

  tryCatch({
    out<- meta::metagen(meanDiff, seDiff)
    return(c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q)))
  },
  error = function(e) print(paste("error",i,sep=": ")))

}

print("Running the parallell loop for all probes to be analysed")
res.meta = foreach(i=1:length(probes), .combine=rbind) %dopar% {
  probeID<-names(probes)[i]
  metaCpG(row=i,probeID=probeID)
}

## setting rownames and column names of result matrix
rownames(res.meta)<-names(probes)
colnames(res.meta)<-c("N_cohorts", "All_Effect_Fixed", "All_Effect_SE_Fixed", "All_P_Fixed", "All_Effect_Random", "All_Effect_SE_Random","All_P_Random", "All_tau", "All_I2", "All_Q", "All_Het P")


# Saving result
print("Saving result")
save(res.meta,file=file.path(workdir,outfile))

print("Done")
