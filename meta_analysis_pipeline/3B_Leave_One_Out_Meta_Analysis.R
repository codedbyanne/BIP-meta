#################################################################################################
# Purpose of script: Sensitivity analysis of significant CpGs
# @author: Anne-Kristin Stavrum
################################################################################################

#####################################################################
# Loading Libraries
#####################################################################
library(foreach)
library(meta)
library(forestplot)
library(dplyr)
library(Cairo)

#####################################################################
# Functions
#####################################################################

# Function for meta-analysing one probe
metaCpG = function(res_i,probeID,fixed){
  
  meanDiff <- sapply(res_i, function(x) x[[1]][match(probeID,rownames(x))])
  seDiff <- sapply(res_i, function(x) x[[2]][match(probeID,rownames(x))])
  
  tryCatch({
    out<- meta::metagen(meanDiff, seDiff)
    if(fixed){
      ci.lower = out$TE.fixed-out$seTE.fixed*qnorm(0.975)
      ci.upper = out$TE.fixed+out$seTE.fixed*qnorm(0.975)
      return(c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed,ci.lower,ci.upper))
    }else{
      ci.lower = out$TE.random-out$seTE.random*qnorm(0.975)
      ci.upper = out$TE.random+out$seTE.random*qnorm(0.975)
      return(c(sum(!is.na(meanDiff)), out$TE.random, out$seTE.random,out$pval.random,ci.lower,ci.upper))
    }
    
  },
  error = function(e) print(paste("error",i,sep=": ")))
  
}

#Iteratively remove one cohort and calling metaCpG for each significan CpG
LOOmeta = function(input, workdir, outfile, fixed){
  
  # Loading Data
  print("Loading data")
  res = get(load(input))
  
  # extracting cpgs for sensitivity analysis 
  cpgs = unique(unlist(lapply(res, rownames)))

  sensitivity = list()
  
  for(cpg in cpgs){
    
    res.meta = foreach(i=1:length(res), .combine=rbind) %do% {
      res_i = res[-i] # removing one cohort
      metaCpG(res_i=res_i,probeID=cpg,fixed=fixed)
    }
    if(fixed){
      colnames(res.meta)<-c("N_cohorts", "mean", "All_Effect_SE_Fixed", "All_P_Fixed","lower","upper")
    }else{
      colnames(res.meta)<-c("N_cohorts", "mean", "All_Effect_SE_Random", "All_P_Random","lower","upper")
    }
    
    rownames(res.meta)<-names(res)
    
    # Repeating the full meta, to include this result (could also have been extracted from original meta-analysis)
    Meta = metaCpG(res_i=res,probeID=cpg,fixed=fixed)
    res.meta = rbind(res.meta,Meta)
    
    res.meta = data.frame(res.meta)
    res.meta$Cohort = rownames(res.meta)
    sensitivity[[cpg]] = res.meta
  }
  
  # Saving result
  print("Saving result")
  save(sensitivity,file=file.path(workdir,outfile))
  
  
  print("Done")
  
  
}
#####################################################################
# Setting working directory and paths and calling LOOmeta
#####################################################################

workdir = "."
input1 = "data_ready_for_meta/Data_sig.random_Females.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/LOOmeta_sig.random_Females.Robj" # name of output file for meta-analysis
fixed=F
LOOmeta(input = input1, workdir = workdir, outfile = outfile, fixed = fixed)

workdir = "."
input1 = "data_ready_for_meta/Data_sig.fixed_Females.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/LOOmeta_sig.fixed_Females.Robj" # name of output file for meta-analysis
fixed=T
LOOmeta(input = input1, workdir = workdir, outfile = outfile, fixed = fixed)


workdir = "."
input1 = "data_ready_for_meta/Data_sig.fixed_Males.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/LOOmeta_sig.fixed_Males.Robj" # name of output file for meta-analysis
fixed=T
LOOmeta(input = input1, workdir = workdir, outfile = outfile, fixed = fixed)


workdir = "."
input1 = "data_ready_for_meta/Data_sig.fixed_SexAgnostic.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/LOOmeta_sig.fixed_SexAgnostic.Robj" # name of output file for meta-analysis
fixed=T
LOOmeta(input = input1, workdir = workdir, outfile = outfile, fixed = fixed)


workdir = "."
input1 = "data_ready_for_meta/Data_sig.random_SexAgnostic.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/LOOmeta_sig.random_SexAgnostic.Robj" # name of output file for meta-analysis
fixed=F
LOOmeta(input = input1, workdir = workdir, outfile = outfile, fixed = fixed)


