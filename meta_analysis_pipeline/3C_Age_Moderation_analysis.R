#################################################################################################
# Purpose of script: Meta regression analysis with age as the moderator
# @author: Anne-Kristin Stavrum
################################################################################################

#####################################################################
# Loading Libraries
#####################################################################
#library(doParallel)
library(foreach)
library(metafor)

# Meta regression for one probe
metaReg = function(res,probeID,moderator){
  
  meanDiff <- sapply(res, function(x) x[[1]][match(probeID,rownames(x))])
  seDiff <- sapply(res, function(x) x[[2]][match(probeID,rownames(x))])
  
  tryCatch({
    out<- metafor::rma(yi=meanDiff, sei=seDiff, mods=~moderator)

    return(c(sum(!is.na(meanDiff)), out$beta[2,1],out$se[2],out$pval[2],out$ci.lb[2],out$ci.ub[2], out$tau2,out$I2,out$R2))
  },
  error = function(e) print(paste("error",i,sep=": ")))
  
}

ageModeration = function(sumstat, ages, workdir, outfile){
  cpgs = unique(unlist(lapply(sumstat, function(x) rownames(x))))
  
  print("Running the loop")
  
  res.meta.age = foreach(i=1:length(cpgs), .combine=rbind) %do% {
    probeID<-cpgs[i]
    metaReg(res=sumstat,probeID=cpgs[i], moderator = ages)
  }
  
  ## create matrix to store results
  rownames(res.meta.age)<-cpgs
  colnames(res.meta.age)<-c("N_cohorts", "Beta", "SE", "Pval","lower","upper","Tau2","I2","R2")
  
  # Saving result
  print("Saving result")
  save(res.meta.age,file=file.path(workdir,outfile))
  
  xlsx::write.xlsx(as.data.frame(res.meta.age),file = file.path(workdir,paste(strsplit(outfile,".Robj")[[1]][1],"xlsx",sep=".")))
  
  print("Done")
  
}

#####################################################################
# Setting working directory and paths and calling the function to run the analysis
#####################################################################

input1 = "data_ready_for_meta/Data_sig.fixed_SexAgnostic.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/AgeModeration_sig.fixed_SexAgnostic.Robj" # name of output file for meta-analysis
fixed=T
age = c(40.29, 34.56, 50.8, 49.15, 42.33, 31.18, 33.27, 39.33,
        39.54, 35.58, 55.5, 49.58, 42.18, 33.00, 31.82, 36.72,
        16.41, 37.26, 37.6, 36.06)
sig = get(load(input1))
ageModeration(sumstat = sig, ages = age, workdir = ".", outfile = outfile)


input1 = "data_ready_for_meta/Data_sig.random_SexAgnostic.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/AgeModeration_sig.random_SexAgnostic.Robj" # name of output file for meta-analysis
fixed=F
age = c(40.29, 34.56, 50.8, 49.15, 42.33, 31.18, 33.27, 39.33,
        39.54, 35.58, 55.5, 49.58, 42.18, 33.00, 31.82, 36.72,
        16.41, 37.26, 37.6, 36.06)
sig = get(load(input1))
ageModeration(sumstat = sig, ages = age, workdir = ".", outfile = outfile)


input1 = "data_ready_for_meta/Data_sig.fixed_Females.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/AgeModeration_sig.fixed_Females.Robj" # name of output file for meta-analysis
fixed=T
age = c(40.29, 34.56, 50.8, 49.15, 42.33, 31.18, 33.27, 39.33)
sig = get(load(input1))
ageModeration(sumstat = sig, ages = age, workdir = ".", outfile = outfile)


input1 = "data_ready_for_meta/Data_sig.random_Females.RData" #list of summary stats for samples to be metaanalysed
outfile = "analyses/meta_ewas/Heterogeneity_analysis/AgeModeration_sig.random_Females.Robj" # name of output file for meta-analysis
fixed=F
age = c(40.29, 34.56, 50.8, 49.15, 42.33, 31.18, 33.27, 39.33)
sig = get(load(input1))
ageModeration(sumstat = sig, ages = age, workdir = ".", outfile = outfile)


