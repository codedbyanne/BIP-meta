#################################################################################################
# Purpose of script: Creating overview table of heterogeneity
# @author: Anne-Kristin Stavrum
################################################################################################

format_data_forest_LOO = function(dat){
  for(i in 1:length(dat)){
    dat[[i]]$CpG = names(dat)[i]
  }
  df = dat[[1]]
  for(i in 2:length(dat)){
    df = rbind(df,dat[[i]])
  }
  return(df)
}

################################################################################################################
# Creating overview table of heterogeneity: which cpg is significant by leaving out a particular cohort
################################################################################################################

load("meta_analysis/Sensitivity_SignDMPs_BIPmeta/LOOmeta_sig.random_SexAgnostic.Robj")
df.loo = format_data_forest_LOO(sensitivity)
df.loo$Cohort = gsub("TOP4","TOP4+",df.loo$Cohort)

loo.cohort=vector(mode = "numeric")
cohorts = unique(df.loo$Cohort)
cohorts = cohorts[-length(cohorts)]
for(c in cohorts){
  df.cohort = df.loo[which(c == df.loo$Cohort),]
  
  loo.cohort[c] = length(which(df.cohort$All_P_Random<=5.5e-8))
  print(c)
  print(df.cohort$CpG[which(df.cohort$All_P_Random>5.5e-8)])
  
}

loo.cohort

loo.cpg=vector(mode = "numeric")
cpgs = unique(df.loo$CpG)
for(cpg in cpgs){
  df.cpg = df.loo[which(cpg == df.loo$CpG),]
  
  loo.cpg[cpg] = length(which(df.cpg$All_P_Random<5.5e-8))-1
  print(c)
  print(df.cpg$CpG[which(df.cpg$All_P_Random>5.5e-8)])
  
}
loo.cpg

loo.res = matrix(nrow=length(cpgs),ncol = length(cohorts))
rownames(loo.res) = cpgs;colnames(loo.res)=cohorts
for(cpg in cpgs){
  df.cpg = df.loo[which(cpg == df.loo$CpG),]
  loo.res[cpg,] = ifelse(df.cpg$All_P_Random[match(cohorts,df.cpg$Cohort)]<5.5e-8,1,0)
}
loo.res = as.data.frame(loo.res)
identical(df.loo[df.loo$Cohort=="Meta",]$CpG, rownames(loo.res))
loo.res$meta_pval =df.loo$All_P_Random[df.loo$Cohort=="Meta"]

loo.res$CpG = rownames(loo.res)
writexl::write_xlsx(loo.res,path ="meta_analysis/Sensitivity_SignDMPs_BIPmeta/LOOmetaAnalysis_CpG_vs_Cohort_withMetaPval.xlsx")




