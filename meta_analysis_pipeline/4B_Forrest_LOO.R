#################################################################################################
# Purpose of script: Forest plot with effect size meta-analysis and individual cohorts when leaving out one cohort. The named cohort is left out.
# @author: Anne-Kristin Stavrum
################################################################################################

library(ggforestplot)
library(tidyverse)
library(ggforce)

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

##########################################
# Forest plot showing effect of leaving out one sample (LOO) on significant CpGs FEMALE cohorts
##########################################
load("meta_analysis/Sensitivity_SignDMPs_BIPmeta/LOOmeta_sig.random_Females.Robj")
df = format_data_forest_LOO(sensitivity)
df$Cohort = gsub("TOP4","TOP4+",df$Cohort)

# Grouping colours
df$group = "SA"
df$group[grep("_M",df$Cohort)]="M"
df$group[grep("_F",df$Cohort)]="F"
df$group[grep("Meta",df$Cohort)]="Meta"
df$group = factor(df$group, levels=c("Meta","F","M","SA"))
colnames(df)[2] = "logFC"

#Creating one plot for each cpg
cpgs = unique(df$CpG)
for(cpg in cpgs){
  print(cpg)
  df1 = df[df$CpG %in% cpg,]
  
  # Vertical line at meta-analysis result
  df2 = df1 %>% filter(Cohort=="Meta")
  
  pdf(file = paste("meta_analysis/Plots/Forest_LOOmeta_Random_Female_",cpg,".pdf",sep=""), width = 7, height= 4)
  print(forestplot(
    df=df1,
    name=Cohort,
    estimate=logFC,
    se=0,
    pvalue=All_P_Random,
    psignif=7.2e-08,
    colour = group
    
  ) + 
    geom_vline(data=df2, mapping=aes(xintercept=logFC))+
    ggtitle(cpg)+
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}



##########################################
# Forest plot showing effect of leaving out one sample (LOO) on significant CpGs Female, male and Sex agnostic cohorts
##########################################

load("meta_analysis/Sensitivity_SignDMPs_BIPmeta/LOOmeta_sig.random_SexAgnostic.Robj")
df = format_data_forest_LOO(sensitivity)
df$Cohort = gsub("TOP4","TOP4+",df$Cohort)

# Grouping colours
df$group = "SA"
df$group[grep("_M",df$Cohort)]="M"
df$group[grep("_F",df$Cohort)]="F"
df$group[grep("Meta",df$Cohort)]="Meta"
df$group = factor(df$group, levels=c("Meta","F","M","SA"))

# Renaming to be consistent with other plots
colnames(df)[grep("mean",colnames(df))]="logFC"
df$Cohort[grep("Meta",df$Cohort)]="MetaAnalysis"

#Creating one plot for each cpg
cpgs = unique(df$CpG)
for(cpg in cpgs){
  print(cpg)
  df1 = df[df$CpG %in% cpg,]
  
  # Vertical line at meta-analysis result
  df2 = df1 %>% filter(Cohort=="MetaAnalysis")
  
  pdf(file = paste("meta_analysis/Plots/Forest_LOOmeta_Random_SexAgnostic_",cpg,".pdf",sep=""), width = 7, height= 4)
  print(forestplot(
    df=df1,
    name=Cohort,
    estimate=logFC,
    se=0,
    pvalue=All_P_Random,
    psignif=7.2e-08,
    colour = group
    
  ) + 
    geom_vline(data=df2, mapping=aes(xintercept=logFC))+
    ggtitle(cpg)+
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

