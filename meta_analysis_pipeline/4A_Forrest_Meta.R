#################################################################################################
# Purpose of script: Forest plot with effect size and standard error of meta-analysis and individual cohorts
# @author: Anne-Kristin Stavrum
################################################################################################

library(ggforestplot)
library(tidyverse)
library(ggforce)


format_data_forest_meta = function(dat){
  for(i in 1:length(dat)){
    dat[[i]]$Cohort = names(dat)[i]
    dat[[i]]$CpG = rownames(dat[[i]])
  }
  df = dat[[1]]
  for(i in 2:length(dat)){
    df = rbind(df,dat[[i]])
  }
  return(df)
}


##########################################
# Forest plot Significant CpGs Female cohorts
##########################################

#SumStat for each cohort
load("meta_analysis/SumStat_significant_BIPmeta/Data_sig.random_Females.RData")

# meta-analysis results
meta_sign = read.csv("meta_analysis/BIPmeta_Female_RandomTOP1e-6_annotated.csv",row.names = 2)
meta_sign = meta_sign[match(rownames(sig.random$FOR2107_F),rownames(meta_sign)),c("All_Effect_Random","All_Effect_SE_Random")]
colnames(meta_sign) = c("logFC","SE")
sig.random[["MetaAnalysis"]] = meta_sign

#Formatting the data
df = format_data_forest_meta(sig.random)
df$Cohort = gsub("TOP4","TOP4+",df$Cohort)

# Group colours
df$group = "SA"
df$group[grep("_M",df$Cohort)]="M"
df$group[grep("_F",df$Cohort)]="F"
df$group[grep("Meta",df$Cohort)]="Meta"
df$group = factor(df$group, levels=c("Meta","F","M","SA"))

#Creating one plot for each cpg
cpgs = unique(df$CpG)
for(cpg in cpgs){
  print(cpg)
  df1 = df[df$CpG %in% cpg,]
  
  # Vertical line at meta-analysis result
  df2 = df1 %>% filter(Cohort=="MetaAnalysis")
  
  pdf(file = paste("meta_analysis/Plots/Forest_Effects_SigRandom_Female_",cpg,".pdf",sep=""), width = 5, height= 4)
  print(forestplot(
    df=df1,
    name=Cohort,
    estimate=logFC,
    se=SE,
    colour = group
    
  ) + 
    geom_vline(data=df2, mapping=aes(xintercept=logFC))+
    ggtitle(cpg)+
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

##########################################
# Forest plot Significant CpGs Sex Agnostic
##########################################

#SumStat for each cohort
load("meta_analysis/SumStat_significant_BIPmeta/Data_sig.random_SexAgnostic.RData")

# meta-analysis results
meta_sign = read.csv("meta_analysis/BIPmeta_SexAgnostic_RandomTOP1e-6.csv",row.names = 1)
meta_sign = meta_sign[match(rownames(sig.random$TOP4_F),rownames(meta_sign)),c("All_Effect_Random","All_Effect_SE_Random")]
colnames(meta_sign) = c("logFC","SE")
sig.random[["MetaAnalysis"]] = meta_sign

#Formatting the data
df = format_data_forest_meta(sig.random)
df$Cohort = gsub("TOP4","TOP4+",df$Cohort)

# Group colours
df$group = "SA"
df$group[grep("_M",df$Cohort)]="M"
df$group[grep("_F",df$Cohort)]="F"
df$group[grep("Meta",df$Cohort)]="Meta"
df$group = factor(df$group, levels=c("Meta","F","M","SA"))

#Creating one plot for each cpg
cpgs = unique(df$CpG)
for(cpg in cpgs){
  print(cpg)
  df1 = df[df$CpG %in% cpg,]
  
  # Vertical line at meta-analysis result
  df2 = df1 %>% filter(Cohort=="MetaAnalysis")
  
  pdf(file = paste("meta_analysis/Plots/Forest_Effects_SigRandom_SexAgnostic_",cpg,".pdf",sep=""), width = 5, height= 4)
  print(forestplot(
    df=df1,
    name=Cohort,
    estimate=logFC,
    se=SE,
    colour = group
    
  ) + 
    geom_vline(data=df2, mapping=aes(xintercept=logFC))+
    ggtitle(cpg)+
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}


