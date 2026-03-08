#################################################################################################
# Purpose of script: Formatting for sensitivity analysis of significant CpGs
# @author: Anne-Kristin Stavrum
################################################################################################

# Females
load("analyses/meta_ewas/Meta_BIP_Female.Robj")
CpGsF_fixed = rownames(res.meta)[res.meta$All_P_Fixed<7.2e-8]
CpGsF_random = rownames(res.meta)[res.meta$All_P_Random<7.2e-8]

load("data_ready_for_meta/Female_BIP.RData")

sig.fixed = lapply(res, function(x){
  overlap=intersect(CpGsF_fixed,rownames(x))
  index = match(overlap,rownames(x))
  data.frame(row.names = rownames(x)[index], logFC=x$logFC[index],SE=x$SE[index])
} )
save(sig.fixed,file="data_ready_for_meta/Data_sig.fixed_Females.RData")

sig.random = lapply(res, function(x){
  overlap=intersect(CpGsF_random,rownames(x))
  index = match(overlap,rownames(x))
  data.frame(row.names = rownames(x)[index], logFC=x$logFC[index],SE=x$SE[index])
} )
save(sig.random,file="data_ready_for_meta/Data_sig.random_Females.RData")

# Males
load("analyses/meta_ewas/Meta_BIP_Male.Robj")
CpGsF_fixed = rownames(res.meta)[res.meta$All_P_Fixed<7.2e-8]
CpGsF_random = rownames(res.meta)[res.meta$All_P_Random<7.2e-8]

load("data_ready_for_meta/Male_BIP.RData")

sig.fixed = lapply(res, function(x){
  overlap=intersect(CpGsF_fixed,rownames(x))
  index = match(overlap,rownames(x))
  data.frame(row.names = rownames(x)[index], logFC=x$logFC[index],SE=x$SE[index])
} )
save(sig.fixed,file="data_ready_for_meta/Data_sig.fixed_Males.RData")


## nothing significant for random effect for males ##


# Sex Agnostic
load("analyses/meta_ewas/Meta_BIP_SexAgnostic.Robj")
CpGsF_fixed = rownames(res.meta)[res.meta$All_P_Fixed<7.2e-8]
CpGsF_random = rownames(res.meta)[res.meta$All_P_Random<7.2e-8]

load("data_ready_for_meta/Sex_agnostic_BIP.RData")

sig.fixed = lapply(res, function(x){
  overlap=intersect(CpGsF_fixed,rownames(x))
  index = match(overlap,rownames(x))
  data.frame(row.names = rownames(x)[index], logFC=x$logFC[index],SE=x$SE[index])
} )
save(sig.fixed,file="data_ready_for_meta/Data_sig.fixed_SexAgnostic.RData")

sig.random = lapply(res, function(x){
  overlap=intersect(CpGsF_random,rownames(x))
  index = match(overlap,rownames(x))
  data.frame(row.names = rownames(x)[index], logFC=x$logFC[index],SE=x$SE[index])
} )
save(sig.random,file="data_ready_for_meta/Data_sig.random_SexAgnostic.RData")



