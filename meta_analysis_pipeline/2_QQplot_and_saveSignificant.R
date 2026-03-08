#################################################################################################
# Untrusted probes were removed in the QC in a platform specific manner. When we meta-analysed their EWAS results, we saw that some of the
# probes that were removed in one platform and not in another appeared in the meta-analysis results. We chose to remove these probes again.
#
# Purpose of script: Removing probes that should have been removed, creating QQ plots and saving top lists
# @author: Anne-Kristin Stavrum
################################################################################################


library(qqman)
library(Cairo)

# LAMBDA FUNCTION
calculate_lambda <- function(pvalues){
  chisq <- qchisq(1-pvalues, 1)
  lambda <- median(chisq)/qchisq(0.5, 1)
  return(lambda)
}

######################################################################################################
# Check cohorts for probes that should have been removed
######################################################################################################
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno.epic2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno.epic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

load(".../Resources/450K_RemoveProbes_ChenProbeIDs.Robj")
remove450 = remove_probes
load(".../Resources/EPICv1_RemoveProbes.Robj")
removeEpic = remove_probes
load(".../Resources/EPICv2_RemoveProbes.RData")
removeEpic2 = remove_probes
removeEpic2 = unique(gsub("_.*$", "", removeEpic2))

load("data_ready_for_meta/Sex_agnostic_BIP.RData")
for(i in 1:length(res)){
  print(names(res)[i])
  cohort = res[[i]]
  print(paste("450K: ",length(intersect(rownames(cohort),remove450))))
  print(paste("epic_v1: ",length(intersect(rownames(cohort),removeEpic))))
  print(paste("epic_v2: ",length(intersect(rownames(cohort),removeEpic2))))
}

######################################################################################################
# View Results Female BIP
######################################################################################################

load("analyses/meta_ewas/backup_meta_analysis_before_removing_probes_that_should_have_been_removed/Meta_BIP_Female.Robj")
res.meta = as.data.frame(res.meta)
res.meta$All_P_Fixed = as.numeric(res.meta$All_P_Fixed)
res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)

if(any(is.na(res.meta$All_P_Random)))
  res.meta = res.meta[-which(is.na(res.meta$All_P_Random)),]
if(length(intersect(rownames(res.meta),remove450))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),remove450)),]
if(length(intersect(rownames(res.meta),removeEpic))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic)),]
if(length(intersect(rownames(res.meta),removeEpic2))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic2)),]

save(res.meta,file = "analyses/meta_ewas/Meta_BIP_Female.Robj")

### QQ plots ##
# Fixed
res.meta = res.meta[with(res.meta, order(All_P_Fixed)),]
lam = calculate_lambda(res.meta$All_P_Fixed)
CairoPNG(file = "meta_analyses/QQ_Female_MetaBIP_FixedEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(res.meta$All_P_Fixed,ylim = c(0,13))
title("Female BIP Meta fixed effects")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

# Random
res.meta = res.meta[with(res.meta, order(All_P_Random)),]
lam = calculate_lambda(na.omit(res.meta$All_P_Random))
CairoPNG(file = "meta_analyses/QQ_Female_MetaBIP_RandomEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Random),ylim = c(0,11))
title("Female BIP Meta random effects")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

# Saving top CpGs with pval < 1e-6
res.meta.top = res.meta[which(res.meta$All_P_Random<1e-6),]
write.csv(res.meta.top,file = "meta_analyses/BIPmeta_Female_RandomTOP1e-6.csv",row.names = T,quote = F)

######################################################################################################
# View Results Male BIP
######################################################################################################
load("analyses/meta_ewas/backup_meta_analysis_before_removing_probes_that_should_have_been_removed/Meta_BIP_Male.Robj")
res.meta = as.data.frame(res.meta)
res.meta$All_P_Fixed = as.numeric(res.meta$All_P_Fixed)
res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)

if(any(is.na(res.meta$All_P_Random)))
  res.meta = res.meta[-which(is.na(res.meta$All_P_Random)),]
if(length(intersect(rownames(res.meta),remove450))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),remove450)),]
if(length(intersect(rownames(res.meta),removeEpic))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic)),]
if(length(intersect(rownames(res.meta),removeEpic2))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic2)),]

save(res.meta,file = "analyses/meta_ewas/Meta_BIP_Male.Robj")

### QQ plots ###
# Fixed
res.meta = res.meta[with(res.meta, order(All_P_Fixed)),]
lam = calculate_lambda(res.meta$All_P_Fixed)
CairoPNG(file = "meta_analyses/QQ_Male_MetaBIP_FixedEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Fixed),ylim = c(0,13))
title("Male BIP Meta fixed effects")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

# Random
res.meta = res.meta[with(res.meta, order(All_P_Random)),]
lam = calculate_lambda(res.meta$All_P_Random)
CairoPNG(file = "meta_analyses/QQ_Male_MetaBIP_RandomEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Random),ylim = c(0,11))
title("Male BIP Meta random effects")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

# Saving top CpGs with pval < 1e-6
res.meta.top = res.meta[which(res.meta$All_P_Random<1e-6),]
write.csv(res.meta.top,file = "meta_analyses/BIPmeta_Male_RandomTOP1e-6.csv",row.names = T,quote = F)

######################################################################################################
# View Results Sex agnostic BIP
######################################################################################################

load("analyses/meta_ewas/backup_meta_analysis_before_removing_probes_that_should_have_been_removed/Meta_BIP_SexAgnostic.Robj")
res.meta = as.data.frame(res.meta)
res.meta$All_P_Fixed = as.numeric(res.meta$All_P_Fixed)
res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)

if(any(is.na(res.meta$All_P_Random)))
  res.meta = res.meta[-which(is.na(res.meta$All_P_Random)),]
if(length(intersect(rownames(res.meta),remove450))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),remove450)),]
if(length(intersect(rownames(res.meta),removeEpic))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic)),]
if(length(intersect(rownames(res.meta),removeEpic2))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic2)),]

save(res.meta,file = "analyses/meta_ewas/Meta_BIP_SexAgnostic.Robj")

### QQ plots ###
# Fixed
lam = calculate_lambda(na.omit(res.meta$All_P_Fixed))
CairoPNG(file = "analyses/meta_ewas/QQ_Sex_Agnostic_MetaBIP_FixedEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Fixed),ylim = c(0,13))
title("Sex agnostic BIP Meta fixed effects")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

# Random
res.meta = res.meta[with(res.meta, order(All_P_Random)),]
lam = calculate_lambda(na.omit(res.meta$All_P_Random))
CairoPNG(file = "analyses/meta_ewas/QQ_Sex_Agnostic_MetaBIP_RandomEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Random),ylim = c(0,11))
title("Sex agnostic BIP Meta random effects")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

# Saving top CpGs with pval < 1e-6

res.meta.top = res.meta[which(res.meta$All_P_Random<1e-6),]
write.csv(res.meta.top,file = "analyses/meta_ewas/BIPmeta_SexAgnostic_RandomTOP1e-6.csv", quote = F,row.names = T)


##############################################################################################
# Preliminary analyses
##############################################################################################
load("meta_analyses/Preliminary_results/Female_MetaBIP.RData")
res.meta = as.data.frame(res.meta)
res.meta$All_P_Fixed = as.numeric(res.meta$All_P_Fixed)
qq(res.meta$All_P_Fixed)
lam = calculate_lambda(na.omit(res.meta$All_P_Fixed))
CairoPNG(file = "meta_analyses/Preliminary_results/QQ_Female_MetaBIP_FixedEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Fixed))
title("Female BIP Meta fixed effects - Preliminary")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)
qq(res.meta$All_P_Random)
lam = calculate_lambda(na.omit(res.meta$All_P_Random))
CairoPNG(file = "meta_analyses/Preliminary_results/QQ_Female_MetaBIP_RandomEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Random))
title("Female BIP Meta random effects - Preliminary")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

load("meta_analyses/Preliminary_results/Male_MetaBIP.RData")
res.meta = as.data.frame(res.meta)
res.meta$All_P_Fixed = as.numeric(res.meta$All_P_Fixed)
qq(res.meta$All_P_Fixed)
lam = calculate_lambda(na.omit(res.meta$All_P_Fixed))
CairoPNG(file = "meta_analyses/Preliminary_results/QQ_Male_MetaBIP_FixedEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Fixed))
title("Male BIP Meta - Preliminary")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)
qq(res.meta$All_P_Random)
lam = calculate_lambda(na.omit(res.meta$All_P_Random))
CairoPNG(file = "meta_analyses/Preliminary_results/QQ_Male_MetaBIP_RandomEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Random))
title("Male BIP Meta random effects - Preliminary")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

load("meta_analyses/Preliminary_results/Sex_Agnostic_MetaBIP.RData")
res.meta = as.data.frame(res.meta)
res.meta$All_P_Fixed = as.numeric(res.meta$All_P_Fixed)
qq(res.meta$All_P_Fixed)
lam = calculate_lambda(na.omit(res.meta$All_P_Fixed))
CairoPNG(file = "meta_analyses/Preliminary_results/QQ_Sex_Agnostic_MetaBIP_FixedEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Fixed))
title("Sex agnostic BIP Meta - Preliminary")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)
qq(res.meta$All_P_Random)
lam = calculate_lambda(na.omit(res.meta$All_P_Random))
CairoPNG(file = "meta_analyses/Preliminary_results/QQ_Sex_Agnostic_MetaBIP_RandomEffect.png", width = 5, height = 5, units = "in", res = 300)
qq(na.omit(res.meta$All_P_Random))
title("Sex agnostic BIP Meta random effects - Preliminary")
mtext(paste("lambda", round(lam, digits=3), sep=" :"), side=1, at=5,line=-1.5)
dev.off()

