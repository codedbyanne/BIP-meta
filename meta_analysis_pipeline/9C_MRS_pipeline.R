#################################################################################################
# Purpose of script: Calculate Methylation Scores
# This script is adapted from a script available from https://pubmed.ncbi.nlm.nih.gov/36908043/
# https://github.com/jche453/Pruning-Thresholding-MRS.git
################################################################################################

#####################################
##############STEP ONE###############
#####################################

#####################################
# Load phenotypes file
pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))

##########################################################
###Run cmr function to generate Co-methylation regions####
##########################################################
### Generate co-methylated regions ###
# install.packages("comeback_0.1.0.tar.gz", repos = NULL, type="source")
#comeback_0.1.0.tar.gz can be found at https://bitbucket.org/flopflip/comeback.
#library(comeback)
source("scripts/MRS_func.R")
source("scripts/comback_mod.R")
#load("analyses/mrs/TOP4_residualised_betas_SexAgeSmoking_2cellPCs_15ctrlPCs_5ancPCs.RData")
#load("analyses/mrs/TOP4_residualised_Mvalues_SexAgeSmoking_2cellPCs_15ctrlPCs_5ancPCs.RData")
load("analyses/mrs/TOP4_residualised_Mvalues_SexAgeSmoking_5cellPCs_15ctrlPCs_5ancPCs.RData")

#Find co-methylated regions
res = t(top4_resid) #column are CpG sites, rows are for samples
cmrs <- cmr(Mdata = res, Iarray="EPICv1v2", corlo = 0.3)
save(cmrs, file="analyses/mrs/CMR_EPICv1v2.TOP4.RData")

# To generate Co-methylation regions dataframe for later analysis, 
# you can either input cmrs calculated from your own dataset,
CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, Overlap = F)

#CoMeRegion is a matrix that assigned a unique number to each co-methylation region
save(CoMeRegion, file = "analyses/mrs/CoMeRegion_EPICv1v2_TOP4.rda")

#####################################
##############STEP TWO###############
#####################################
### Calculate MRS
library(dplyr)    

# Load real phenotype and methylation data for individuals we want to calculated different MS for
#DNAm = get(load("analyses/mrs/TOP4_residualised_Mvalues_SexAgeSmoking_2cellPCs_15ctrlPCs_5ancPCs.RData"))
DNAm = get(load("analyses/mrs/TOP4_residualised_Mvalues_SexAgeSmoking_5cellPCs_15ctrlPCs_5ancPCs.RData"))
rm(top4_resid);gc()
DNAm = t(DNAm)
#DNAm = res
DNAm[1:6,1:6] #check DNAm file

########################################################################################################################
### Calculation of MS for BIP
########################################################################################################################

ss <- get(load("analyses/mrs/sumstat/Meta_Sex_agnostic_BIP_LeavingOutTOP4.Robj"))
ssn = as.data.frame(apply(ss,2,as.numeric))
rownames(ssn) = rownames(ss)
ssn = ssn[which(!is.na(ssn$All_Effect_Fixed)),]
ssn$Marker = rownames(ssn)
rownames(ssn) = 1:length(ssn$Marker)
SS = ssn[,c(12,5,6,7)]
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)
rm(ss,ssn)
#Get the smallest p-value
minpvalue = min(na.omit(SS$Pvalue[SS$Pvalue != 0]))
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("CoMeRegion_EPICv1v2_TOP4.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(DNAm, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "analyses/mrs/MRS_BIP_pvalueinfo_EPICv1v2.csv", row.names = F)
write.csv(MRS$MRS, "analyses/mrs/MRS_BIP_EPICv1v2.csv", row.names = F)

########################################################################################################################
### Calculation of MS for PTSD 
########################################################################################################################
ss = readxl::read_xlsx("analyses/mrs/sumstat/PTSD_Meta_analysis_23_Cohorts_13073_2024_1417_MOESM2_ESM.xlsx")
SS = ss[,c(1,5,6,8)]
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)

#Get the smallest p-value
minpvalue = min(na.omit(SS$Pvalue[SS$Pvalue != 0]))
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("CoMeRegion_EPICv1v2_TOP4.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(DNAm, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "analyses/mrs/MRS_PTSD_pvalueinfo_EPICv1v2.csv", row.names = F)
write.csv(MRS$MRS, "analyses/mrs/MRS_PTSD_EPICv1v2.csv", row.names = F)


########################################################################################################################
### Calculation of MS for Schizophrenia
########################################################################################################################
ss = read.csv("analyses/mrs/sumstat/SCZ_Meta2_PC3_Imput.EWAS_Meta_Analysis_Results.csv")
SS = ss[,c(1,6,7,8)]
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)

#Get the smallest p-value
minpvalue = min(na.omit(SS$Pvalue[SS$Pvalue != 0]))
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("CoMeRegion_EPICv1v2_TOP4.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(DNAm, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "analyses/mrs/MRS_SCZ_pvalueinfo_EPICv1v2.csv", row.names = F)
write.csv(MRS$MRS, "analyses/mrs/MRS_SCZ_EPICv1v2.csv", row.names = F)


########################################################################################################################
### Calculation of MS for MDD
########################################################################################################################
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
mdd = readxl::read_xlsx("analyses/mrs/ScottishMDD_Meta_Table_2.xlsx")
mdd$cg = sapply(1:length(mdd$CHR), function(x) anno$Name[which(anno$chr==paste("chr", mdd$CHR[x],sep="") & anno$pos==mdd$BP[x])])

SS = mdd[,c(13,3,4,5)]
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)

#Get the smallest p-value
minpvalue = min(na.omit(SS$Pvalue[SS$Pvalue != 0]))
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("CoMeRegion_EPICv1v2_TOP4.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(DNAm, SS, Pthred, CoMeRegion, CoMeBack = F, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "analyses/mrs/MRS_MDD_pvalueinfo_EPICv1v2.csv", row.names = F)
write.csv(MRS$MRS, "analyses/mrs/MRS_MDD_EPICv1v2.csv", row.names = F)



