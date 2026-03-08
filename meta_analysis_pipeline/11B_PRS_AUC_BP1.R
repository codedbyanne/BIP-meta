###############################################################################################################
## Associations between PC1 of PRS calculated across different p-value thresholds and Case Control in BP1 sub-diagnosis of TOP+
# @author: Anne-Kristin Stavrum
##############################################################################################################

source("scripts/Util_functions.R")

#################################################################################################################################
################ BIPOLAR SUB-TYPE BD1 #######################
#################################################################################################################################

# Load phenotype data
pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))
bpsub = read.csv("Bipolar_Subtype_TOP4_Cases.csv")
pheno$BP_subtype = sapply(pheno$Subj_ID, function(x) ifelse(pheno$Case_Control[match(x,pheno$Subj_ID)]=="Case",bpsub$DIAG[match(x,bpsub$Subj_ID)],"Control"))
bp1 = pheno[pheno$BP_subtype %in% c("BDI","Control"),]

# Load PRS-PC1
pca = as.data.frame(read.csv("analyses/prs/Subj_IDs_TOP4_plink_daner_pgc4_bd_eur_no23andMe_noTOP_HRCfrq.sumstats.gz.all_score_PCA",sep=" "))
pca$SubjID = sapply(pca$IID, function(x) strsplit(x,"_")[[1]][3] )

# SubjID for BP1 samples with PRS
overlap = intersect(bp1$Subj_ID,pca$SubjID)

pca = pca[match(overlap,pca$SubjID),]
bp1 = bp1[match(overlap,bp1$Subj_ID),]

# Load genetic PCs
gPC = as.data.frame(read.csv("analyses/prs/Subj_IDs_TOP4_PCs.txt",sep="\t"))
gPC$SubjID = sapply(gPC$X.IID, function(x) strsplit(x,"_")[[1]][3] )
gPC = gPC[match(overlap,gPC$SubjID),]

# Extract PRS for different thresholds and find PC1 for BP1 samples
prs = as.data.frame(pca[,c(15,3:13)])
colnames(prs)[1]="ID"
PRS_PC1 = prs.pc(prs,"BP1")$data[,c("ID","BP1.prs.pc")]

# Extract MS for different thresholds and find PC1 for BP1 samples
MRS = read.csv("analyses/mrs/MRS_BIP_EPICv1v2.csv")
MRS$ID = gsub("X","",MRS$ID)
MRS = MRS[match(bp1$Basename,MRS$ID),]
MS_PC1 = prs.pc(MRS,"BIP")$data[,c("ID","BIP.prs.pc")]

# Prepare case control phenotypes
TOP4_case_control = factor(bp1$BP_subtype)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

#################################################################################################################################
# Regression
#################################################################################################################################

#Check if PC1 and MS have a positive correlation. Multiply PC1 by -1 if the correlation is negative
if(cor(MS_PC1$BIP.prs.pc,MRS$P5e.08)<0){
  MS_PC1$BIP.prs.pc = MS_PC1$BIP.prs.pc*-1
}

#Check if PC1 and PRS have a positive correlation. Multiply PC1 by -1 if the correlation is negative
if(cor(PRS_PC1$BP1.prs.pc,prs$Pt_5e.08)<0){
  PRS_PC1$BP1.prs.pc = PRS_PC1$BP1.prs.pc*-1
}

# Null model
logReg_null = glm(TOP4_case_control ~ gPC$PC1_AVG+gPC$PC2_AVG+gPC$PC3_AVG+gPC$PC4_AVG+gPC$PC5_AVG, family=binomial) # null model


# PRS-PC1
logReg = glm(TOP4_case_control ~ PRS_PC1$BP1.prs.pc + gPC$PC1_AVG+gPC$PC2_AVG+gPC$PC3_AVG+gPC$PC4_AVG+gPC$PC5_AVG, family=binomial) 
summary(logReg)
# Calculation of NagelkerkeR2 (NKv)
N=length(TOP4_case_control)
LL1 = logLik(logReg)
LL0 = logLik(logReg_null)
CSv = 1-exp((2/N)*(LL0[1]-LL1[1]))
NKv = CSv/(1-exp((2/N)*LL0[1]))

# PRS-PC1 + MS-PC1
logReg = glm(TOP4_case_control ~ MS_PC1$BIP.prs.pc+PRS_PC1$BP1.prs.pc + gPC$PC1_AVG+gPC$PC2_AVG+gPC$PC3_AVG+gPC$PC4_AVG+gPC$PC5_AVG , family=binomial) 
summary(logReg)
# Calculation of NagelkerkeR2 (NKv)
N=length(TOP4_case_control)
LL1 = logLik(logReg)
LL0 = logLik(logReg_null)
CSv = 1-exp((2/N)*(LL0[1]-LL1[1]))
NKv = CSv/(1-exp((2/N)*LL0[1]))

# p-value of the combined effect of PRS-PC1 + MS-PC1
anova(logReg_null,logReg, test = "Chisq")$`Pr(>Chi)`[2]


#################################################################################################################################
# AUC
#################################################################################################################################
library(pROC)

data = data.frame(bp1$Basename,bp1$Subj_ID,TOP4_case_control,pca$PCA_PRS,MS_PC1$BIP.prs.pc,gPC$PC1_AVG,gPC$PC2_AVG,gPC$PC3_AVG,gPC$PC4_AVG,gPC$PC5_AVG)
colnames(data) = c("Basename","Subj_ID","TOP4_case_control","PC1_PRS","PC1_MS","gPC1","gPC2","gPC3","gPC4","gPC5")

aucs_null = vector(mode="numeric",length = 1000)
aucs_prs = vector(mode="numeric",length = 1000)
aucs_ms = vector(mode="numeric",length = 1000)
aucs_prs_ms = vector(mode="numeric",length = 1000)

for(i in 1:1000){
  index = sample(nrow(pheno),nrow(pheno)*0.2)
  train = data[-index,]
  test = data[index,]
  
  logReg = glm(TOP4_case_control ~ gPC1+gPC2+gPC3+gPC4+gPC5, family=binomial, data = train) # association between TOP4 CC and MS_PC1
  predicted = predict(logReg, test, type="response")
  aucs_null[i] = auc(test$TOP4_case_control,predicted)
  
  logReg = glm(TOP4_case_control ~ PC1_PRS + gPC1+gPC2+gPC3+gPC4+gPC5, family=binomial, data = train) # association between TOP4 CC and MS_PC1
  predicted = predict(logReg, test, type="response")
  aucs_prs[i] = auc(test$TOP4_case_control,predicted)
  
  logReg = glm(TOP4_case_control ~ PC1_MS + gPC1+gPC2+gPC3+gPC4+gPC5, family=binomial, data = train) # association between TOP4 CC and MS_PC1
  predicted = predict(logReg, test, type="response")
  aucs_ms[i] = auc(test$TOP4_case_control,predicted)
  
  logReg = glm(TOP4_case_control ~ PC1_PRS+PC1_MS + gPC1+gPC2+gPC3+gPC4+gPC5, family=binomial, data = train) # association between TOP4 CC and MS_PC1
  predicted = predict(logReg, test, type="response")
  aucs_prs_ms[i] = auc(test$TOP4_case_control,predicted)
}

mean(aucs_null)
mean(aucs_prs)
mean(aucs_ms)
mean(aucs_prs_ms)
