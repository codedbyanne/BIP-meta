###############################################################################################################
## Purpose of script: Associations analysis between PC1 of MS calculated accross different p-value thresholds and diagnoses
## Diagnoses: BIP, PTSD, SCZ and MDD
# @author: Anne-Kristin Stavrum
###############################################################################################################

source("scripts/Util_functions.R")

#Load phenotype
pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))
pheno = pheno[,c("Basename","Case_Control")]
head(pheno)
pheno$ID = paste("X",pheno$Basename,sep="")
colnames(pheno) = c("Basename", "Case_Control","ID")



################ BIPOLAR #######################
MRS = read.csv("analyses/mrs/MRS_BIP_EPICv1v2.csv")
MS_PC1 = prs.pc(MRS,"BIP")$data[,c("ID","BIP.prs.pc")]
save(MS_PC1, file="analyses/mrs/PC1_BIP.Robj")

identical(MS_PC1$ID,pheno$ID)
TOP4_case_control = factor(pheno$Case_Control)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

if(cor(MS_PC1$BIP.prs.pc,MRS$P5e.10)<0){
  pca$PCA_PRS = pca$PCA_PRS*-1
}

logReg = glm(TOP4_case_control ~ MS_PC1$BIP.prs.pc, family=binomial) # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)


### Figure MS vs Case Control
colours = pheno$Case_Control
colours[colours=="Case"]="Blue"
colours[colours=="Control"]="Green"

par(mfrow=c(2,4))
boxplot(MRS$P5e.10 ~ pheno$Case_Control, main="P<10e-10", xlab = "")
boxplot(MRS$P5e.09 ~ pheno$Case_Control, main="P<10e-09", xlab = "")
boxplot(MRS$P5e.08 ~ pheno$Case_Control, main="P<10e-08", xlab = "")
boxplot(MRS$P5e.07 ~ pheno$Case_Control, main="P<10e-07", xlab = "")
boxplot(MRS$P5e.06 ~ pheno$Case_Control, main="P<10e-06", xlab = "")
boxplot(MRS$P5e.05 ~ pheno$Case_Control, main="P<10e-05", xlab = "")
boxplot(MRS$P5e.04 ~ pheno$Case_Control, main="P<10e-04", xlab = "")
boxplot(MS_PC1$BIP.prs.pc ~ pheno$Case_Control, main="MS - PC1", xlab = "")


################ Schisophrenia #######################
MRS = read.csv("analyses/mrs/MRS_SCZ_EPICv1v2.csv")
MS_PC1 = prs.pc(MRS,"SCZ")$data[,c("ID","SCZ.prs.pc")]
save(MS_PC1, file="analyses/mrs/PC1_SCZ.Robj")

identical(MS_PC1$ID,pheno$ID)
TOP4_case_control = factor(pheno$Case_Control)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

if(cor(MS_PC1$SCZ.prs.pc,MRS$P5e.08)<0){
  pca$PCA_PRS = pca$PCA_PRS*-1
}


logReg = glm(TOP4_case_control ~ MS_PC1$SCZ.prs.pc, family=binomial) # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

################ PTSD #######################
MRS = read.csv("analyses/mrs/MRS_PTSD_EPICv1v2.csv")
MS_PC1 = prs.pc(MRS,"PTSD")$data[,c("ID","PTSD.prs.pc")]
save(MS_PC1, file="analyses/mrs/PC1_PTSD.Robj")

identical(MS_PC1$ID,pheno$ID)
TOP4_case_control = factor(pheno$Case_Control)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

if(cor(MS_PC1$PTSD.prs.pc,MRS$P5e.09)<0){
  pca$PCA_PRS = pca$PCA_PRS*-1
}

logReg = glm(TOP4_case_control ~ MS_PC1$PTSD.prs.pc, family=binomial) # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)


################ MDD #######################
MRS = read.csv("analyses/mrs/MRS_MDD_EPICv1v2.csv")
MRS = MRS[,1:4] # The rest of the columns are identical to column 4, so removing them
MS_PC1 = prs.pc(MRS,"MDD")$data[,c("ID","MDD.prs.pc")]
save(MS_PC1, file="analyses/mrs/PC1_MDD.Robj")

identical(MS_PC1$ID,pheno$ID)
TOP4_case_control = factor(pheno$Case_Control,)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

logReg = glm(TOP4_case_control ~ MS_PC1$MDD.prs.pc, family=binomial) # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)



par(mfrow=c(2,2))
load("analyses/mrs/PC1_BIP.Robj")
boxplot(MS_PC1$BIP.prs.pc ~ TOP4_case_control, main="BIP", xlab = "")
load("analyses/mrs/PC1_SCZ.Robj")
boxplot(MS_PC1$SCZ.prs.pc ~ TOP4_case_control, main="SCZ", xlab = "")
load("analyses/mrs/PC1_MDD.Robj")
boxplot(MS_PC1$MDD.prs.pc ~ TOP4_case_control, main="MDD", xlab = "")
load("analyses/mrs/PC1_PTSD.Robj")
boxplot(MS_PC1$PTSD.prs.pc ~ TOP4_case_control, main="PTSD", xlab = "")



################ BIPOLAR SUB-TYPE BD1 #######################
pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))

bpsub = read.csv("Bipolar_Subtype_TOP4_Cases.csv")
pheno$BP_subtype = sapply(pheno$Subj_ID, function(x) ifelse(pheno$Case_Control[match(x,pheno$Subj_ID)]=="Case",bpsub$DIAG[match(x,bpsub$Subj_ID)],"Control"))

bp1 = pheno[pheno$BP_subtype %in% c("BDI","Control"),]

MRS = read.csv("analyses/mrs/MRS_BIP_EPICv1v2.csv")
MRS$ID = gsub("X","",MRS$ID)
MRS = MRS[match(bp1$Basename,MRS$ID),]
MS_PC1 = prs.pc(MRS,"BIP")$data[,c("ID","BIP.prs.pc")]

TOP4_case_control = factor(bp1$BP_subtype)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

if(cor(MS_PC1$BIP.prs.pc,MRS$P5e.8)<0){
  pca$PCA_PRS = pca$PCA_PRS*-1
}

logReg = glm(TOP4_case_control ~ MS_PC1$BIP.prs.pc, family=binomial) # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

################ BIPOLAR SUB-TYPE BD2 #######################
pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))

bpsub = read.csv("Bipolar_Subtype_TOP4_Cases.csv")
pheno$BP_subtype = sapply(pheno$Subj_ID, function(x) ifelse(pheno$Case_Control[match(x,pheno$Subj_ID)]=="Case",bpsub$DIAG[match(x,bpsub$Subj_ID)],"Control"))

bp1 = pheno[pheno$BP_subtype %in% c("BDII","Control"),]

MRS = read.csv("analyses/mrs/MRS_BIP_EPICv1v2.csv")
MRS$ID = gsub("X","",MRS$ID)
MRS = MRS[match(bp1$Basename,MRS$ID),]
MS_PC1 = prs.pc(MRS,"BIP")$data[,c("ID","BIP.prs.pc")]

TOP4_case_control = factor(bp1$BP_subtype)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

if(cor(MS_PC1$BIP.prs.pc,MRS$P5e.8)<0){
  pca$PCA_PRS = pca$PCA_PRS*-1
}

logReg = glm(TOP4_case_control ~ MS_PC1$BIP.prs.pc, family=binomial) # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

################ BIPOLAR SUB-TYPE BDNOS #######################
pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))

bpsub = read.csv("Bipolar_Subtype_TOP4_Cases.csv")
pheno$BP_subtype = sapply(pheno$Subj_ID, function(x) ifelse(pheno$Case_Control[match(x,pheno$Subj_ID)]=="Case",bpsub$DIAG[match(x,bpsub$Subj_ID)],"Control"))

bp1 = pheno[pheno$BP_subtype %in% c("BDNOS","Control"),]

MRS = read.csv("analyses/mrs/MRS_BIP_EPICv1v2.csv")
MRS$ID = gsub("X","",MRS$ID)
MRS = MRS[match(bp1$Basename,MRS$ID),]
MS_PC1 = prs.pc(MRS,"BIP")$data[,c("ID","BIP.prs.pc")]

TOP4_case_control = factor(bp1$BP_subtype)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

if(cor(MS_PC1$BIP.prs.pc,MRS$P5e.8)<0){
  pca$PCA_PRS = pca$PCA_PRS*-1
}

logReg = glm(TOP4_case_control ~ MS_PC1$BIP.prs.pc, family=binomial) # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

