###############################################################################################################
## Associations between PC1 of PRS calculated across different p-value thresholds and Case Control in TOP+
# @author: Anne-Kristin Stavrum
##############################################################################################################

source("scripts/Util_functions.R")


#################################################################################################################################
# Load data
#################################################################################################################################

pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))

# get PC1-PRS and format for matching to phenotypes
pca = as.data.frame(read.csv("analyses/prs/Subj_IDs_TOP4_plink_daner_pgc4_bd_eur_no23andMe_noTOP_HRCfrq.sumstats.gz.all_score_PCA",sep=" "))
pca$SubjID = sapply(pca$IID, function(x) strsplit(x,"_")[[1]][3] )
pheno = pheno[match(pca$SubjID,pheno$Subj_ID),]

# Prepare case control phenotypes
TOP4_case_control = factor(pheno$Case_Control)
TOP4_case_control = relevel(TOP4_case_control, ref="Control")

# Get genetic PCs
gPC = as.data.frame(read.csv("analyses/prs/Subj_IDs_TOP4_PCs.txt",sep="\t"))
gPC$SubjID = sapply(gPC$X.IID, function(x) strsplit(x,"_")[[1]][3] )
identical(gPC$SubjID,pca$SubjID)

#################################################################################################################################
# Effect of PRS-BIP
#################################################################################################################################

#Check if PC1 and PRS have a positive correlation. Multiply PC1 by -1 if the correlation is negative
if(cor(pca$PCA_PRS,pca$Pt_5e.08)<0){
  pca$PCA_PRS = pca$PCA_PRS*-1
}

logReg = glm(TOP4_case_control ~ pca$PCA_PRS+gPC$PC1_AVG+gPC$PC2_AVG+gPC$PC3_AVG+gPC$PC4_AVG+gPC$PC5_AVG, family=binomial) # association between TOP4 CC and PRS_PC1
summary(logReg)
logReg_null = glm(TOP4_case_control ~ gPC$PC1_AVG+gPC$PC2_AVG+gPC$PC3_AVG+gPC$PC4_AVG+gPC$PC5_AVG, family=binomial) # null model

# Calculation of NagelkerkeR2 (NKv)
N=length(pheno$Subj_ID)
LL1 = logLik(logReg)
LL0 = logLik(logReg_null)
CSv = 1-exp((2/N)*(LL0[1]-LL1[1]))
NKv = CSv/(1-exp((2/N)*LL0[1]))

# Plot to show association between PC1-PRS and PRS calculated for different thresholds
boxplot(pca$PCA_PRS ~ TOP4_case_control, main="PRS-PC1", xlab = "", ylab="PC1")
par(mfrow=c(2,3))
boxplot(pca$Pt_5e.08 ~ TOP4_case_control, main="PRS", xlab = "", ylab="pval: 5E-8")
boxplot(pca$Pt_1e.07 ~ TOP4_case_control, main="PRS", xlab = "", ylab="pval: 1E-7")
boxplot(pca$Pt_1e.06 ~ TOP4_case_control, main="PRS", xlab = "", ylab="pval: 1E-6")
boxplot(pca$Pt_1e.05 ~ TOP4_case_control, main="PRS", xlab = "", ylab="pval: 1E-5")
boxplot(pca$Pt_0.0001 ~ TOP4_case_control, main="PRS", xlab = "", ylab="pval: 1E-4")
boxplot(pca$Pt_0.001 ~ TOP4_case_control, main="PRS", xlab = "", ylab="pval: 1E-3")

#################################################################################################################################
# Combined effect of PRS-BIP + MS-BIP
#################################################################################################################################
load("analyses/mrs/PC1_BIP.Robj")
MS_PC1$ID = gsub("X","",MS_PC1$ID)
MS_PC1 = MS_PC1[match(pheno$Basename,MS_PC1$ID),]

logReg_null = glm(TOP4_case_control ~ gPC$PC1_AVG+gPC$PC2_AVG+gPC$PC3_AVG+gPC$PC4_AVG+gPC$PC5_AVG, family=binomial) # null model

logReg = glm(TOP4_case_control ~ pca$PCA_PRS+MS_PC1$BIP.prs.pc+gPC$PC1_AVG+gPC$PC2_AVG+gPC$PC3_AVG+gPC$PC4_AVG+gPC$PC5_AVG, family=binomial)
summary(logReg)
# Calculation of NagelkerkeR2 (NKv)
N=length(pheno$Subj_ID)
LL1 = logLik(logReg)
LL0 = logLik(logReg_null)
CSv = 1-exp((2/N)*(LL0[1]-LL1[1]))
NKv = CSv/(1-exp((2/N)*LL0[1]))

# p-value of the combined effect of PRS-PC1 + MS-PC1
anova(logReg_null,logReg, test = "Chisq")$`Pr(>Chi)`[2]

#Plot MS-BIP-PC1 vs PRS-BIP-PC1
plot(MS_PC1$BIP.prs.pc,pca$PCA_PRS)
abline(lm( pca$PCA_PRS ~ MS_PC1$BIP.prs.pc))
summary(lm( pca$PCA_PRS ~ MS_PC1$BIP.prs.pc))

#################################################################################################################################
# AUC
#################################################################################################################################
library(pROC)

data = data.frame(pheno$Basename,pheno$Subj_ID,TOP4_case_control,pca$PCA_PRS,MS_PC1$BIP.prs.pc,gPC$PC1_AVG,gPC$PC2_AVG,gPC$PC3_AVG,gPC$PC4_AVG,gPC$PC5_AVG)
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


