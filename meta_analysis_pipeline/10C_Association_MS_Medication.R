
###############################################################################################################
## Associations between PC1 of MS calculated accross different p-value thresholds and medication
## Medication groups tested: AP, AD, AE and lithium
# @author: Anne-Kristin Stavrum
###############################################################################################################

top4 = readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx")

table(top4$Case_Control,top4$Project)

medication = haven::read_sav(".../medic.sav")

medication = medication[medication$Subj_ID %in% top4$Subj_ID,]

medication$AP1 = factor(rio::factorize(medication$AP1))
medication$AP2 = factor(rio::factorize(medication$AP2))
medication$AP3 = factor(rio::factorize(medication$AP3))

medication$AD1 = factor(rio::factorize(medication$AD1_name))
medication$AD2 = factor(rio::factorize(medication$AD2_name))

medication$AE1 = factor(rio::factorize(medication$AE1_name))
medication$AE2 = factor(rio::factorize(medication$AE2_name))
medication$AE3 = factor(rio::factorize(medication$AE3_name))

medication$LIT = factor(rio::factorize(medication$LIT_name))

print("number of samples with medication information in TOP4")
length(intersect(medication$Subj_ID,top4$Subj_ID))

table(c(medication$AP1,medication$AP2,medication$AP3))
table(c(medication$AD1,medication$AD2))
table(c(medication$AE1,medication$AE2,medication$AE3))
table(medication$LIT)

top4 = top4[match(medication$Subj_ID,top4$Subj_ID),]

top4$AP = 0
top4$AP[which(!is.na(medication$AP1))]=1
top4$AD = 0
top4$AD[which(!is.na(medication$AD1))]=1
top4$AE = 0
top4$AE[which(!is.na(medication$AE1))]=1
top4$Lithium = 0
top4$Lithium[which(!is.na(medication$LIT))]=1

table(top4$AP)
table(top4$AD)
table(top4$AE)
table(top4$Lithium)


load("analyses/mrs/PC1_BIP.Robj")
MS_PC1$ID = gsub("X","",MS_PC1$ID)
MS_PC1 = MS_PC1[match(top4$Basename,MS_PC1$ID),]


logReg = glm(factor(top4$AP) ~ MS_PC1$BIP.prs.pc, family="binomial") # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

logReg = glm(factor(top4$AD) ~ MS_PC1$BIP.prs.pc, family="binomial") # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

logReg = glm(factor(top4$AE) ~ MS_PC1$BIP.prs.pc, family="binomial") # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

logReg = glm(factor(top4$Lithium) ~ MS_PC1$BIP.prs.pc, family="binomial") # association between TOP4 CC and MS_PC1
summary(logReg)
fmsb::NagelkerkeR2(logReg)

