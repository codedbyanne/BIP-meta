#################################################################################################
# Purpose of script: Associations between PC1 of MS calculated across different p-value thresholds and smoking
# @author: Anne-Kristin Stavrum
################################################################################################

source("scripts/Util_functions.R")

#############################################
# Association smoking
#############################################

betas = get(load(".../phaseIV/Bipolar_pipeline/output/RData/beta_final_autosomes.RData"));rm(beta)
smokescores = calc_smokescore(betas)
save(smokescores,file="analyses/mrs/smokescores.Robj")

pdf("analyses/mrs/plots/BoxplotSmokingScoreSplitByCaseControl.pdf")
boxplot(smokescores ~ factor(pheno$Case_Control), ylab = "DNAm smoking score", main = "Current smoking status")
dev.off()
pdf("analyses/mrs/plots/BoxplotSmokingScoreSplitBySmokerNonsmoker_.pdf")
boxplot(smokescores[-which(pheno$royker_snuser_daglig==999)] ~ factor(pheno$royker_snuser_daglig[-which(pheno$royker_snuser_daglig==999)]), ylab = "DNAm smoking score", main = "", names = c("Non-smoker", "Smoker"))
dev.off()

load("analyses/mrs/PC1_BIP.Robj")
load("analyses/mrs/smokescores.Robj")

reg = glm(smokescores ~ MS_PC1$BIP.prs.pc, family="gaussian") # association between TOP4 CC and MS_PC1
summary(reg)
fmsb::NagelkerkeR2(reg)






