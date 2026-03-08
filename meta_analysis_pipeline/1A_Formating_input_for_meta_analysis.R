#################################################################################################
# Purpose of script: Load summary statistics from individual studies and put in a list formatted for meta-analysis
# @author: Anne-Kristin Stavrum
################################################################################################

################################################################################################
# Formatting for full meta-analysis
################################################################################################

# Females
load("summary_statistics/FOR2107/DMPs_BaselineResponse_Females_10CtrlPCs_2AncPCs_limma.RData")
FOR2107= data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/UTHealth/DMPs_BaselineResponse_Females_10CtrlPCs_5AncPCs_limma.RData")
UTHealth = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/UNICA_BD/DMPs_BaselineResponse_Females_10CtrlPCs_2AncPCs_limma.RData")
UNICA = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/HALIFAX-CAGLIARI/DMPs_BaselineResponse_Females_15CtrlPCs_2AncPCs_limma.RData")
HALICAG = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/BIPOGENT_IPM/DMPs_BaselineResponse_Females_10CtrlPCs_10AncPCs_limma.RData")
BIPOGENT = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/TOP/PhaseI_DMPs_BaselineResponse_Females_10CtrlPCs_2AncPCs_limma.RData")
TOP1 = data.frame(row.names = rownames(DMPs),logFC=DMPs$bacon_corr_EFFECT,SE=DMPs$bacon_corr_SE)
load("summary_statistics/TOP/PhaseIII_DMPs_BaselineResponse_Females_15CtrlPCs_2AncPCs_limma.RData")
TOP3 = data.frame(row.names = rownames(DMPs),logFC=DMPs$bacon_corr_EFFECT,SE=DMPs$bacon_corr_SE)
load("summary_statistics/TOP/PhaseIV_DMPs_BaselineResponse_Females_15CtrlPCs_5AncPCs_limma.RData")
TOP4 = data.frame(row.names = rownames(DMPs),logFC=DMPs$bacon_corr_EFFECT,SE=DMPs$bacon_corr_SE)
res = list(FOR2107, UTHealth, UNICA, HALICAG, BIPOGENT, TOP1, TOP3, TOP4)
names(res) = c("FOR2107_F", "UTHealth_F", "UNICA_F","HALICAG_F","BIPOGENT_F", "TOP1_F", "TOP3_F", "TOP4_F")
save(res, file="data_ready_for_meta/Female_BIP.RData")

# Males
load("summary_statistics/FOR2107/DMPs_BaselineResponse_Males_10CtrlPCs_2AncPCs_limma.RData")
FOR2107= data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/UTHealth/DMPs_BaselineResponse_Males_10CtrlPCs_5AncPCs_limma.RData")
UTHealth = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/UNICA_BD/DMPs_BaselineResponse_Males_5CtrlPCs_2AncPCs_limma.RData")
UNICA = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/HALIFAX-CAGLIARI/DMPs_BaselineResponse_Males_15CtrlPCs_2AncPCs_limma.RData")
HALICAG = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/BIPOGENT_IPM/DMPs_BaselineResponse_Males_10CtrlPCs_10AncPCs_limma.RData")
BIPOGENT = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/TOP/PhaseI_DMPs_BaselineResponse_Males_10CtrlPCs_2AncPCs_limma.RData")
TOP1 = data.frame(row.names = rownames(DMPs),logFC=DMPs$bacon_corr_EFFECT,SE=DMPs$bacon_corr_SE)
load("summary_statistics/TOP/PhaseIII_DMPs_BaselineResponse_Males_15CtrlPCs_2AncPCs_limma.RData")
TOP3 = data.frame(row.names = rownames(DMPs),logFC=DMPs$bacon_corr_EFFECT,SE=DMPs$bacon_corr_SE)
load("summary_statistics/TOP/PhaseIV_DMPs_BaselineResponse_Males_15CtrlPCs_5AncPCs_limma.RData")
TOP4 = data.frame(row.names = rownames(DMPs),logFC=DMPs$bacon_corr_EFFECT,SE=DMPs$bacon_corr_SE)
load("summary_statistics/PDMH/DMPs_BaselineResponse_Males_10CtrlPCs_10AncPCs_limma.RData")
PDMH = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
res = list(FOR2107, UTHealth, UNICA, HALICAG, BIPOGENT, TOP1, TOP3, TOP4, PDMH)
names(res) = c("FOR2107_M", "UTHealth_M", "UNICA_M","HALICAG_M","BIPOGENT_M", "TOP1_M", "TOP3_M", "TOP4_M","PDMH_M")
save(res, file="data_ready_for_meta/Male_BIP.RData")

# Sex-agnostic
load("summary_statistics/UNSW_BIP/transfer_3278297_files_a2bab723/DMPs_BaselineResponse_Sex_adjusted_5CtrlPCs_2AncPCs_limma.RData")
UNSW_BIP = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/IGP/OneDrive_3_12-2-2024/DMPs_BaselineResponse_Sex_adjusted_2CtrlPCs_2AncPCs_limma_Ficoll2.RData")
IGP_F = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/IGP/OneDrive_3_12-2-2024/DMPs_BaselineResponse_Sex_adjusted_15CtrlPCs_2AncPCs_limma_Blood2 1.RData")
IGP_B = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
load("summary_statistics/PDMH/DMPs_BaselineResponse_Sex_adjusted_10CtrlPCs_10AncPCs_limma.RData")
PDMH = data.frame(row.names = rownames(DMPs),logFC=DMPs$logFC,SE=DMPs$SE)
sex_agnostic = list(UNSW_BIP,IGP_F,IGP_B,PDMH)
names(sex_agnostic) = c("UNSW_BIP","IGP_Fi","IGP_B", "PDMH")

resF = get(load("data_ready_for_meta/Female_BIP.RData"));rm(res)
resM = get(load("data_ready_for_meta/Male_BIP.RData"));rm(res)
resM = resM[-grep("PDMH",names(resM))]
names(resM)
res = c(resM,resF,sex_agnostic);rm(resM,resF)

save(res, file="data_ready_for_meta/Sex_agnostic_BIP.RData")



################################################################################################
# Formatting for Sex agnostic meta-analysis leaving out TOP4
################################################################################################

load("data_ready_for_meta/Sex_agnostic_BIP.RData")
names(res)
res = res[-grep("TOP4",names(res))]
names(res)
save(res, file="data_ready_for_meta/Sex_agnostic_BIP_LeavingOutTOP4.RData")


