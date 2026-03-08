#################################################################################################
# Purpose of script: Extract and prepare phenotypes for residualsation of TOP4
# Most of this has been copied from the QC and preprocessing pipeline scripts
################################################################################################

#### LOAD DATA ####
betas = get(load(".../phaseIV/Bipolar_pipeline/output/RData/beta_final_autosomes.RData"))
pheno = as.data.frame(readxl::read_xlsx(".../phaseIV/Bipolar_pipeline/output/Tables/samplesheet_after_QC.xlsx"))
sample_names <- intersect(pheno$Basename, colnames(betas))

load(".../phaseIV/Bipolar_pipeline/output/RData/Positive_ctrlprobe_intensities.RData") #ctrl probes
load(".../phaseIV/Bipolar_pipeline/output/RData/Cellproportions.RData") 
load(".../phaseIV/Bipolar_pipeline/output/RData/comb_SNPs.RData")

#############################################

# VARIABLES
Age <- scale(as.numeric(pheno$Age))[, 1]
smoking <- scale(betas["cg05575921", ])[, 1]
Sex <- factor(pheno$Sex)

cellcounts_df <- as.data.frame(cellcounts)
# Cell proportion PC
cellcounts_df = cellcounts_df[rownames(cellcounts_df) %in% sample_names, ]
cell_PCs = as.data.frame(prcomp(cellcounts_df)$x)
colnames(cell_PCs) <- paste0("cell_", colnames(cell_PCs))
cell_PCs <- sapply(cell_PCs, scale)

# CTRL PROBE PCs
ctrl <- ctrl[, colnames(ctrl) %in% sample_names]
pca <- prcomp(na.omit(t(ctrl))) # run pca
ctrlprobe_PCAscores <- as.data.frame(pca$x) #extract PCA scores
colnames(ctrlprobe_PCAscores) <- paste0("Ctrl_", colnames(ctrlprobe_PCAscores))
ctrlprobe_PCAscores <- sapply(ctrlprobe_PCAscores, scale)

#### ANCESTRY PCs ###
# PCA
comb_SNPs_red <- comb_SNPs[, colnames(comb_SNPs) %in% sample_names]
pc <- prcomp(t(comb_SNPs_red))
anc_PCs <- as.data.frame(pc$x)
colnames(anc_PCs) <- paste0("anc_", colnames(anc_PCs))
anc_PCs <- sapply(anc_PCs, scale)

# MAKE DF WITH VARIABLES FOR MODEL
variables_df <- data.frame(Sex = Sex, Age = Age, smoking = smoking)
variables_df <- cbind(variables_df, cell_PCs[,1:5], ctrlprobe_PCAscores[,1:15], anc_PCs[,1:5])  

save(variables_df,file="analyses/mrs/Covariates_TOP4ewas_for_use_in_residualisation_of_TOP4.RData")
