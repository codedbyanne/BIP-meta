library(missMethyl)

#############################################
# Loading data
#############################################

load("UPDATED_ANNOTATION_v1v2.RData")

load("meta_analysis/Meta_BIP_SexAgnostic.Robj")
res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)
res.meta = res.meta[with(res.meta, order(All_P_Random)),]
cpgs = rownames(res.meta)[res.meta$All_P_Random<1e-5]
cpgs = rownames(res.meta)[1:(length(rownames(res.meta))*0.005)]
length(intersect(cpgs,annov1v2$Name))


#############################################
# PW analysis
#############################################
#anno_diff = anno[which(anno$Name %in% setdiff(anno$Name,anno_red$Name)),]
pw = gometh(sig.cpg = cpgs, anno=annov1v2, sig.genes = T, collection = "GO", array.type="EPIC")
pw6.sorted = pw[with(pw,order(pw$P.DE)),]
save(pw6.sorted, file="pathway/GO_cpg_topHALFpercent_Encode_Refseq.Robj")

pw = gometh(sig.cpg = cpgs, anno=annov1v2, sig.genes = T, collection = "KEGG")
pw6.sorted = pw[with(pw,order(pw$P.DE)),]
save(pw6.sorted, file="pathway/KEGG_cpg_topHALFpercent_Encode_Refseq.Robj")

########
# MSigDB Collection 2 Formatting
msig = read.csv("pathway/c2.all.v2024.1.Hs.entrez.gmt",header = F)
msig = sapply(msig, function(x)strsplit(x,"\t"))
names(msig) = sapply(msig, function(x) x[1])
msig = sapply(msig, function(x) x[3:length(x)])

pw = gsameth(sig.cpg = cpgs, anno=annov1v2, collection = msig, sig.genes = T)
pw6_msig.sorted = pw[with(pw,order(pw$P.DE)),]
save(pw6_msig.sorted, file="pathway/MsigDB_C2_cpg_topHALFpercent_Encode_Refseq.Robj")