
#############################################
# Format annotation
#############################################
anno1 = get(load("UPDATED_ANNOTATION.RData"))
anno1 = anno1[,c("Name","CHR","MAPINFO", "Gene","RegulRegion")]
anno2 = get(load("S1_Updated_EPICv2_annotation.RData"));rm(anno_upd)
anno2 = anno2[,c("Name","CHR","MAPINFO", "Gene","RegulRegion")]
anno2 = anno2[which(!anno2$CHR=="chr0"),]
anno1_diff = anno1[anno1$Name %in% setdiff(anno1$Name,anno2$Name),]
annov1v2 = rbind(anno2,anno1_diff)
colnames(annov1v2) = c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")
rownames(annov1v2) = annov1v2$Name
lengths = sapply(annov1v2$UCSC_RefGene_Name, function(x) length(strsplit(x,";")[[1]]))
lengths2 = sapply(annov1v2$UCSC_RefGene_Group, function(x) length(strsplit(x,";")[[1]]))
diffs = sapply(1:length(lengths), function(x) ifelse(lengths[x]==lengths2[x], T,F))

annov1v2$UCSC_RefGene_Group[which(!diffs)] = sapply(which(!diffs), function(x) {
  paste(annov1v2$UCSC_RefGene_Group[x],stringr::str_dup("; ",(lengths[x]-1)),sep="")
})
lengths3 = sapply(annov1v2$UCSC_RefGene_Group, function(x) length(strsplit(x,";")[[1]]))
diffs2 = sapply(1:length(lengths), function(x) ifelse(lengths[x]==lengths3[x], T,F))
length(which(!diffs2))

save(annov1v2, file="UPDATED_ANNOTATION_v1v2.RData")