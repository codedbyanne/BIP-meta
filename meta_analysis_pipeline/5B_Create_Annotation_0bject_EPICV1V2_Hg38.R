
########
# Create annotation object for v1 and v2
######
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19,quietly = TRUE)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38,quietly = TRUE)
anno_v1 = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_v2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_v2_aggregated = aggregate_to_probes(anno_v2);rm(anno_v2)
anno_v1_specific = anno_v1[setdiff(rownames(anno_v1),rownames(anno_v2_aggregated)),];rm(anno_v1)

anno_v1_specific$coordinate = paste(anno_v1_specific$chr,paste(anno_v1_specific$pos,anno_v1_specific$pos,sep="-"),sep=":")
anno_v2_aggregated$coordinate = paste(anno_v2_aggregated$chr,paste(anno_v2_aggregated$pos,anno_v2_aggregated$pos,sep="-"),sep=":")

########
# Adding Hg38 coordinates to EPICv1 specific probes
######
load("Probe_coordinates_Hg38_EPICV1.Robj")
length(intersect(epicv1$hg19,anno_v1_specific$coordinate))
missing_anno_v1_specific = anno_v1_specific[which(anno_v1_specific$coordinate %in% setdiff(anno_v1_specific$coordinate,epicv1$hg19)),] # failed lift
anno_v1_specific = anno_v1_specific[anno_v1_specific$coordinate %in% epicv1$hg19,]
anno_v1_specific$chr_hg38 = epicv1$chr_hg38[match(anno_v1_specific$coordinate, epicv1$hg19)]
anno_v1_specific$pos_hg38 = epicv1$pos_hg38[match(anno_v1_specific$coordinate, epicv1$hg19)]

########
# Creating a new dataframe with Hg38 coordinates for all CpGs from EPICv1 and EPICv2
######
anno = data.frame(row.names = union(rownames(anno_v1_specific),rownames(anno_v2_aggregated)))
anno_v2_aggregated = anno_v2_aggregated[match(rownames(anno),rownames(anno_v2_aggregated)),]
anno_v1_specific = anno_v1_specific[match(rownames(anno),rownames(anno_v1_specific)),]
anno$annotation[which(rownames(anno) %in% rownames(anno_v2_aggregated))] = "EPIC_V2"
anno$annotation[which(rownames(anno) %in% rownames(anno_v1_specific))] = "EPIC_V1"
anno$chr[which(rownames(anno) %in% rownames(anno_v2_aggregated))] = anno_v2_aggregated$chr[which(rownames(anno) %in% rownames(anno_v2_aggregated))]
anno$pos[which(rownames(anno) %in% rownames(anno_v2_aggregated))] = anno_v2_aggregated$pos[which(rownames(anno) %in% rownames(anno_v2_aggregated))]
anno$chr[which(rownames(anno) %in% rownames(anno_v1_specific))] = anno_v1_specific$chr_hg38[which(rownames(anno) %in% rownames(anno_v1_specific))]
anno$pos[which(rownames(anno) %in% rownames(anno_v1_specific))] = anno_v1_specific$pos_hg38[which(rownames(anno) %in% rownames(anno_v1_specific))]
anno$pos = as.integer(anno$pos)
save(anno, file="Probe_coordinates_Hg38_EPIC_V1_V2.Robj")
