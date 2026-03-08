library(liftOver)
library(GenomicRanges)
library(rtracklayer)

######
# Preparing epic v1 annotation for liftover with UCSC genome browser http://genome.ucsc.edu/cgi-bin/hgLiftOver
#####

manifest_epicv1.b5 = read.csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7, sep=",")
manifest_epicv1.b5$Name[which(is.na(manifest_epicv1.b5$MAPINFO))]
manifest_epicv1.b5 = manifest_epicv1.b5[which(!is.na(manifest_epicv1.b5$MAPINFO)),]
coordinates = paste(paste("chr",manifest_epicv1.b5$CHR,sep=""),paste(manifest_epicv1.b5$MAPINFO,manifest_epicv1.b5$MAPINFO,sep="-"),sep=":")
write.csv(coordinates, file = "Epic_V1_Coordinates.csv", quote = F, row.names = F)


########
# Combining liftover results with original coordinates in a dataframe
######
coordinates =read.csv("Epic_V1_Coordinates.csv")
epicv1.hg38 = read.csv("EpicV1_liftedToHg38.bed",header = F)
failed.lift = read.csv("Epic_v1_failedTheLiftOverToHg38.txt")
failed.lift = failed.lift[-grep("#",failed.lift$X.Deleted.in.new),]
coordinates = coordinates$x[-which(coordinates$x %in% intersect(coordinates$x,failed.lift))]
epicv1 = data.frame(hg19=coordinates)#,hg38=epicv1.hg38$V1)
epicv1$chr_hg19 = sapply(epicv1$hg19, function(x) strsplit(x,":")[[1]][1])
epicv1$pos_hg19 = sapply(epicv1$hg19, function(x) strsplit(strsplit(x,":")[[1]][2],"-")[[1]][1])
epicv1$hg38=epicv1.hg38$V1
epicv1$chr_hg38 = sapply(epicv1$hg38, function(x) strsplit(x,":")[[1]][1])
epicv1$pos_hg38 = sapply(epicv1$hg38, function(x) strsplit(strsplit(x,":")[[1]][2],"-")[[1]][1])

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19,quietly = TRUE)
anno_v1 = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_v1$coordinate = paste(anno_v1$chr,paste(anno_v1$pos,anno_v1$pos,sep="-"),sep=":")
length(intersect(epicv1$hg19,anno_v1$coordinate))
epicv1$Name = anno_v1$Name[match(epicv1$hg19,anno_v1$coordinate)]
rownames(epicv1) = epicv1$Name
any(is.na(epicv1$Name))
save(epicv1,file="Probe_coordinates_Hg38_EPICV1.Robj")
