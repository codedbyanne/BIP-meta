#!/usr/bin/Rscript

#################################################################################################
# Purpose of script: Format EWAS results before running Comb-p
# @author: Anne-Kristin Stavrum
################################################################################################

print("loading libraries")
load("Probe_coordinates_Hg38_EPIC_V1_V2.Robj")

args =commandArgs(trailingOnly = TRUE)

dir = args[1]

files = list.files(path=dir, pattern = ".Robj")


# print("Converting files...")
for(i in 1:length(files)){
  dmps = get(load(file.path(dir,files[i])))
  length(intersect(rownames(anno),rownames(dmps)))
  dmps = dmps[rownames(dmps) %in% rownames(anno),]
  anno = anno[match(rownames(dmps),rownames(anno)),]
  #bed = data.frame(chrom = gsub("chr","",anno$chr),start=anno$pos,end=anno$pos,pval=dmps$All_P_Fixed)
  bed = data.frame(chrom = gsub("chr","",anno$chr),start=anno$pos,end=anno$pos,pval=dmps$All_P_Random)
  bed = bed[with(bed, order(chrom,start)),]
  outfile = paste(c(strsplit(files[i],"\\.R")[[1]][1],"_Random.bed"),collapse = "")
  print(outfile)
  write.table(bed,file=file.path(dir,outfile),sep="\t",row.names = FALSE,quote = FALSE)
  
}





