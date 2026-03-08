#!/usr/bin/Rscript

print("loading libraries")
anno = get(load("UPDATED_ANNOTATION.RData"));rm(anno_upd)
colnames(anno)[grep("Gene",colnames(anno))] = "UCSC_RefGene_Name"

args =commandArgs(trailingOnly = TRUE)

#dir1 = args[1] # directory with .regions-p.bed files
dir1 ="DMRTEST/"

dmr.files = list.files(path=dir1,pattern = ".regions-p.bed")


################

for(f in 1:length(dmr.files)){
  print(dmr.files[f])

  dmrs = read.table(file.path(dir1,dmr.files[f]), header = F)
  colnames(dmrs) = c("chrom","start","end","min_p","n_probes","z_p","z_sidak_p")
  
  outfile = file.path(dir1,paste(c(strsplit(dmr.files[f],"\\.bed")[[1]][1],"_SignificantAndAnnotated_EncodeRefseq.Robj"),collapse = ""))
  outfile2 = file.path(dir1,paste(c(strsplit(dmr.files[f],"\\.bed")[[1]][1],"_SignificantAndAnnotated_EncodeRefseq.csv"),collapse = ""))
  
  dmrs.sig = dmrs[which(dmrs$z_sidak_p<0.05 & dmrs$n_probes>2),]
  
   if(length(dmrs.sig[,1])>0){
    dmr.probes=list()

    for(i in 1:length(dmrs.sig[,1])){
      dmr.probes[[i]] = c(anno$Name[which(gsub("chr","",anno$chr)==dmrs.sig$chrom[i] & anno$pos>= dmrs.sig$start[i] & anno$pos<= dmrs.sig$end[i])])
    }

    genes = list()
    for(i in 1:length(dmr.probes)){
      genes[[i]] = unique(anno$UCSC_RefGene_Name[which(anno$Name %in% dmr.probes[[i]])])
      genes[[i]] = unique(unlist(sapply(genes[[i]],function(x) strsplit(x,";"))))
    }
    genes.unlisted = unlist(genes)
    genes.unlisted

    genes.concatenated = list()
    for(i in 1:length(dmr.probes)){
      genes.concatenated[i] = paste(genes[[i]],collapse=";")
    }
    genes.concatenated = as.character(unlist(genes.concatenated))
    dmrs.sig$genes=genes.concatenated

    dmrs.sig = dmrs.sig[with(dmrs.sig, order(z_sidak_p)),]
    save(dmrs.sig,file=outfile)
    write.csv(dmrs.sig,file=outfile2, row.names = F)
  }
  
}

