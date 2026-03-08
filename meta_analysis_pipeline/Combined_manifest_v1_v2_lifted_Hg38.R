library(comeback)

########################################
# Array CpG positions
########################################
anno_lifted = get(load("analyses/mrs/Probe_coordinates_Hg38.Robj"));rm(anno)

# Checing format of the dataset 
head(comeback:::init_data$EPIC_Manifest)

manifest_v1v2_Hg38 = data.frame(row.names = rownames(anno_lifted),CHR = gsub("chr",replacement ="" ,anno_lifted$chr),MAPINFO=anno_lifted$pos)

# Checking format of the created dataset
head(manifest_v1v2_Hg38)

########################################
# Get genomic CpG positions
########################################
library(BSgenome)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

genome = BSgenome.Hsapiens.UCSC.hg38

find_cpg_sites <- function(sequence){
  cpg_positions <- as.numeric(gregexpr("CG", sequence)[[1]])
  return(cpg_positions)
}

cpg_sites_list = lapply(seqnames(genome)[1:24], function(chrom){
  seq <- getSeq(genome, chrom)
  cpg_positions <- find_cpg_sites(as.character(seq))
  GRanges(seqnames = chrom, ranges = IRanges(start=cpg_positions, width=2))
})

cpg_sites <- do.call(c, cpg_sites_list)

chr_seq_CpGpos_hg38 = lapply(cpg_sites_list, function(chrom){
  chrom@ranges@start
})
names(chr_seq_CpGpos_hg38)=seqnames(genome)[1:24]

########################################
# Saving array and genomic CpG positions to a list to be used by comback_mod.R
########################################
init_data_hg38=list(chr_seq_GpCpos=chr_seq_CpGpos_hg38,EPICv1v2=manifest_v1v2_Hg38)

save(init_data_hg38,file="analyses/mrs/Init_Data_For_Comeback_EPICv1v2_Hg38.RData")
