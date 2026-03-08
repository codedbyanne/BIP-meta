library(ggplot2)
library(reshape2)
library(pheatmap)
library(gridExtra)
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# anno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
load("UPDATED_ANNOTATION_v1v2.RData")
anno = annov1v2;rm(annov1v2)
##############################################################################################################################
# Functions
##############################################################################################################################
jaccard_distance <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(union(a, b))
  return (1 - intersection / union)
}
clusterGL <- function(genelists){
  distances = matrix(nrow=length(genelists),ncol=length(genelists))
  for(i in 1:length(genelists)){
    for(j in 1:length(genelists)){
      #distances[i,j] = proxy::dist(genelists[[i]],genelists[[j]], method="Jaccard")
      distances[i,j] = jaccard_distance(genelists[[i]], genelists[[j]])
    }
    
  }
  
  dist_obj <- as.dist(distances, diag = TRUE)
  hclust_result <- hclust(dist_obj, method = "complete")
  plot(hclust_result)
  return(hclust_result)
}
#getClusters <- function(hclust_result, k, cl.size = 5){
getClusters <- function(hclust_result, k){
  cl = cutree(hclust_result, k=k)
  cluster_members <- split(seq_along(cl), cl)
  cluster_size = sapply(cluster_members, function(x) length(x))

  #cluster_members_n = cluster_members[which(cluster_size>cl.size)]

  return(cluster_members)
}
findCommonMembers <- function(genelists, cluster_members){
  
  clusterIntersect = list()
  for(i in 1:length(cluster_members)){
    clusterIntersect[[i]] = Reduce(intersect, genelists[cluster_members[[i]]])
  }
  return(clusterIntersect)
  
}
getClusterGenes <- function(genelists, cluster_members){
  
  allMembers = list()
  for(i in 1:length(cluster_members)){
    allMembers[[i]] = Reduce(union, genelists[cluster_members[[i]]])
  }
  return(allMembers)
  
}

createHeatMap <- function(limma_res, genelists, cluster_members_i, cl.genes_i){  
  dmps.cl = as.data.frame(t(sapply(cl.genes_i,function(x) limma_res[grep(x,limma_res$gene),][1,] ) ))
  if(any(is.na(dmps.cl$logFC))){
    remind = which(is.na(dmps.cl$logFC))
    dmps.cl = dmps.cl[-remind,]
    cl.genes_i = cl.genes_i[-remind]
  }
  
  data = matrix(0,nrow = length(names(genelists[cluster_members_i])), ncol = length(dmps.cl$logFC))
  rownames(data) = names(genelists[cluster_members_i]);colnames(data) = rownames(dmps.cl);
  
  for(i in 1:length(data[,1])){
    for(j in 1:length(data[1,])){
      g = cl.genes_i[j]
      if(g %in% genelists[cluster_members_i][[i]]){
        data[i,j] = as.numeric(dmps.cl$logFC[match(g,rownames(dmps.cl))][[1]])
      }
      
    }
  }
   clrsp <- colorRampPalette(c("blue", "white", "red"))   
   clrs <- clrsp(200) 
   
   breaks1 <- seq(-0.01, 0.01, length.out = 200)
  
  return(pheatmap(data, angle_col = 45, breaks = breaks1, color = clrs))
  #return(pheatmap(data, angle_col = 45))
}

createHeatMap_meta <- function(meta_res, genelists, cluster_members_i, cl.genes_i){  
  #dmps001 = limma_res[limma_res$P.Value<pval,]
  #dmps001.cl = as.data.frame(t(sapply(cl.genes_i,function(x) dmps001[grep(x,dmps001$gene),][1,] ) ))
  meta_res$gene = anno$UCSC_RefGene_Name[match(rownames(meta_res),anno$Name)]
  dmps.cl = as.data.frame(t(sapply(cl.genes_i,function(x) meta_res[grep(x,meta_res$gene),][1,] ) ))
  if(any(is.na(dmps.cl$All_Effect_Fixed))){
    remind = which(is.na(dmps.cl$All_Effect_Fixed))
    dmps.cl = dmps.cl[-remind,]
    cl.genes_i = cl.genes_i[-remind]
  }
  
  data = matrix(0,nrow = length(names(genelists[cluster_members_i])), ncol = length(dmps.cl$All_Effect_Fixed))
  rownames(data) = names(genelists[cluster_members_i]);colnames(data) = rownames(dmps.cl);
  
  for(i in 1:length(data[,1])){
    for(j in 1:length(data[1,])){
      g = cl.genes_i[j]
      if(g %in% genelists[cluster_members_i][[i]]){
        data[i,j] = as.numeric(dmps.cl$All_Effect_Fixed[match(g,rownames(dmps.cl))][[1]])
      }
      
    }
  }
  clrsp <- colorRampPalette(c("blue", "white", "red"))   
  clrs <- clrsp(200) 
  
  breaks1 <- seq(-0.01, 0.01, length.out = 200)
  
  return(pheatmap(data, angle_col = 45, breaks = breaks1, color = clrs))

}

  
##############################################################################################################################
# GSA summary Sex Agnostic Meta - GO
##############################################################################################################################

load("meta_analysis/Meta_BIP_SexAgnostic.Robj")

gsa = get(load("pathway/GO_cpg_topHALFpercent_Encode_Refseq.Robj"))
gsa = gsa[which(gsa$N>15 & gsa$N<500),]
gsa.sig = gsa[gsa$FDR<0.05,]

genelists = sapply(gsa.sig$SigGenesInSet,function(x) strsplit(x,","))
names(genelists) = gsa.sig$TERM
hclust_result = clusterGL(genelists=genelists)
cluster_members = getClusters(hclust_result, k=4)
overlap = findCommonMembers(genelists,cluster_members)
cl.genes = getClusterGenes(genelists,cluster_members)
cl.genes.length =sapply(cl.genes,length)
rem_cl = which(cl.genes.length<2)
if(length(rem_cl)>0){
  cl.genes  = cl.genes[-rem_cl]
  cluster_members = cluster_members[-rem_cl]
}
length(cluster_members)

# for(i in 1:length(cluster_members)){
#   createHeatMap_meta(res.meta.sorted, genelists,cluster_members[[i]],cl.genes[[i]])
# }

plot_list = list()
for(i in 1:length(cluster_members)){
  plot_list[[i]] = createHeatMap_meta(res.meta, genelists,cluster_members[[i]],cl.genes[[i]])$gtable
}
do.call(grid.arrange,plot_list)


##############################################################################################################################
# GSA summary Sex Agnostic Meta - MSIGDB
##############################################################################################################################

load("meta_analysis/Meta_BIP_SexAgnostic.Robj")

gsa = get(load("pathway/MsigDB_C2_cpg_topHALFpercent_Encode_Refseq.Robj"))
gsa = gsa[which(gsa$N>15 & gsa$N<500),]
#gsa$FDR = p.adjust(gsa$P.DE, method = "fdr")
gsa.sig = gsa[gsa$FDR<0.05,]

genelists = sapply(gsa.sig$SigGenesInSet,function(x) strsplit(x,","))
#names(genelists) = strtrim(gsa.sig$TERM,40)
names(genelists) = rownames(gsa.sig)
hclust_result = clusterGL(genelists=genelists)
cluster_members = getClusters(hclust_result, k=21)
overlap = findCommonMembers(genelists,cluster_members)
cl.genes = getClusterGenes(genelists,cluster_members)
cl.genes.length =sapply(cl.genes,length)
rem_cl = which(cl.genes.length<2)
if(length(rem_cl)>0){
  cl.genes  = cl.genes[-rem_cl]
  cluster_members = cluster_members[-rem_cl]
}
length(cluster_members)

# for(i in 1:length(cluster_members)){
#   createHeatMap_meta(res.meta.sorted, genelists,cluster_members[[i]],cl.genes[[i]])
# }

plot_list = list()
for(i in 1:length(cluster_members)){
  plot_list[[i]] = createHeatMap_meta(res.meta, genelists,cluster_members[[i]],cl.genes[[i]])$gtable
}
do.call(grid.arrange,plot_list)


##############################################################################################################################
# GSA summary Sex Agnostic Meta - KEGG
##############################################################################################################################

load("meta_analysis/Meta_BIP_SexAgnostic.Robj")

gsa = get(load("pathway/KEGG_cpg_topHALFpercent_Encode_Refseq.Robj"))
gsa = gsa[which(gsa$N>15 & gsa$N<500),]
#gsa$FDR = p.adjust(gsa$P.DE, method = "fdr")
gsa.sig = gsa[gsa$FDR<0.05,]

genelists = sapply(gsa.sig$SigGenesInSet,function(x) strsplit(x,","))
#names(genelists) = strtrim(gsa.sig$TERM,40)
names(genelists) = gsa.sig$Description
hclust_result = clusterGL(genelists=genelists)
cluster_members = getClusters(hclust_result, k=4)
overlap = findCommonMembers(genelists,cluster_members)
cl.genes = getClusterGenes(genelists,cluster_members)
cl.genes.length =sapply(cl.genes,length)
rem_cl = which(cl.genes.length<2)
if(length(rem_cl)>0){
  cl.genes  = cl.genes[-rem_cl]
  cluster_members = cluster_members[-rem_cl]
}
length(cluster_members)

# for(i in 1:length(cluster_members)){
#   createHeatMap_meta(res.meta.sorted, genelists,cluster_members[[i]],cl.genes[[i]])
# }

plot_list = list()
for(i in 1:length(cluster_members)){
  plot_list[[i]] = createHeatMap_meta(res.meta, genelists,cluster_members[[i]],cl.genes[[i]])$gtable
}
do.call(grid.arrange,plot_list)

