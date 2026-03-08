################################################################################################
# This script is adapted from a script available from https://pubmed.ncbi.nlm.nih.gov/36908043/
# https://github.com/jche453/Pruning-Thresholding-MRS.git
################################################################################################

load("./Init_Data_For_Comeback_EPICv1v2_Hg38.RData")
  
cmr = function (Mdata, meds = NULL, Iarray = c("450K", "EPIC","EPICv1v2"), cormethod = c("pearson", 
                                                                        "kendall", "spearman"), corlo = NULL, corhi = corlo, corlodst = 400, 
          corhidst = 1, maxprbdst = 2000, maxdlvl = Inf, verbose = TRUE) 
{
  cormethod <- match.arg(cormethod)
  Iarray <- match.arg(Iarray)
  if (is.null(corlo) && is.null(corhi)) {
    corlo = corhi = ccutN_f(nrow(Mdata))
    print(paste0("No correlation cut-off specifified, using ad-hoc sample-size based cut-off ", 
                 corlo))
  }
  if (verbose) 
    print("Getting CpG info, about to start estimation any day now.")
  
  if (Iarray == "450K"){
    EPIC_Manifest = comeback:::init_data$I450K_Manifest
    chr_seq_GpCpos = comeback:::init_data$chr_seq_GpCpos
  } 
    
  else if(Iarray == "EPIC"){
    EPIC_Manifest = comeback:::init_data$EPIC_Manifest
    chr_seq_GpCpos = comeback:::init_data$chr_seq_GpCpos
  }
  else if(Iarray == "EPICv1v2"){
    EPIC_Manifest = init_data_hg38$EPICv1v2
    chr_seq_GpCpos = init_data_hg38$chr_seq_GpCpos
  } 
  
  EPIC_OK_prb = rownames(EPIC_Manifest)
  prb_MMat = colnames(Mdata)
  wkprbs = intersect(EPIC_OK_prb, prb_MMat)
  numprbs = length(wkprbs)
  EPIC_Manifest = EPIC_Manifest[wkprbs, ]
  EPIC_Manifest$Name = wkprbs
  EPIC_Manifest = EPIC_Manifest[order(EPIC_Manifest$CHR, EPIC_Manifest$MAPINFO),]
  wkprbs_crd = EPIC_Manifest$MAPINFO
  names(wkprbs_crd) = rownames(EPIC_Manifest)
  chrom_nams = unique(EPIC_Manifest$CHR)
  chrom_nams = chrom_nams[order(chrom_nams)]
  chroms = length(chrom_nams)
  cmr_ac = vector(mode = "list", length = chroms)
  names(cmr_ac) = paste0("chr", chrom_nams)
  if (verbose) 
    print(paste("Found ", numprbs, " probes from ", chroms, 
                " out of 24 chromosomes", collapse = ""))
  if ((maxdlvl == 1) | (is.null(meds))) {
    for (i in 1:chroms) {
      cni = chrom_nams[[i]]
      manchr = EPIC_Manifest[EPIC_Manifest$CHR == cni, 
                             c("Name", "MAPINFO")]
      manchr = manchr[order(manchr$MAPINFO), ]
      x = manchr$Name
      x_prbcrd = manchr$MAPINFO
      cni = paste0("chr", chrom_nams[[i]])
      cpgi = match(x_prbcrd, chr_seq_GpCpos[[cni]])
      cpgi_tf = !is.na(cpgi)
      cpgi = cpgi[cpgi_tf]
      x = x[cpgi_tf]
      x_prbcrd = x_prbcrd[cpgi_tf]
      ccut_tf = sapply(2:length(cpgi), function(y) {
        res = FALSE
        ym1 = y - 1
        dcpg = chr_seq_GpCpos[[cni]][(cpgi[[ym1]] + 1):cpgi[[y]]] - 
          chr_seq_GpCpos[[cni]][cpgi[[ym1]]:(cpgi[[y]] - 
                                               1)]
        maxdcpg = max(dcpg)
        if ((maxdcpg <= corlodst) & (abs(wkprbs_crd[[x[[(ym1)]]]] - 
                                         wkprbs_crd[[x[[(y)]]]]) <= maxprbdst)) {
          corcut = corhi - (maxdcpg - corhidst) * (corhi - 
                                                     corlo)/(corlodst - corhidst)
          res = (cor(x = Mdata[, x[[(ym1)]]], y = Mdata[, 
                                                        x[[y]]], use = "pairwise.complete.obs", method = cormethod) > 
                   corcut)
        }
        return(res)
      })
      ccut_i = which(ccut_tf)
      lpc = length(ccut_i)
      cmr1dcr = NULL
      if (lpc > 0) {
        cc = c(ccut_i[[1]])
        cp = 1
        res = list()
        cct = 0
        if (lpc > 1) {
          while (cp < lpc) {
            cp1 = cp + 1
            if ((ccut_i[[cp1]] - ccut_i[[cp]]) == 1) {
              cc = c(cc, ccut_i[[cp1]])
            }
            else {
              cct = cct + 1
              res[[cct]] = cc
              cc = c(ccut_i[[cp1]])
              if (cp1 == lpc) 
                res[[cct + 1]] = cc
            }
            cp = cp + 1
          }
        }
        else {
          res[[1]] = cc
        }
        cmr1dcr = lapply(res, function(z) {
          lz = length(z)
          if (lz == 0) 
            print(res)
          if (lz == 1) {
            cl = x[z:(z + 1)]
          }
          else {
            cl = x[z[[1]]:(z[[lz]] + 1)]
          }
          return(cl)
        })
      }
      cmr_ac[[cni]] = cmr1dcr
      cmr_ac[[cni]] = cmr_ac[[cni]][sapply(cmr_ac[[cni]], 
                                           length) > 0]
      if (verbose) 
        print(paste("Done chr", as.character(cni), "with", 
                    length(cmr_ac[[i]]), "cmrs"))
    }
  }
  else {
    print(sprintf("Applying adjacent probe filter: max level difference %f", 
                  maxdlvl))
    names(meds) = colnames(Mdata)
    for (i in 1:chroms) {
      cni = chrom_nams[[i]]
      manchr = EPIC_Manifest[EPIC_Manifest$CHR == cni, 
                             c("Name", "MAPINFO")]
      manchr = manchr[order(manchr$MAPINFO), ]
      x = manchr$Name
      x_prbcrd = manchr$MAPINFO
      cni = paste0("chr", chrom_nams[[i]])
      cpgi = match(x_prbcrd, chr_seq_GpCpos[[cni]])
      cpgi_tf = !is.na(cpgi)
      cpgi = cpgi[cpgi_tf]
      x = x[cpgi_tf]
      x_prbcrd = x_prbcrd[cpgi_tf]
      ccut_tf = sapply(2:length(cpgi), function(y) {
        res = FALSE
        ym1 = y - 1
        dcpg = chr_seq_GpCpos[[cni]][(cpgi[[ym1]] + 1):cpgi[[y]]] - 
          chr_seq_GpCpos[[cni]][cpgi[[ym1]]:(cpgi[[y]] - 
                                               1)]
        maxdcpg = max(dcpg)
        if ((maxdcpg <= corlodst) & (abs(wkprbs_crd[[x[[(ym1)]]]] - 
                                         wkprbs_crd[[x[[(y)]]]]) <= maxprbdst)) {
          corcut = corhi - (maxdcpg - corhidst) * (corhi - 
                                                     corlo)/(corlodst - corhidst)
          res = ((cor(x = Mdata[, x[[(ym1)]]], y = Mdata[, 
                                                         x[[y]]], use = "pairwise.complete.obs") > 
                    corcut) & (abs(meds[[x[[(ym1)]]]] - meds[[x[[y]]]]) < 
                                 maxdlvl))
        }
        return(res)
      })
      ccut_i = which(ccut_tf)
      lpc = length(ccut_i)
      cmr1dcr = NULL
      if (lpc > 0) {
        cc = c(ccut_i[[1]])
        cp = 1
        res = list()
        cct = 0
        if (lpc > 1) {
          while (cp < lpc) {
            cp1 = cp + 1
            if ((ccut_i[[cp1]] - ccut_i[[cp]]) == 1) {
              cc = c(cc, ccut_i[[cp1]])
            }
            else {
              cct = cct + 1
              res[[cct]] = cc
              cc = c(ccut_i[[cp1]])
              if (cp1 == lpc) 
                res[[cct + 1]] = cc
            }
            cp = cp + 1
          }
        }
        else {
          res[[1]] = cc
        }
        cmr1dcr = lapply(res, function(z) {
          lz = length(z)
          if (lz == 0) 
            print(res)
          if (lz == 1) {
            cl = x[z:(z + 1)]
          }
          else {
            cl = x[z[[1]]:(z[[lz]] + 1)]
          }
          return(cl)
        })
      }
      cmr_ac[[cni]] = cmr1dcr
      cmr_ac[[cni]] = cmr_ac[[cni]][sapply(cmr_ac[[cni]], 
                                           length) > 0]
      if (verbose) 
        print(paste("Done chr", as.character(cni), "with", 
                    length(cmr_ac[[cni]]), "cmrs"))
    }
  }
  return(cmr_ac)
}
