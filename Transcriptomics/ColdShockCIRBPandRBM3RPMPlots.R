###Author: Harry Fischl
###Mean nuclear and cytoplasmic RNA levels in reads per million (RPM) of known cold-induced genes CIRBP and RBM3 for AC16 cells kept at 37°C or transferred from 37°C to 28°C, 18°C, or 8°C for 24h. 
###Figure EV2B
###############################Work directories

TWD = getwd()
TWD1 = "./Rplots/"

load("./CountTableGeneration/CountTables/coldRefSeqGeneRPMCountTable.Rda")
rct1 = coldRefSeqGeneRPMCountTable

genes = c("CIRBP","RBM3")

rct2 = rct1[row.names(rct1) %in% genes,]

############################################################################
############################################################################
###########

conditions = c("AC16_37d_","AC16_28d24h_","AC16_18d24h_","AC16_8d24h_")
colours1 = c("gold","firebrick1","dodgerblue3","green3")

colSems = function(dataframe){
  semvec = rep(NA, length(dataframe))
  for(a1 in 1:length(dataframe)){
    vec = dataframe[,a1]
    vec1 = vec[!is.na(vec)]
    semvec[a1] = sd(vec1)/sqrt(length(vec1))
  }
  return(semvec)
}

boxplotter = function(mat1){
  par(mar=c(1,3,1,1))
  cm = colMeans(mat1, na.rm=T)
  cs = colSems(data.frame(mat1))
  par(lwd=2)
  
  par(lwd=2)
  mp = barplot(cm, ylim=c(0,1.1*max(mat1, na.rm=T)), col=colours1,
               ylab = "", las=1,cex.axis=1.5, border="black")
  par(lwd=1)
  
  subseg1 = cm-cs
  subseg2 = cm+cs
  segments(mp, subseg1, mp, subseg2, lwd=1)
  segments(mp-0.2, subseg1, mp+0.2, subseg1, lwd=1)
  segments(mp-0.2, subseg2, mp+0.2, subseg2, lwd=1)
  
  par(lwd=2)
  for(a1 in 1:length(mat1[1,])){
    vec = mat1[,a1]
    vec = sort(vec[!is.na(vec)])
    pos = rep(mp[a1],length(vec))
    pos = pos + c(-0.05,0.05)
    points(pos,vec, col=colours1[a1], pch=20)
    points(pos,vec, col="black", pch=21)
  }
  box(lwd=2)
}

boxprepper = function(genename){
  ccond = paste(conditions, "C", sep="")   
  ncond = paste(conditions, "N", sep="")   
  
  cmat = matrix(NA, nrow = 6, ncol = length(conditions))
  nmat = matrix(NA, nrow = 7, ncol = length(conditions))
  
  for(a1 in 1:length(conditions)){
    ###############################
    tc1 = ccond[a1]
    gene1 = as.numeric(rct2[genename,grep(tc1, colnames(rct2))])
    cmat[1:length(gene1),a1] = gene1
    ##############################
    tn1 = ncond[a1]
    gene1 = as.numeric(rct2[genename,grep(tn1, colnames(rct2))])
    nmat[1:length(gene1),a1] = gene1
  }
  
  #######################################################################
  cytname = "RPM_000_cyt.png"
  nucname = "RPM_000_nuc.png"
  cytname1 = sub("000", genename, cytname)
  nucname1 = sub("000", genename, nucname)
  
  png(cytname1, width = 3, height = 4, units = "in",res = 600)
  boxplotter(cmat)
  dev.off()
  
  png(nucname1, width = 3, height = 4, units = "in",res = 600)
  boxplotter(nmat)
  dev.off()
}

setwd(TWD1)
boxprepper("CIRBP")
boxprepper("RBM3")
