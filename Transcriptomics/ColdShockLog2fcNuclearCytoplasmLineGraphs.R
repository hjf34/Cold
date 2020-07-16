###Author: Harry Fischl
###Log2 fold change in cytoplasmic (black line, left axis) and nuclear (red line, right axis) RNA level of key genes at time points during the transfer of AC16 cells and U2OS cells from 37°C to 18°C for 24h and then back to 37°C for 24h relative to cells kept at 37°C. 
###Figure 3B and Figure EV3F-G
###############################Work directories

TWD = getwd()
TWD1 = "./DESeq/DifferentialGeneExpressionGeneLists/"
TWD2 = "./Rplots/"
setwd(TWD1)

library(scales)

lf1 = list.files(pattern = "AC16_18d")
lf1a = lf1[!grepl("90min", lf1)]
lf2 = as.vector(sapply(lf1a, function(n) strsplit(n, "vs_")[[1]][2]))
lf3 = as.vector(sapply(lf2, function(n) strsplit(n, "_DES")[[1]][1]))

ulf1 = list.files(pattern = "U2OS_18d")
ulf2 = as.vector(sapply(ulf1, function(n) strsplit(n, "vs_")[[1]][2]))
ulf3 = as.vector(sapply(ulf2, function(n) strsplit(n, "_DES")[[1]][1]))

##############################################
dslist1 = list()
for(a1 in 1:length(lf1a)){
  dslist1[[a1]] = read.csv(lf1a[a1], stringsAsFactors = F)
}

names(dslist1) = lf3
##############################################
dslist2 = list()
for(a1 in 1:length(ulf1)){
  dslist2[[a1]] = read.csv(ulf1[a1], stringsAsFactors = F)
}

names(dslist2) = ulf3
##############################################

cyt1 = dslist1[c(7,1,2,5,6,3,4)]
nuc1 = dslist1[c(14,8,9,12,13,10,11)]

cyt2 = dslist2[c(3,1,2)]
nuc2 = dslist2[c(6,4,5)]

lineplotter = function(genename, cyt, nuc){
  #######Collate DESEQ data
  cmat = as.data.frame(matrix(NA, nrow=length(cyt), ncol=7))
  nmat = as.data.frame(matrix(NA, nrow=length(nuc), ncol=7))
  for(a1 in 1:length(cyt)){
    dst = cyt[[a1]]
    cmat[a1,] = dst[dst$X == genename,]
    dst = nuc[[a1]]
    nmat[a1,] = dst[dst$X == genename,]
  }
  names(cmat) = names(dst)
  names(nmat) = names(dst)
  row.names(cmat) = names(cyt)
  row.names(nmat) = names(nuc)
  allmat = rbind(nmat,cmat)
  ###########################################
  ####PLOT CYT IN BLACK
  par(mar=c(1,3,1,3))
  m1 = c(0,cmat$log2FoldChange) ##l2fc
  s1 = c(0,cmat$lfcSE) ###Standard error
  ylimhigher = max(m1+s1)*1.1
  ylimlower = min(m1-s1)*1.05
  if(ylimhigher < 1){ylimhigher = 1}
  if(ylimlower > -1){ylimlower = -1}
  ############################################################################
  #rectangle dimensions
  r1 = grep("_18d",row.names(cmat))
  r2 = grep("h37d",row.names(cmat))
  r3 = c(r1,r2)
  r4 = r3[!(duplicated(r3)|duplicated(r3, fromLast=TRUE))]
  cold = c(r4[1]+0.5,ylimlower*10,r4[length(r4)]+1.5,ylimhigher*10)
  warma = c(0,ylimlower*10,r4[1]+0.5,ylimhigher*10)
  warmb = c(r4[length(r4)]+1.5,ylimlower*10,r2[length(r2)]+2,ylimhigher*10)
  ############################################################################
  xlim1 = 1
  xlim2 = length(m1)
  if(xlim2 < 5){xlim1 = 0.7; xlim2 = xlim2+0.3}
  plot(m1, type="o", pch=20, las=1, ylim=c(ylimlower,ylimhigher), xlim = c(xlim1,xlim2),
       main = "", cex.axis = 1, xlab="",ylab="", xaxt="n", col="white")
  
  rect(warma[1],warma[2],warma[3],warma[4],col=alpha("yellow",0.2), border=NA)
  rect(cold[1],cold[2],cold[3],cold[4],col=alpha("skyblue",0.2), border=NA)
  rect(warmb[1],warmb[2],warmb[3],warmb[4],col=alpha("yellow",0.2), border=NA)
  
  points(m1, type="o", pch=20, las=1, ylim=c(ylimlower,ylimhigher), 
         main = "", cex.axis = 1, xlab="",ylab="", xaxt="n")
  mp = 1:length(m1)
  subseg1 = m1-s1
  subseg2 = m1+s1
  segments(mp, subseg1, mp, subseg2, lwd=1)
  segments(mp-0.2, subseg1, mp+0.2, subseg1, lwd=1)
  segments(mp-0.2, subseg2, mp+0.2, subseg2, lwd=1)
  axis(1, labels=rep("",length(m1)), at = mp)
  ########################################
  par(new=T)
  ###Plot Nuclear
  par(mar=c(1,3,1,3))
  m1 = c(0,nmat$log2FoldChange)
  s1 = c(0,nmat$lfcSE)
  ylimhigher = max(m1+s1)*1.1
  ylimlower = min(m1-s1)*1.05
  if(ylimhigher < 1){ylimhigher = 1}
  if(ylimlower > -1){ylimlower = -1}
  
  plot(m1, type="o", pch=1, las=1, ylim=c(ylimlower,ylimhigher), xlim = c(xlim1,xlim2),
       main = "", xlab="",ylab="", xaxt="n",yaxt="n", col="red")
  mp = 1:length(m1)
  subseg1 = m1-s1
  subseg2 = m1+s1
  segments(mp, subseg1, mp, subseg2, lwd=1, col="red")
  segments(mp-0.1, subseg1, mp+0.1, subseg1, lwd=1, col="red")
  segments(mp-0.1, subseg2, mp+0.1, subseg2, lwd=1, col="red")
  z = axis(4, col="red", cex.axis = 1, labels=F)
  mtext(side = 4, text=z,at=z, col="red", las=1, line=1, cex = 1)
  
  box(lwd=2)
  return(allmat)
}

ac1 = lineplotter("NR1D1", cyt1, nuc1)
ac2 = lineplotter("CRY1", cyt1, nuc1)
ac3 = lineplotter("CRY2", cyt1, nuc1)
ac4 = lineplotter("PER1", cyt1, nuc1)
ac5 = lineplotter("PER2", cyt1, nuc1)
ac6 = lineplotter("ARNTL", cyt1, nuc1)
ac7 = lineplotter("CLOCK", cyt1, nuc1)
ac8 = lineplotter("TP53", cyt1, nuc1)

uo1 = lineplotter("NR1D1", cyt2, nuc2)
uo2 = lineplotter("CRY1", cyt2, nuc2)
uo3 = lineplotter("CRY2", cyt2, nuc2)
uo4 = lineplotter("PER1", cyt2, nuc2)
uo5 = lineplotter("PER2", cyt2, nuc2)
uo6 = lineplotter("ARNTL", cyt2, nuc2)
uo7 = lineplotter("CLOCK", cyt2, nuc2)

u2osdata = rbind(uo1,uo2,uo3,uo4,uo5,uo6,uo7)
ac16data = rbind(ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8)

setwd(TWD)
setwd(TWD2)
write.csv(ac16data, "l2fcSE_DESeqData_AC16_for_LineGraphs.csv")
write.csv(u2osdata, "l2fcSE_DESeqData_U2OS_for_LineGraphs.csv")

height1 = 2.5
width1 = 3.5

png("l2fcSE_CRY1_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("CRY1", cyt1, nuc1)
dev.off()

png("l2fcSE_CRY2_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("CRY2", cyt1, nuc1)
dev.off()

png("l2fcSE_PER1_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("PER1", cyt1, nuc1)
dev.off()

png("l2fcSE_PER2_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("PER2", cyt1, nuc1)
dev.off()

png("l2fcSE_ARNTL_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("ARNTL", cyt1, nuc1)
dev.off()

png("l2fcSE_CLOCK_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("CLOCK", cyt1, nuc1)
dev.off()

png("l2fcSE_TP53_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("TP53", cyt1, nuc1)
dev.off()

png("l2fcSE_NR1D1_AC16NucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("NR1D1", cyt1, nuc1)
dev.off()

###################################################
##########################################################################

height1 = 2.5
width1 = 2.5

png("l2fcSE_CRY1_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("CRY1", cyt2, nuc2)
dev.off()

png("l2fcSE_CRY2_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("CRY2", cyt2, nuc2)
dev.off()

png("l2fcSE_PER1_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("PER1", cyt2, nuc2)
dev.off()

png("l2fcSE_PER2_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("PER2", cyt2, nuc2)
dev.off()

png("l2fcSE_ARNTL_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("ARNTL", cyt2, nuc2)
dev.off()

png("l2fcSE_CLOCK_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("CLOCK", cyt2, nuc2)
dev.off()

png("l2fcSE_TP53_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("TP53", cyt2, nuc2)
dev.off()

png("l2fcSE_NR1D1_U2OSNucAndCyt.png", height = height1, width = width1, units = "in", res = 600)
lineplotter("NR1D1", cyt2, nuc2)
dev.off()

###################################################
