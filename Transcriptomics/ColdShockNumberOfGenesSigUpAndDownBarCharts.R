###Author: Harry Fischl
###Bar charts of number of genes showing significant (padj < 0.05) differential transcript levels (up and down regulation) at each cold temperature condition v 37 degrees
###Figure 2A and Figure EV2F
###############################Work directories

TWD = getwd()
TWD1 = "./DESeq/DifferentialGeneExpressionGeneLists/"
TWD2 = "./Rplots/"
setwd(TWD1)
################################################################################

ctemp1 = read.csv("AC16_37d_C_vs_AC16_28d24h_C_DESEQ.csv", stringsAsFactors = F) #A
ctemp2 = read.csv("AC16_37d_C_vs_AC16_18d24h_C_DESEQ.csv", stringsAsFactors = F) #B
ctemp3 = read.csv("AC16_37d_C_vs_AC16_8d24h_C_DESEQ.csv", stringsAsFactors = F) #C
ntemp1 = read.csv("AC16_37d_N_vs_AC16_28d24h_N_DESEQ.csv", stringsAsFactors = F) #A
ntemp2 = read.csv("AC16_37d_N_vs_AC16_18d24h_N_DESEQ.csv", stringsAsFactors = F) #B
ntemp3 = read.csv("AC16_37d_N_vs_AC16_8d24h_N_DESEQ.csv", stringsAsFactors = F) #C

###24h at 28, 18 or 8 degrees C
ctemplist = list(ctemp1,ctemp2,ctemp3)
ntemplist = list(ntemp1,ntemp2,ntemp3)

ctime1 = read.csv("AC16_37d_C_vs_AC16_18d5h_C_DESEQ.csv", stringsAsFactors = F) #D
ctime2 = read.csv("AC16_37d_C_vs_AC16_18d10h_C_DESEQ.csv", stringsAsFactors = F) #E
ntime1 = read.csv("AC16_37d_N_vs_AC16_18d5h_N_DESEQ.csv", stringsAsFactors = F) #D
ntime2 = read.csv("AC16_37d_N_vs_AC16_18d10h_N_DESEQ.csv", stringsAsFactors = F) #E

###5, 10 , or 24h at 18 degrees C
ctimelist = list(ctime1,ctime2,ctemp2)
ntimelist = list(ntime1,ntime2,ntemp2)

#######################################################################################
setwd(TWD)
setwd(TWD2)

upregC = rep(NA, 3)
dnregC = rep(NA, 3)
upregN = rep(NA, 3)
dnregN = rep(NA, 3)
for(a1 in 1:length(ctemplist)){
  c1 = ctemplist[[a1]]
  c2 = c1[!is.na(c1$padj),]
  c3 = c2[c2$padj < 0.05,]
  upregC[a1] = length(c3[c3$log2FoldChange > 0,]$X)
  dnregC[a1] = length(c3[c3$log2FoldChange < 0,]$X)
  #########################################
  n1 = ntemplist[[a1]]
  n2 = n1[!is.na(n1$padj),]
  n3 = n2[n2$padj < 0.05,]
  upregN[a1] = length(n3[n3$log2FoldChange > 0,]$X)
  dnregN[a1] = length(n3[n3$log2FoldChange < 0,]$X)
}

png("DifferentTempsBarChartCytNucStacked_withNumbers.png", height = 3, width = 3, units = "in", res = 600)
layout(matrix(c(1,2), 2, byrow = TRUE))
q1 = c(upregN,upregC)
q2 = -c(dnregN,dnregC)
par(mar=c(0,4,1,1))
mp = barplot(q1,ylim=c(0,1600), las=1, col=c("firebrick1","dodgerblue3","green3"))
abline(v=mean(mp),lty=2,lwd=2, col="grey")
text(x = as.vector(mp), y=q1+100, labels = q1, cex=1)
box(lwd=2)
par(mar=c(1,4,0,1))
mp = barplot(q2,ylim=c(-1600,0), las=1, angle=45, density=35, col=c("firebrick1","dodgerblue3","green3"))
abline(v=mean(mp),lty=2,lwd=2, col="grey")
text(x = as.vector(mp), y=q2-100, labels = abs(q2), cex=1)
box(lwd=2)
dev.off()
#######################################################################################
#######################################################################################
#######################################################################################

upregC = rep(NA, 3)
dnregC = rep(NA, 3)
upregN = rep(NA, 3)
dnregN = rep(NA, 3)
for(a1 in 1:length(ctimelist)){
  c1 = ctimelist[[a1]]
  c2 = c1[!is.na(c1$padj),]
  c3 = c2[c2$padj < 0.05,]
  upregC[a1] = length(c3[c3$log2FoldChange > 0,]$X)
  dnregC[a1] = length(c3[c3$log2FoldChange < 0,]$X)
  #########################################
  n1 = ntimelist[[a1]]
  n2 = n1[!is.na(n1$padj),]
  n3 = n2[n2$padj < 0.05,]
  upregN[a1] = length(n3[n3$log2FoldChange > 0,]$X)
  dnregN[a1] = length(n3[n3$log2FoldChange < 0,]$X)
}

png("DifferentTimes18dBarChartCytNucStacked_withNumbers.png", height = 3, width = 3, units = "in", res = 600)
layout(matrix(c(1,2), 2, byrow = TRUE))
q1 = c(upregN,upregC)
q2 = -c(dnregN,dnregC)
par(mar=c(0,4,1,1))
mp = barplot(q1,ylim=c(0,1600), las=1, col=c("skyblue","brown","dodgerblue3"))
abline(v=mean(mp),lty=2,lwd=2, col="grey")
text(x = as.vector(mp), y=q1+100, labels = q1, cex=1)
box(lwd=2)
par(mar=c(1,4,0,1))
mp = barplot(q2,ylim=c(-1600,0), las=1, angle=45, density=35, col=c("skyblue","brown","dodgerblue3"))
abline(v=mean(mp),lty=2,lwd=2, col="grey")
text(x = as.vector(mp), y=q2-100, labels = abs(q2), cex=1)
box(lwd=2)
dev.off()

