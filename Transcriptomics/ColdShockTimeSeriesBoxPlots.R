###Author: Harry Fischl
###Boxplots of standardized nuclear and cytoplasmic RNA levels at time points during the transfer of AC16 or U2OS cells from 37°C to 18°C or 28°C for 24h and then back to 37°C for 24h or 2h for either the group of genes showing significant up or genes showing significant downregulation in nuclear RNA levels at the 18°C or 28°C for 24h time point.
###Figure 3A and Figure EV3A-E
###############################Work directories

#####Boxplot Time Series

TWD = getwd()
TWD1 = "./DESeq/DifferentialGeneExpressionGeneLists/"
TWD2 = "./Rplots/"

library(DESeq2)

setwd(TWD)
load("./CountTableGeneration/CountTables/coldRefSeqGeneCountTable.Rda")

################################################################################
################################################################################
################################################################################
#########Regularized Log Transformation Of the Data

rlogtransformer = function(xmat1){
  #Creating ddsMat object
  c1 = colnames(xmat1)
  c2 = as.vector(sapply(c1, function(n) strsplit(n, "R")[[1]][1]))
  conditions = c2
  coldata <- data.frame(conditions)
  row.names(coldata) = c1
  ########################################################################
  ddsMat <- DESeqDataSetFromMatrix(countData = xmat1,
                                   colData = coldata,
                                   design = ~ conditions)
  rld <- assay(rlogTransformation(ddsMat, blind=F))
}

#########Mean of Regularized Log Transformed Data Repeats

repeataverager = function(dfxrld){
  nam = names(dfxrld)
  nam1 = as.vector(sapply(nam, function(x) strsplit(x, "R")[[1]][1]))
  n1 = unique(nam1)
  adf = data.frame(matrix(NA,nrow = dim(dfxrld)[1], ncol = length(n1)))
  names(adf) = n1
  row.names(adf) = row.names(dfxrld)
  for(a1 in 1:length(n1)){
    nx = n1[a1]
    cols1 = dfxrld[,grep(nx, names(dfxrld))]
    adf[,a1] = rowMeans(cols1)
  }
  return(adf)
}

#########Standardised across time series

timeseriesstandardiser = function(avdf){
  avmat = as.matrix(avdf)
  return((avmat - rowMeans(avmat))/apply(avmat,1, function(n) sd(n)))
}

################################################################################
################################################################################
################################################################################
################################################################################
#########AC16 18d time series

mct1 = coldRefSeqGeneCountTable
cn1 = colnames(mct1)
cn2 = cn1[c(grep("AC16_37d",cn1),grep("AC16_18d", cn1))]
cn2 = cn2[!grepl("90min", cn2)]

mat1 = mct1[,cn2]
cmat1 = mat1[,grep("CR", colnames(mat1))]
nmat1 = mat1[,grep("NR", colnames(mat1))]

crld = rlogtransformer(cmat1)
nrld = rlogtransformer(nmat1)

################################################################################
##################################################################################
######Averaging of repeats

dfcrld = as.data.frame(crld)
dfnrld = as.data.frame(nrld)

avcyt = repeataverager(dfcrld)
avnuc = repeataverager(dfnrld)

tscyt = timeseriesstandardiser(avcyt)
tsnuc = timeseriesstandardiser(avnuc)

##############################################################################
##################################################################################
######Get output from DESEQ
setwd(TWD)
setwd(TWD1)

nfx = read.csv("AC16_37d_N_vs_AC16_18d24h_N_DESEQ.csv", stringsAsFactors = F)
nfx1 = na.omit(nfx)

##########################################################################
##########################################################################
####Genes Up/Dn in Nucleus 18d24h

nfx1up = as.vector(nfx1[nfx1$padj < 0.05 & nfx1$log2FoldChange > 0,]$X)
nfx1dn = as.vector(nfx1[nfx1$padj < 0.05 & nfx1$log2FoldChange < 0,]$X)

cytup = tscyt[row.names(tscyt) %in% nfx1up,]
nucup = tsnuc[row.names(tsnuc) %in% nfx1up,]
cytdn = tscyt[row.names(tscyt) %in% nfx1dn,]
nucdn = tsnuc[row.names(tsnuc) %in% nfx1dn,]

#####################################################################
setwd(TWD)
setwd(TWD2)

colour1 = c("coral1","dodgerblue","dodgerblue","dodgerblue","coral1","coral1","coral1","coral1")
ylab1 = "Standardised RNA levels"

png("TimeSeriesBoxPlot_CytoBehaviour_UpNuc18d24h.png", height = 5, width = 5, units = "in", res = 600)
par(mar=c(1,4,3,1))
mp = boxplot(cytup, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Cytoplasm", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=2.5, cex=1.5)
box(lwd=2)
dev.off()

png("TimeSeriesBoxPlot_NucBehaviour_UpNuc18d24h.png", height = 5, width = 5, units = "in", res = 600)
par(mar=c(1,4,3,1))
mp = boxplot(nucup, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Nucleus", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=2.5, cex=1.5)
box(lwd=2)
dev.off()

#####################################################################

png("TimeSeriesBoxPlot_CytoBehaviour_DnNuc18d24h.png", height = 5, width = 5, units = "in", res = 600)
par(mar=c(1,4,3,1))
mp = boxplot(cytdn, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Cytoplasm", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=2.5, cex=1.5)
box(lwd=2)
dev.off()

png("TimeSeriesBoxPlot_NucBehaviour_DnNuc18d24h.png", height = 5, width = 5, units = "in", res = 600)
par(mar=c(1,4,3,1))
mp = boxplot(nucdn, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Nucleus", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=2.5, cex=1.5)
box(lwd=2)
dev.off()

#####################################################################
#####################################################################
#########paired t-test

cytup1 = cytup[,c(1,4,5)]
cytdn1 = cytdn[,c(1,4,5)]
nucup1 = nucup[,c(1,4,5)]
nucdn1 = nucdn[,c(1,4,5)]

#####################################################################

cut1 = unlist(t.test(cytup1[,1], cytup1[,2], paired = T))
cut2 = unlist(t.test(cytup1[,2], cytup1[,3], paired = T))
cut3 = unlist(t.test(cytup1[,1], cytup1[,3], paired = T))

cdt1 = unlist(t.test(cytdn1[,1], cytdn1[,2], paired = T))
cdt2 = unlist(t.test(cytdn1[,2], cytdn1[,3], paired = T))
cdt3 = unlist(t.test(cytdn1[,1], cytdn1[,3], paired = T))

nut1 = unlist(t.test(nucup1[,1], nucup1[,2], paired = T))
nut2 = unlist(t.test(nucup1[,2], nucup1[,3], paired = T))
nut3 = unlist(t.test(nucup1[,1], nucup1[,3], paired = T))

ndt1 = unlist(t.test(nucdn1[,1], nucdn1[,2], paired = T))
ndt2 = unlist(t.test(nucdn1[,2], nucdn1[,3], paired = T))
ndt3 = unlist(t.test(nucdn1[,1], nucdn1[,3], paired = T))

#####################################################################

AC16_T18d24h_pairedttest = rbind(cut1,cut2,cut3,cdt1,cdt2,cdt3,nut1,nut2,nut3,ndt1,ndt2,ndt3)
write.csv(AC16_T18d24h_pairedttest, file="TimeSeries_AC16_T18d24h_pairedttest_data.csv")

################################################################################
################################################################################
################################################################################
################################################################################
#########AC16 28d time series

mct1 = coldRefSeqGeneCountTable
cn1 = colnames(mct1)
cn2 = cn1[c(grep("AC16_37d",cn1),grep("AC16_28d", cn1))]

mat1 = mct1[,cn2]
cmat1 = mat1[,grep("CR", colnames(mat1))]
nmat1 = mat1[,grep("NR", colnames(mat1))]

crld = rlogtransformer(cmat1)
nrld = rlogtransformer(nmat1)

################################################################################
##################################################################################
######Averaging of repeats

dfcrld = as.data.frame(crld)
dfnrld = as.data.frame(nrld)

avcyt = repeataverager(dfcrld)
avnuc = repeataverager(dfnrld)

tscyt = timeseriesstandardiser(avcyt)
tsnuc = timeseriesstandardiser(avnuc)

##############################################################################
##################################################################################
######Get output from DESEQ

setwd(TWD)
setwd(TWD1)

nfx = read.csv("AC16_37d_N_vs_AC16_28d24h_N_DESEQ.csv", stringsAsFactors = F)
nfx1 = na.omit(nfx)

##########################################################################
##########################################################################
####Genes Up/Dn in Nucleus 28d24h

nfx1up = as.vector(nfx1[nfx1$padj < 0.05 & nfx1$log2FoldChange > 0,]$X)
nfx1dn = as.vector(nfx1[nfx1$padj < 0.05 & nfx1$log2FoldChange < 0,]$X)

cytup = tscyt[row.names(tscyt) %in% nfx1up,]
nucup = tsnuc[row.names(tsnuc) %in% nfx1up,]
cytdn = tscyt[row.names(tscyt) %in% nfx1dn,]
nucdn = tsnuc[row.names(tsnuc) %in% nfx1dn,]

#####################################################################

setwd(TWD)
setwd(TWD2)

colour1 = c("coral1","dodgerblue","coral1")
ylab1 = "Standardised RNA levels"

png("TimeSeriesBoxPlot_CytoBehaviour_UpNuc28d24h.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(cytup, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Cytoplasm", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()

png("TimeSeriesBoxPlot_NucBehaviour_UpNuc28d24h.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(nucup, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Nucleus", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()

#####################################################################
#####################################################################
#####################################################################
#####################################################################

png("TimeSeriesBoxPlot_CytoBehaviour_DnNuc28d24h.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(cytdn, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Cytoplasm", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()

png("TimeSeriesBoxPlot_NucBehaviour_DnNuc28d24h.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(nucdn, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Nucleus", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()


#####################################################################
#####################################################################
#########paired t-test

cytup1 = cytup
cytdn1 = cytdn
nucup1 = nucup
nucdn1 = nucdn

#####################################################################

cut1 = unlist(t.test(cytup1[,1], cytup1[,2], paired = T))
cut2 = unlist(t.test(cytup1[,2], cytup1[,3], paired = T))
cut3 = unlist(t.test(cytup1[,1], cytup1[,3], paired = T))

cdt1 = unlist(t.test(cytdn1[,1], cytdn1[,2], paired = T))
cdt2 = unlist(t.test(cytdn1[,2], cytdn1[,3], paired = T))
cdt3 = unlist(t.test(cytdn1[,1], cytdn1[,3], paired = T))

nut1 = unlist(t.test(nucup1[,1], nucup1[,2], paired = T))
nut2 = unlist(t.test(nucup1[,2], nucup1[,3], paired = T))
nut3 = unlist(t.test(nucup1[,1], nucup1[,3], paired = T))

ndt1 = unlist(t.test(nucdn1[,1], nucdn1[,2], paired = T))
ndt2 = unlist(t.test(nucdn1[,2], nucdn1[,3], paired = T))
ndt3 = unlist(t.test(nucdn1[,1], nucdn1[,3], paired = T))

#####################################################################

AC16_T28d24h_pairedttest = rbind(cut1,cut2,cut3,cdt1,cdt2,cdt3,nut1,nut2,nut3,ndt1,ndt2,ndt3)
write.csv(AC16_T28d24h_pairedttest, file="TimeSeries_AC16_T28d24h_pairedttest_data.csv")


################################################################################
################################################################################
################################################################################
################################################################################
#########U2OS 18d time series

mct1 = coldRefSeqGeneCountTable
cn1 = colnames(mct1)
cn2 = cn1[c(grep("U2OS_37d",cn1),grep("U2OS_18d", cn1))]

mat1 = mct1[,cn2]
cmat1 = mat1[,grep("CR", colnames(mat1))]
nmat1 = mat1[,grep("NR", colnames(mat1))]

crld = rlogtransformer(cmat1)
nrld = rlogtransformer(nmat1)

################################################################################
##################################################################################
######Averaging of repeats

dfcrld = as.data.frame(crld)
dfnrld = as.data.frame(nrld)

avcyt = repeataverager(dfcrld)
avnuc = repeataverager(dfnrld)

tscyt = timeseriesstandardiser(avcyt)
tsnuc = timeseriesstandardiser(avnuc)

##############################################################################
##################################################################################
######Get output from DESEQ

setwd(TWD)
setwd(TWD1)

nfx = read.csv("U2OS_37d_N_vs_U2OS_18d24h_N_DESEQ.csv", stringsAsFactors = F)
nfx1 = na.omit(nfx)

##########################################################################
##########################################################################
####Genes Up/Dn in Nucleus 18d24h U2OS cells

nfx1up = as.vector(nfx1[nfx1$padj < 0.05 & nfx1$log2FoldChange > 0,]$X)
nfx1dn = as.vector(nfx1[nfx1$padj < 0.05 & nfx1$log2FoldChange < 0,]$X)

cytup = tscyt[row.names(tscyt) %in% nfx1up,]
nucup = tsnuc[row.names(tsnuc) %in% nfx1up,]
cytdn = tscyt[row.names(tscyt) %in% nfx1dn,]
nucdn = tsnuc[row.names(tsnuc) %in% nfx1dn,]

#####################################################################

setwd(TWD)
setwd(TWD2)

colour1 = c("coral1","dodgerblue","dodgerblue","coral1")
ylab1 = "Standardised RNA levels"

png("TimeSeriesBoxPlot_CytoBehaviour_UpNuc18d24hU2OS.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(cytup, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Cytoplasm", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()

png("TimeSeriesBoxPlot_NucBehaviour_UpNuc18d24hU2OS.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(nucup, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Nucleus", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()

#####################################################################
#####################################################################
#####################################################################

png("TimeSeriesBoxPlot_CytoBehaviour_DnNuc18d24hU2OS.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(cytdn, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Cytoplasm", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()

png("TimeSeriesBoxPlot_NucBehaviour_DnNuc18d24hU2OS.png", height = 5, width = 3, units = "in", res = 600)
par(mar=c(1,5,3,1))
mp = boxplot(nucdn, range = 0, notch=T, lwd=2, xaxt="n", 
             las=1, col = colour1, yaxt="n", main = "Nucleus", cex.main = 2)
axis(1, at=1:length(mp$n), labels=rep("", length(mp$n)))
axis(2, las=2, cex.axis=1.5)
mtext(ylab1, side=2, line=3.5, cex=1.5)
box(lwd=2)
dev.off()

#####################################################################
#####################################################################
#########paired t-test

cytup1 = cytup[,c(1,3,4)]
cytdn1 = cytdn[,c(1,3,4)]
nucup1 = nucup[,c(1,3,4)]
nucdn1 = nucdn[,c(1,3,4)]

#####################################################################

cut1 = unlist(t.test(cytup1[,1], cytup1[,2], paired = T))
cut2 = unlist(t.test(cytup1[,2], cytup1[,3], paired = T))
cut3 = unlist(t.test(cytup1[,1], cytup1[,3], paired = T))

cdt1 = unlist(t.test(cytdn1[,1], cytdn1[,2], paired = T))
cdt2 = unlist(t.test(cytdn1[,2], cytdn1[,3], paired = T))
cdt3 = unlist(t.test(cytdn1[,1], cytdn1[,3], paired = T))

nut1 = unlist(t.test(nucup1[,1], nucup1[,2], paired = T))
nut2 = unlist(t.test(nucup1[,2], nucup1[,3], paired = T))
nut3 = unlist(t.test(nucup1[,1], nucup1[,3], paired = T))

ndt1 = unlist(t.test(nucdn1[,1], nucdn1[,2], paired = T))
ndt2 = unlist(t.test(nucdn1[,2], nucdn1[,3], paired = T))
ndt3 = unlist(t.test(nucdn1[,1], nucdn1[,3], paired = T))

#####################################################################

U2OS_T18d24h_pairedttest = rbind(cut1,cut2,cut3,cdt1,cdt2,cdt3,nut1,nut2,nut3,ndt1,ndt2,ndt3)
write.csv(U2OS_T18d24h_pairedttest, file="TimeSeries_U2OS_T18d24h_pairedttest_data.csv")
