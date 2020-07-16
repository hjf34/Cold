###Author: Harry Fischl
###Relationship between the log2 fold change in RNA level for all genes upon transfer of AC16 cells from 37°C to 28°C or 18°C, 28°C or 8°C, and 18°C or 8°C for 24h in both the nucleus and cytoplasm
###Figure 2C
###############################Work directories

TWD = getwd()
TWD1 = "./DESeq/DifferentialGeneExpressionGeneLists/"
TWD2 = "./Rplots/"
setwd(TWD1)

#######Input of DESEQ output files
xc1 = read.csv("AC16_37d_C_vs_AC16_28d24h_C_DESEQ.csv", stringsAsFactors = F)
yc1 = read.csv("AC16_37d_C_vs_AC16_18d24h_C_DESEQ.csv", stringsAsFactors = F)
zc1 = read.csv("AC16_37d_C_vs_AC16_8d24h_C_DESEQ.csv", stringsAsFactors = F)

xn1 = read.csv("AC16_37d_N_vs_AC16_28d24h_N_DESEQ.csv", stringsAsFactors = F)
yn1 = read.csv("AC16_37d_N_vs_AC16_18d24h_N_DESEQ.csv", stringsAsFactors = F)
zn1 = read.csv("AC16_37d_N_vs_AC16_8d24h_N_DESEQ.csv", stringsAsFactors = F)

cdf = cbind(xc1[,c("X","log2FoldChange","padj")],yc1[,c("log2FoldChange","padj")],zc1[,c("log2FoldChange","padj")])
ndf = cbind(xn1[,c("X","log2FoldChange","padj")],yn1[,c("log2FoldChange","padj")],zn1[,c("log2FoldChange","padj")])

colnames(cdf) = c("gene","fc28","p28","fc18","p18","fc8","p8")
colnames(ndf) = c("gene","fc28","p28","fc18","p18","fc8","p8")

##########Exclusion of genes with NA values for either condition being compared
tcdf1 = na.omit(cdf[,c(1,2,3,4,5)])
tcdf2 = na.omit(cdf[,c(1,2,3,6,7)])
tcdf3 = na.omit(cdf[,c(1,4,5,6,7)])

tndf1 = na.omit(ndf[,c(1,2,3,4,5)])
tndf2 = na.omit(ndf[,c(1,2,3,6,7)])
tndf3 = na.omit(ndf[,c(1,4,5,6,7)])

##############################
###Gives different values for genes significantly differentially expressed in one condition or both conditions

valuer = function(tdf){
  tdfx = tdf
  val1 = rep(NA, length(tdfx$gene))
  for(a1 in 1:length(val1)){
    vec1 = tdfx[a1,]
    if(vec1[3] >= 0.05 & vec1[5] >= 0.05){val1[a1] = 4; next}
    if(vec1[3] < 0.05 & vec1[5] < 0.05 & vec1[2]/vec1[4] < 0){val1[a1] = 0; next}
    if(vec1[3] < 0.05 & vec1[5] < 0.05 & vec1[2]/vec1[4] >= 0){val1[a1] = 1; next}
    if(vec1[3] < 0.05 & vec1[5] >= 0.05 ){val1[a1] = 2; next}
    if(vec1[3] >= 0.05 & vec1[5] < 0.05){val1[a1] = 3; next}
  }
  tdfx$val = val1
  fdf1 = tdfx[order(-tdfx$val),]
  return(fdf1)
}
##############

vcdf1 = valuer(tcdf1)
vcdf2 = valuer(tcdf2)
vcdf3 = valuer(tcdf3)

vndf1 = valuer(tndf1)
vndf2 = valuer(tndf2)
vndf3 = valuer(tndf3)

##############

densityplotter = function(valueroutput){
  q1 = valueroutput
  
  x <- densCols(q1[,2],q1[,4], nbin=700,colramp=colorRampPalette(c("black", "white")))
  q1$dens <- col2rgb(x)[1,] + 1L
  cols <-  colorRampPalette(c("gray70","black"))(256)
  q1$col <- cols[q1$dens]
  q2 = q1[order(q1$dens),]
  
  ylim1 = c(-3,3)
  xlim1 = c(-3,3)
  
  par(xpd=F)
  par(mar=c(1,1,1,1))
  plot(q2[,2],q2[,4], ylim=ylim1, xlim=xlim1,col=q2$col,pch=20,cex=0.3)
  return(q2)
}

cortdf = data.frame(matrix(NA, ncol =4, nrow = 6))
row.names(cortdf) = c("C28d24hvC18d24h","C28d24hvC8d24h","C18d24hvC8d24","N28d24hvN18d24h","N28d24hvN8d24h","N18d24hvN8d24h")
colnames(cortdf) = c("t","doff","p","r")
#############################################################################
####CYT
setwd(TWD)
setwd(TWD2)
#28 v 18

png("TemperatureComparison_28v18cyt.png", width = 3, height = 3, units = "in",res = 200,pointsize=7)
q2 = densityplotter(vcdf1)

points(q2[q2$val == 2,2],q2[q2$val == 2,4], col="firebrick1", lwd=2, cex=0.75)
points(q2[q2$val == 3,2],q2[q2$val == 3,4], col="dodgerblue3", lwd=2, cex=0.75)
points(q2[q2$val == 1,2],q2[q2$val == 1,4], col="orange", lwd=2, cex=1.5, pch=2)
points(q2[q2$val == 0,2],q2[q2$val == 0,4], col="brown", lwd=2, cex=1.5, pch=2)
box(lwd=2)
dev.off()
#############################################################################

cort = cor.test(q2[,2],q2[,4], exact = T)
cortvec = c(cort[[1]],cort[[2]],cort[[3]],cort[[4]])
cortdf[1,] = cortvec
#############################################################################


#28 v 8
png("TemperatureComparison_28v8cyt.png",  width = 3, height = 3, units = "in",res = 200,pointsize=7)
q2 = densityplotter(vcdf2)

points(q2[q2$val == 2,2],q2[q2$val == 2,4], col="firebrick1", lwd=2, cex=0.75)
points(q2[q2$val == 3,2],q2[q2$val == 3,4], col="green3", lwd=2, cex=0.75)
points(q2[q2$val == 1,2],q2[q2$val == 1,4], col="orange", lwd=2, cex=1.5, pch=2)
points(q2[q2$val == 0,2],q2[q2$val == 0,4], col="brown", lwd=2, cex=1.5, pch=2)
box(lwd=2)
dev.off()
#############################################################################

cort = cor.test(q2[,2],q2[,4], exact = T)
cortvec = c(cort[[1]],cort[[2]],cort[[3]],cort[[4]])
cortdf[2,] = cortvec
#############################################################################

#18 v 8
png("TemperatureComparison_18v8cyt.png", width = 3, height = 3, units = "in",res = 200,pointsize=7)
q2 = densityplotter(vcdf3)

points(q2[q2$val == 2,2],q2[q2$val == 2,4], col="dodgerblue3", lwd=2, cex=0.75)
points(q2[q2$val == 3,2],q2[q2$val == 3,4], col="green3", lwd=2, cex=0.75)
points(q2[q2$val == 1,2],q2[q2$val == 1,4], col="orange", lwd=2, cex=1.5, pch=2)
points(q2[q2$val == 0,2],q2[q2$val == 0,4], col="brown", lwd=2, cex=1.5, pch=2)
box(lwd=2)
dev.off()
#############################################################################

cort = cor.test(q2[,2],q2[,4], exact = T)
cortvec = c(cort[[1]],cort[[2]],cort[[3]],cort[[4]])
cortdf[3,] = cortvec
#############################################################################

#############################################################################
####NUC
#28 v 18

png("TemperatureComparison_28v18nuc.png", width = 3, height = 3, units = "in",res = 200,pointsize=7)
q2 = densityplotter(vndf1)

points(q2[q2$val == 2,2],q2[q2$val == 2,4], col="firebrick1", lwd=2, cex=0.75)
points(q2[q2$val == 3,2],q2[q2$val == 3,4], col="dodgerblue3", lwd=2, cex=0.75)
points(q2[q2$val == 1,2],q2[q2$val == 1,4], col="orange", lwd=2, cex=1.5, pch=2)
points(q2[q2$val == 0,2],q2[q2$val == 0,4], col="brown", lwd=2, cex=1.5, pch=2)

box(lwd=2)
dev.off()
#############################################################################

cort = cor.test(q2[,2],q2[,4], exact = T)
cortvec = c(cort[[1]],cort[[2]],cort[[3]],cort[[4]])
cortdf[4,] = cortvec
#############################################################################

#28 v 8
png("TemperatureComparison_28v8nuc.png",  width = 3, height = 3, units = "in",res = 200,pointsize=7)
q2 = densityplotter(vndf2)

points(q2[q2$val == 2,2],q2[q2$val == 2,4], col="firebrick1", lwd=2, cex=0.75)
points(q2[q2$val == 3,2],q2[q2$val == 3,4], col="green3", lwd=2, cex=0.75)
points(q2[q2$val == 1,2],q2[q2$val == 1,4], col="orange", lwd=2, cex=1.5, pch=2)
points(q2[q2$val == 0,2],q2[q2$val == 0,4], col="brown", lwd=2, cex=1.5, pch=2)

box(lwd=2)
dev.off()
#############################################################################

cort = cor.test(q2[,2],q2[,4], exact = T)
cortvec = c(cort[[1]],cort[[2]],cort[[3]],cort[[4]])
cortdf[5,] = cortvec
#############################################################################

#18 v 8
png("TemperatureComparison_18v8nuc.png", width = 3, height = 3, units = "in",res = 200,pointsize=7)
q2 = densityplotter(vndf3)

points(q2[q2$val == 2,2],q2[q2$val == 2,4], col="dodgerblue3", lwd=2, cex=0.75)
points(q2[q2$val == 3,2],q2[q2$val == 3,4], col="green3", lwd=2, cex=0.75)
points(q2[q2$val == 1,2],q2[q2$val == 1,4], col="orange", lwd=2, cex=1.5, pch=2)
points(q2[q2$val == 0,2],q2[q2$val == 0,4], col="brown", lwd=2, cex=1.5, pch=2)

box(lwd=2)
dev.off()
#############################################################################

cort = cor.test(q2[,2],q2[,4], exact = T)
cortvec = c(cort[[1]],cort[[2]],cort[[3]],cort[[4]])
cortdf[6,] = cortvec
#############################################################################

####Statistics for correlation analysis
colnames(cortdf) = c("t_statistic","degrees_of_freedom","p_value","pearsons_r")
write.csv(cortdf, file="TemperatureComparison_CorrelationStatistics.csv")
