###Author: Harry Fischl
###Density scatter graphs of the relationship between the log2 fold change in RNA level in the nucleus and that in the cytoplasm for all genes upon transfer of AC16 cells from 37°C to different cold temperatures for 24h.
###Figure 2B and Figure EV2D-E,G-H
###############################Work directories

TWD = getwd()
TWD1 = "./DESeq/DifferentialGeneExpressionGeneLists/"
TWD2 = "./Rplots/"
setwd(TWD)
setwd(TWD1)

lf1 = list.files(pattern=".csv")

################################################################################
################################################################################

corplot = function(csvfile1,csvfile2){
  setwd(TWD)
  setwd(TWD1)
  r1 = read.csv(csvfile1, stringsAsFactors = F)
  s1 = read.csv(csvfile2, stringsAsFactors = F)
  setwd(TWD)
  setwd(TWD2)
  r2 = data.frame(r1)
  s2 = data.frame(s1)
  r3 = r2[r2$baseMean >= 30,]
  s3 = s2[s2$baseMean >= 30,]
  
  r4 = r3[row.names(r3) %in% row.names(s3),]
  s4 = s3[row.names(s3) %in% row.names(r3),]
  
  xlabX = as.vector(sapply(csvfile1, function(n) strsplit(n, "_")[[1]]))
  xlab1 = xlabX[2]
  xlab2 = xlabX[6]
  xfract = xlabX[3]
  xlab = paste(xlab1, xlab2, sep=" v ")
  xlab = paste(xlabX[1], xlab, sep=" ")
  labelx = paste(xlab,sub("000",xfract,  "(000)"), sep=" ")
  
  ylabX = as.vector(sapply(csvfile2, function(n) strsplit(n, "_")[[1]]))
  ylab1 = ylabX[2]
  ylab2 = ylabX[6]
  yfract = ylabX[3]
  ylab = paste(ylab1, ylab2, sep=" v ")
  ylab = paste(ylabX[1], ylab, sep=" ")
  labely = paste(ylab,sub("000",yfract,  "(000)"), sep=" ")
  
  q1 = data.frame(gene = r4$X, l2fc1 = r4$log2FoldChange, l2fc2 = s4$log2FoldChange)
  x <- densCols(q1$l2fc1,q1$l2fc2, nbin=700,colramp=colorRampPalette(c("black", "white")))
  q1$dens <- col2rgb(x)[1,] + 1L
  cols <-  colorRampPalette(c("black","blue","purple","red","orange","yellow"))(256)
  q1$col <- cols[q1$dens]
  
  ylim1 = c(-3,3)
  xlim1 = c(-3,3)
  
  par(xpd=F)
  par(mar=c(5,5,1,1))
  plot(l2fc2~l2fc1,data= q1[order(q1$dens),], ylim=ylim1, xlim=xlim1,col=col,pch=20,cex=0.3,
       xlab=labelx,ylab=labely, las=1, cex.axis=1.5, cex.lab=1)
  
  box(lwd=2)
}

################################################################################
################################################################################

png("NucVCyt_AC16_37d_v_18d5h_C_v_AC16_37d_v_18d5h_N.png",
    height = 5, width = 5, units = "in", res = 600)
corplot(lf1[grep("AC16_37d_C_vs_AC16_18d5h_C", lf1)],
        lf1[grep("AC16_37d_N_vs_AC16_18d5h_N", lf1)])
dev.off()

png("NucVCyt_AC16_37d_v_18d10h_C_v_AC16_37d_v_18d10h_N.png",
    height = 5, width = 5, units = "in", res = 600)
corplot(lf1[grep("AC16_37d_C_vs_AC16_18d10h_C", lf1)],
        lf1[grep("AC16_37d_N_vs_AC16_18d10h_N", lf1)])
dev.off()

png("NucVCyt_AC16_37d_v_18d24h_C_v_AC16_37d_v_18d24h_N.png",
    height = 5, width = 5, units = "in", res = 600)
corplot(lf1[grep("AC16_37d_C_vs_AC16_18d24h_C", lf1)],
        lf1[grep("AC16_37d_N_vs_AC16_18d24h_N", lf1)])
dev.off()

png("NucVCyt_AC16_37d_v_28d24h_C_v_AC16_37d_v_28d24h_N.png",
    height = 5, width = 5, units = "in", res = 600)
corplot(lf1[grep("AC16_37d_C_vs_AC16_28d24h_C", lf1)],
        lf1[grep("AC16_37d_N_vs_AC16_28d24h_N", lf1)])
dev.off()

png("NucVCyt_AC16_37d_v_8d24h_C_v_AC16_37d_v_8d24h_N.png",
    height = 5, width = 5, units = "in", res = 600)
corplot(lf1[grep("AC16_37d_C_vs_AC16_8d24h_C", lf1)],
        lf1[grep("AC16_37d_N_vs_AC16_8d24h_N", lf1)])
dev.off()

png("NucVCyt_U2OS_37d_v_18d24h_C_v_U2OS_37d_v_18d24h_N.png",
    height = 5, width = 5, units = "in", res = 600)
corplot(lf1[grep("U2OS_37d_C_vs_U2OS_18d24h_C", lf1)],
        lf1[grep("U2OS_37d_N_vs_U2OS_18d24h_N", lf1)])
dev.off()

png("NucVCyt_AC16Spiked_37d_v_18d24h_C_v_AC16Spiked_37d_v_18d24h_N.png",
    height = 5, width = 5, units = "in", res = 600)
corplot(lf1[grep("AC16Spiked_37d_C_vs_AC16Spiked_18d24h_C", lf1)],
        lf1[grep("AC16Spiked_37d_N_vs_AC16Spiked_18d24h_N", lf1)])
dev.off()



##################################################################################
##################################################################################
##################################################################################
