###Author: Harry Fischl
###Mean ratio of D. melanogaster (fly) aligned polyA site (PAS) reads to human aligned PAS reads for nuclear and cytoplasmic samples prepared from AC16 cells harvested after exposure to the three different conditions shown and spiked with a known ratio of D. melanogaster Schneider 2 (S2) cells.
###Figure EV2C
###############################Work directories

TWD = getwd()
TWD1 = "./Rplots/"
  
load("./CountTableGeneration/CountTables/coldCountTableRaw.Rda")
load("./CountTableGeneration/Drosophila/CountTable/FlyCountTableRaw1.Rda")

cct = coldCountTableRaw
fct = FlyCountTableRaw1

cct1 = cct[,names(fct)]

humancounts = colSums(cct1[,7:length(cct1)])
flycounts = colSums(fct[,7:length(fct)])

ratios = flycounts/humancounts

###Repeats 1 and 2 were spiked in with 5x fewer Drosophila cells to AC16 cells than Repeat 3 and 4
ratios[grep("R1", names(ratios))] = ratios[grep("R1", names(ratios))]*5
ratios[grep("R2", names(ratios))] = ratios[grep("R2", names(ratios))]*5

rats = ratios

ratsmat = matrix(NA, nrow=4, ncol=6)
ratsmat[,1] = rats[17:20]
ratsmat[,2] = rats[1:4]
ratsmat[,3] = rats[9:12]
ratsmat[,4] = rats[21:24]
ratsmat[,5] = rats[5:8]
ratsmat[,6] = rats[13:16]

cyt= ratsmat[,1:3]
nuc= ratsmat[,4:6]

#########################################################################
colSems = function(dataframe){
  semvec = rep(NA, length(dataframe))
  for(a1 in 1:length(dataframe)){
    vec = dataframe[,a1]
    vec1 = vec[!is.na(vec)]
    semvec[a1] = sd(vec1)/sqrt(length(vec1))
  }
  return(semvec)
}

setwd(TWD1)

png("SpikeInRatio_Cyt_BarChart.png", height = 4, width = 3, units = "in", res = 600)
par(mar=c(1,4,1,1))
m1 = colMeans(data.frame(cyt))
s1 = colSems(data.frame(cyt))
mp = barplot(m1, main="", col=c("white"),
             ylab = "", las=1,cex.axis=1.5, ylim= c(0,1.2*max(cyt)), border=F)
mp3 = rep(mp[3],4)
mp3[1] = mp3[1]-0.05
mp3[4] = mp3[4]+0.05

par(lwd=2)
mp = barplot(m1, main="", col=c("gold","dodgerblue3","purple"), angle=45, density=20,
             ylab = "", las=1,cex.axis=1.5, ylim= c(0,1.2*max(cyt)), border="black")
par(lwd=1)

subseg1 = m1-s1
subseg2 = m1+s1
segments(mp, subseg1, mp, subseg2, lwd=1)
segments(mp-0.2, subseg1, mp+0.2, subseg1, lwd=1)
segments(mp-0.2, subseg2, mp+0.2, subseg2, lwd=1)

points(x=rep(mp[1],4), y= cyt[,1], pch=20, col="gold", lwd =2)
points(x=rep(mp[1],4), y= cyt[,1], pch=21, lwd =2)
points(x=rep(mp[2],4), y= cyt[,2], pch=20, col="dodgerblue3", lwd =2)
points(x=rep(mp[2],4), y= cyt[,2], pch=21, lwd =2)
points(x=mp3, y= cyt[,3], pch=20, col="purple", lwd =2)
points(x=mp3, y= cyt[,3], pch=21, lwd =2)

box(lwd=2)

dev.off()

###########################################################

png("SpikeInRatio_Nuc_BarChart.png", height = 4, width = 3, units = "in", res = 600)
par(mar=c(1,4,1,1))
m1 = colMeans(data.frame(nuc))
s1 = colSems(data.frame(nuc))
mp = barplot(m1, main="", col=c("white"),
             ylab = "", las=1,cex.axis=1.5, ylim= c(0,1.2*max(nuc)), border=F)
mp3 = rep(mp[3],4)
mp3[1] = mp3[1]-0.05
mp3[4] = mp3[4]+0.05

par(lwd=2)
mp = barplot(m1, main="", col=c("gold","dodgerblue3","purple"), angle=45, density=20,
             ylab = "", las=1,cex.axis=1.5, ylim= c(0,1.2*max(nuc)), border="black")
par(lwd=1)

subseg1 = m1-s1
subseg2 = m1+s1
segments(mp, subseg1, mp, subseg2, lwd=1)
segments(mp-0.2, subseg1, mp+0.2, subseg1, lwd=1)
segments(mp-0.2, subseg2, mp+0.2, subseg2, lwd=1)

points(x=rep(mp[1],4), y= nuc[,1], pch=20, col="gold", lwd =2)
points(x=rep(mp[1],4), y= nuc[,1], pch=21, lwd =2)
points(x=rep(mp[2],4), y= nuc[,2], pch=20, col="dodgerblue3", lwd =2)
points(x=rep(mp[2],4), y= nuc[,2], pch=21, lwd =2)
points(x=mp3, y= nuc[,3], pch=20, col="purple", lwd =2)
points(x=mp3, y= nuc[,3], pch=21, lwd =2)

box(lwd=2)

dev.off()


condition = as.factor(c(rep("A1",4),rep("A2",4),rep("A3",4)))
ratvec1 = c(cyt[,1],cyt[,2],cyt[,3])
ratvec2 = c(nuc[,1],nuc[,2],nuc[,3])

summary(aov(ratvec1 ~ condition))[[1]][["Pr(>F)"]][1]
summary(aov(ratvec2 ~ condition))[[1]][["Pr(>F)"]][1]

cytoratios = data.frame(cyt)
nuclratios = data.frame(nuc)
names(cytoratios) = c("T37d","T18d24h","T18d24h37d2h")
names(nuclratios) = c("T37d","T18d24h","T18d24h37d2h")

write.csv(cytoratios, file="SpikeInRatios_cyt.csv")
write.csv(nuclratios, file="SpikeInRatios_nuc.csv")

write.csv(data.frame(summary(aov(ratvec1 ~ condition))[[1]]), file="SpikeInRatios_anovastats_cyt.csv")
write.csv(data.frame(summary(aov(ratvec2 ~ condition))[[1]]), file="SpikeInRatios_anovastats_nuc.csv")
