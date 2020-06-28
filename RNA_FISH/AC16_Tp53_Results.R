setwd("/path/to/Results.csv/directory/")
lf = list.files()
lf1 = lf[grep("Results.csv", lf)]
filenamestem = as.vector(sapply(lf1, function(n) strsplit(n, ".dv")[[1]][1]))

Nr1d1NuclearDots = rep(NA, length(lf1))
Nr1d1CytoplasmicDots = rep(NA, length(lf1))
NuclearArea = rep(NA, length(lf1))
WholeCellArea = rep(NA, length(lf1))

for(a1 in 1:length(lf1)){
  csv1 = read.csv((lf1[a1]), stringsAsFactors = F)
  NuclearArea[a1] = csv1$Area[1]
  WholeCellArea[a1] = csv1$Area[2]
  Nr1d1NuclearDots[a1] = csv1$Count[3]
  Nr1d1CytoplasmicDots[a1] = csv1$Count[4]
}

condition_repeat = as.vector(sapply(lf1, function(n) strsplit(n, "_Poor")[[1]][1]))
condition = as.vector(sapply(lf1, function(n) strsplit(n, "_Rep")[[1]][1]))

df = data.frame(condition_repeat, condition, Nr1d1NuclearDots,Nr1d1CytoplasmicDots,NuclearArea,WholeCellArea)

df$CytoplasmicArea = df$WholeCellArea - df$NuclearArea
df$cba = df$Nr1d1CytoplasmicDots/df$CytoplasmicArea
df$nba = df$Nr1d1NuclearDots/df$NuclearArea
df1 = df
row.names(df1) = filenamestem
colnames(df1)[c(8,9)] = c("DotsPerCyt","DotsPerNuc")

repsimaged = data.frame(images=c(table(df1$condition_repeat)))
write.csv(repsimaged, file="AC16_Tp53_RepsImaged.csv")
write.csv(df1, file="AC16_Tp53_Dots_Per_Image.csv")

##########################################################################
##########################################################################

dba = df[,c("condition","cba","nba")]
dba$condition = as.vector(dba$condition)

dba1 = dba
dba1$condition1 = factor(dba1$condition, levels=names(table(dba1$condition))[c(3,1,2)])

###########################################################################
########PER IMAGE PLOT AND STATS

mx1 = aggregate(dba1[,2:3], by=list(dba1$condition1), function(n) mean(n))
sx1 = aggregate(dba1[,2:3], by=list(dba1$condition1), function(n) sd(n)/sqrt(length(n)))

m1 = mx1[,3]; m2 = mx1[,2]
sem1 = sx1[,3]; sem2 = sx1[,2]

#################################
#################################

barplotter = function(means,sems,colours){
  mx = means
  smx = sems 
  par(mar = c(1,5,1,1))
  mp = barplot(mx, cex.axis=1.5,col=colours,
               ylim = c(0,1.2*max(mx+smx, na.rm=T)),xaxt="n", las=1, main = "", ylab = "")
  
  subseg1 = as.numeric(mx) -as.numeric(smx)
  subseg2 = as.numeric(mx) +as.numeric(smx)
  subseg1[subseg1 < 0] = 0
  segments(mp, subseg1, mp, subseg2, lwd=2)
  segments(mp-0.2, subseg1, mp+0.2, subseg1, lwd=2)
  segments(mp-0.2, subseg2, mp+0.2, subseg2, lwd=2)
  box(lwd = 2)
}

colours1 = c("gold","dodgerblue3","purple")

#########################NUCLEAR
png("AC16_Tp53NucPerImageMeanSEM.png", height = 3, width = 2.5, units = "in", res = 600)
barplotter(m1,sem1,colours1)
dev.off()

#########################CYTOPLASMIC
png("AC16_Tp53CytPerImageMeanSEM.png", height = 3, width = 2.5, units = "in", res = 600)
barplotter(m2,sem2,colours1)
dev.off()

###########################################################################
###########################################################################
###########################################################################

summary(aov(cba~condition1, data=dba1))
cTHSD = TukeyHSD(aov(cba~condition1, data=dba1))$condition1
summary(aov(nba~condition1, data=dba1))
nTHSD = TukeyHSD(aov(nba~condition1, data=dba1))$condition1

aovdata = rbind(data.frame(summary(aov(nba~condition1, data=dba1))[[1]]),data.frame(summary(aov(cba~condition1, data=dba1))[[1]]))
row.names(aovdata) = c("Nuc_Condition","Nuc_Residuals","Cyt_Condition","Cyt_Residuals")

THSD = rbind(nTHSD,cTHSD)
row.names(THSD) = c(paste("Nuc_",row.names(THSD)[1:3], sep=""),paste("Cyt_",row.names(THSD)[4:6], sep=""))

write.csv(THSD, file="AC16_Tp53_TukeyHSDdata.csv")
write.csv(aovdata, file="AC16_Tp53_ANOVAdata.csv")

#############################################################################
#############################################################################

