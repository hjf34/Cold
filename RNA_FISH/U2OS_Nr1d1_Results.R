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

condition_repeat = as.vector(sapply(lf1, function(n) strsplit(n, "_Nr1d1")[[1]][1]))
condition = as.vector(sapply(lf1, function(n) strsplit(n, "_Rep")[[1]][1]))

df = data.frame(condition_repeat, condition, Nr1d1NuclearDots,Nr1d1CytoplasmicDots,NuclearArea,WholeCellArea)

df$CytoplasmicArea = df$WholeCellArea - df$NuclearArea
df$cba = df$Nr1d1CytoplasmicDots/df$CytoplasmicArea
df$nba = df$Nr1d1NuclearDots/df$NuclearArea
df1 = df
row.names(df1) = filenamestem
colnames(df1)[c(8,9)] = c("DotsPerCyt","DotsPerNuc")

repsimaged = data.frame(images=c(table(df1$condition_repeat)))
write.csv(repsimaged, file="U2OS_Nr1d1_RepsImaged.csv")
write.csv(df1, file="U2OS_Nr1d1_Dots_Per_Image.csv")

##########################################################################
##########################################################################

dba = df[,c("condition","cba","nba")]
dba$condition = as.vector(dba$condition)

dba1 = dba[!dba$condition == "U2OSNr1d1KO_T18d24h",]
dba1$condition1 = factor(dba1$condition, levels=names(table(dba1$condition))[c(5,3,4,1,2)])

dba2 = dba[dba$condition %in% c("U2OSNr1d1KO_T18d24h","U2OS_T18d24h"),]
dba2$condition1 = factor(dba2$condition, levels=names(table(dba2$condition)))

###########################################################################
########PER IMAGE PLOT AND STATS

mx1 = aggregate(dba1[,2:3], by=list(dba1$condition1), function(n) mean(n))
sx1 = aggregate(dba1[,2:3], by=list(dba1$condition1), function(n) sd(n)/sqrt(length(n)))

m1 = mx1[,3]; m2 = mx1[,2]
sem1 = sx1[,3]; sem2 = sx1[,2]

mx2 = aggregate(dba2[,2:3], by=list(dba2$condition1), function(n) mean(n))
sx2 = aggregate(dba2[,2:3], by=list(dba2$condition1), function(n) sd(n)/sqrt(length(n)))

m3 = mx2[,3]; m4 = mx2[,2]
sem3 = sx2[,3]; sem4 = sx2[,2]

#################################
#################################

barplotter = function(means,sems,colours){
  mx = means
  smx = sems 
  par(mar = c(1,4,1,1))
  mp = barplot(mx, cex.axis=1.5,col=colours,
               ylim = c(0,1.1*max(mx+smx, na.rm=T)),xaxt="n", las=1, main = "", ylab = "")
  
  subseg1 = as.numeric(mx) -as.numeric(smx)
  subseg2 = as.numeric(mx) +as.numeric(smx)
  subseg1[subseg1 < 0] = 0
  segments(mp, subseg1, mp, subseg2, lwd=2)
  segments(mp-0.2, subseg1, mp+0.2, subseg1, lwd=2)
  segments(mp-0.2, subseg2, mp+0.2, subseg2, lwd=2)
  box(lwd = 2)
}

colours1 = c("gold","skyblue","pink","dodgerblue3","purple")
colours2 = c("grey","white")

#########################NUCLEAR
png("U2OS_Nr1d1NucPerImageMeanSEM.png", height = 3, width = 2.5, units = "in", res = 600)
barplotter(m1,sem1,colours1)
dev.off()

#########################CYTOPLASMIC
png("U2OS_Nr1d1CytPerImageMeanSEM.png", height = 3, width = 2.5, units = "in", res = 600)
barplotter(m2,sem2,colours1)
dev.off()

#########################NUCLEAR
png("U2OS_WTvKO_Nr1d1NucPerImageMeanSEM.png",height = 3, width = 2, units = "in", res = 600)
barplotter(m3,sem3,colours2)
dev.off()

#########################CYTOPLASMIC
png("U2OS_WTvKO_Nr1d1CytPerImageMeanSEM.png",height = 3, width = 2, units = "in", res = 600)
barplotter(m4,sem4,colours2)
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
row.names(THSD) = c(paste("Nuc_",row.names(THSD)[1:10], sep=""),paste("Cyt_",row.names(THSD)[11:20], sep=""))

write.csv(THSD, file="U2OS_Nr1d1_TukeyHSDdata.csv")
write.csv(aovdata, file="U2OS_Nr1d1_ANOVAdata.csv")

#############################################################################
#############################################################################

nvec1 = dba2[dba2$condition1=="U2OS_T18d24h",3]
nvec2 = dba2[dba2$condition1=="U2OSNr1d1KO_T18d24h",3]

cvec1 = dba2[dba2$condition1=="U2OS_T18d24h",2]
cvec2 = dba2[dba2$condition1=="U2OSNr1d1KO_T18d24h",2]

ttestpval = data.frame(comparison=c("Nuc_WTvKO","Cyt_WTvKO"),pvalue = c(t.test(nvec1,nvec2)[[3]],t.test(cvec1,cvec2)[[3]]))
write.csv(ttestpval, file="U2OS_Nr1d1_ttestpval_WTvKO.csv")

