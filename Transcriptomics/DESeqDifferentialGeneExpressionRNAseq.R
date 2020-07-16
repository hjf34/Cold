###Author: Harry Fischl
###Differential expression analysis of RNA-seq data using DESeq2 R package 
###############################Work directories

TWD1 = "./DESeq/DifferentialGeneExpressionGeneLists/"

###Load count table of RNA-seq counts per gene for all samples
load("./CountTableGeneration/CountTables/coldRefSeqGeneCountTable.Rda")
gct = coldRefSeqGeneCountTable

library(DESeq2)

##################################################################################
#### Select two temperatures to compare and whether to compare cytoplasmic (C) or nuclear (N) fractions 
deseqer = function(temperature1, temperature2, CorN, l2fcThreshold=1){
  tx = as.character(c(temperature1, temperature2))
  tx = paste(tx, CorN, sep="_")
  temperature_options = colnames(gct)
  temperature_options_split = as.vector(sapply(temperature_options, function(n) strsplit(n, "R")[[1]][1]))
  ################################################################
  temp1_x = unique(temperature_options_split[grep(tx[1],temperature_options_split)])
  temp2_x = unique(temperature_options_split[grep(tx[2],temperature_options_split)])
  ################################################################
  ###Get column numbers of all repeats of each temperature condition
  colnum1 = grep(tx[1], temperature_options_split)
  colnum2 = grep(tx[2], temperature_options_split)
  colnums = c(colnum1,colnum2)
  mat1 = gct[,colnums]
  #Creating ddsMat object for DESeq
  c1 = colnames(mat1)
  c2 = as.vector(sapply(c1, function(n) strsplit(n, "R")[[1]][1]))
  conditions = c2
  coldata <- data.frame(conditions)
  row.names(coldata) = c1
  ########################################################################
  ddsMat <- DESeqDataSetFromMatrix(countData = mat1,
                                   colData = coldata,
                                   design = ~ conditions)
  #DESEQ results table
  res = results(DESeq(ddsMat,quiet=T), 
                contrast=c("conditions",temp1_x,temp2_x),
                lfcThreshold=log2(l2fcThreshold),altHypothesis="greaterAbs")
  return(res)
}

##################################################################################
#####Write DESEQ results table to csv file
csvwriter = function(deseqeroutput){
  df = data.frame(deseqeroutput)
  df1 = df
  df1$log2FoldChange = -df$log2FoldChange
  m1 = mcols(deseqeroutput)[[2]][2]
  c1vc2 = paste(strsplit(m1, " ")[[1]][6], strsplit(m1, " ")[[1]][8], sep="_vs_")
  c1vc2filename = paste(c1vc2, "_DESEQ.csv", sep="")
  write.csv(df1, file= c1vc2filename)
}

setwd(TWD1)
csvwriter(deseqer("AC16_37d","AC16_18d5h","C"))
csvwriter(deseqer("AC16_37d","AC16_18d5h37d90min","C"))
csvwriter(deseqer("AC16_37d","AC16_18d10h","C"))
csvwriter(deseqer("AC16_37d","AC16_18d24h","C"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d2h","C"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d5h","C"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d10h","C"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d24h","C"))
csvwriter(deseqer("AC16_37d","AC16_8d24h","C"))
csvwriter(deseqer("AC16_37d","AC16_28d24h","C"))
csvwriter(deseqer("AC16Spiked_37d","AC16Spiked_18d24h","C"))
csvwriter(deseqer("U2OS_37d","U2OS_18d5h","C"))
csvwriter(deseqer("U2OS_37d","U2OS_18d24h","C"))
csvwriter(deseqer("U2OS_37d","U2OS_18d24h37d2h","C"))

csvwriter(deseqer("AC16_37d","AC16_18d5h","N"))
csvwriter(deseqer("AC16_37d","AC16_18d5h37d90min","N"))
csvwriter(deseqer("AC16_37d","AC16_18d10h","N"))
csvwriter(deseqer("AC16_37d","AC16_18d24h","N"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d2h","N"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d5h","N"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d10h","N"))
csvwriter(deseqer("AC16_37d","AC16_18d24h37d24h","N"))
csvwriter(deseqer("AC16_37d","AC16_8d24h","N"))
csvwriter(deseqer("AC16_37d","AC16_28d24h","N"))
csvwriter(deseqer("AC16Spiked_37d","AC16Spiked_18d24h","N"))
csvwriter(deseqer("U2OS_37d","U2OS_18d5h","N"))
csvwriter(deseqer("U2OS_37d","U2OS_18d24h","N"))
csvwriter(deseqer("U2OS_37d","U2OS_18d24h37d2h","N"))
