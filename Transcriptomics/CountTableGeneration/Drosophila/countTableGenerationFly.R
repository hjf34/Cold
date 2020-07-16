###Author: Harry Fischl
###Counts number of reads in positive and negative strand narrowed bigwig files of 3' end RNA-seq data mapped to dm6 drosophila genome overlapping with fly poly A site annotations
###############################Work directories

TWD = getwd()
TWD1 = "./dm6bigwigFiles/"
TWD2 = "./CountTable/"

###############################Library input

library(GenomicRanges)
library(GenomicAlignments)
library("BSgenome.Dmelanogaster.UCSC.dm6")
library(rtracklayer)

seq_info_object1 = SeqinfoForBSGenome(Dmelanogaster)

#########polyA Sites from doi:10.1261/rna.062661.117
###Transcription Elongation Rate Has a Tissue-Specific Impact on Alternative Cleavage and Polyadenylation in Drosophila melanogaster. RNA, 2017. Xiaochuan Liu, Jaime Freitas, Dinghai Zheng, Marta S Oliveira, Mainul Hoque, Torcato Martins, Telmo Henriques, Bin Tian, Alexandra Moreira
##Each PAS was extended 20 nt 3’ and 200 nt 5’ from the site of cleavage and those that overlapped on the same strand after extension were combined into a single PAS annotation.
##These extended and combined PAS are saved as a GenomicRanges object in file "reducedFlyPas.Rda"

load("./polyASites/reducedFlyPas.Rda")
posReducedFlyPas = reducedFlyPas[strand(reducedFlyPas) == "+"]
negReducedFlyPas = reducedFlyPas[strand(reducedFlyPas) == "-"]

countTable = data.frame(c(posReducedFlyPas,negReducedFlyPas))

##################################################################################################
########COUNTS TABLE FROM BIGWIGS
#################################################
####Set directory to file containing bigwigs of narrowed reads aligned to dm6 processed from bamfiles
####pos.bw and neg.bw for all samples can be downloaded from GEO (https://www.ncbi.nlm.nih.gov/geo/) under accession code GSE137003

setwd(TWD)
setwd(TWD1)

lf = list.files(pattern="pos.bw")

for(a1 in 1:length(lf)){
  posfile = lf[a1]
  filestem = strsplit(posfile, "pos.bw")[[1]][1]
  negfile = paste0(filestem, "neg.bw")
  ########################Import positive and negative strand bigwigs
  p1 = import(posfile)
  n1 = import(negfile)
  ########################Assign levels and seq info to bigwigs
  seqlevels(p1) = seqlevels(seq_info_object1)
  seqinfo(p1) = seq_info_object1
  seqlevels(n1) = seqlevels(seq_info_object1)
  seqinfo(n1) = seq_info_object1
  strand(p1) = "+"
  strand(n1) = "-"
  ########################Intersect with reducedFlyPas Granges
  ###################PositiveStrand
  pcov = coverage(p1, weight=p1$score)
  v1 = Views(pcov, posReducedFlyPas)
  s1 = lapply(v1, function(n) sum(n))
  pc = as.vector(unlist(s1))
  ###################NegativeStrand
  ncov = coverage(n1, weight=n1$score)
  v1 = Views(ncov, negReducedFlyPas)
  s1 = lapply(v1, function(n) sum(n))
  nc = as.vector(unlist(s1))
  ##############################################
  countTable$new1 = c(pc,nc)
  names(countTable)[which(names(countTable) == "new1")] = filestem
}

setwd(TWD)
setwd(TWD2)
FlyCountTableRaw1 = countTable
#save(FlyCountTableRaw1, file="FlyCountTableRaw1.Rda")


