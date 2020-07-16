###Author: Harry Fischl
###Counts number of reads in positive and negative strand narrowed bigwig files of 3' end RNA-seq data overlapping with poly A site annotations
###############################Work directories

TWD = getwd()
TWD1 = "./bigwigFiles/"
TWD2 = "./CountTables/"

###############################Library input

library(GenomicRanges)
library(GenomicAlignments)
library("BSgenome.Hsapiens.UCSC.hg19")
library(rtracklayer)

seq_info_object1 = SeqinfoForBSGenome(Hsapiens)

#########polyA sites (PAS) downloaded from PolyA_DB3 
##PolyA_DB 3 Catalogs Cleavage and Polyadenylation Sites Identified by Deep Sequencing in Multiple Genomes. Nucleic Acids Research, 2018. Ruijia Wang, Ram Nambiar, Dinghai Zheng, Bin Tian
##Each PAS was extended 20 nt 3’ and 200 nt 5’ from the site of cleavage and those that overlapped on the same strand after extension were combined into a single PAS annotation.
##These extended and combined PAS are saved as a GenomicRanges object in file "reducedHumanPas.Rda"

load("./polyASites/reducedHumanPas.Rda")
posReducedHumanPas = reducedHumanPas[strand(reducedHumanPas) == "+"]
negReducedHumanPas = reducedHumanPas[strand(reducedHumanPas) == "-"]

countTable = data.frame(c(posReducedHumanPas,negReducedHumanPas))

##################################################################################################
########COUNTS TABLE FROM BIGWIGS
#################################################
####Set directory to file containing bigwigs of narrowed reads aligned to hg19 processed from bamfiles
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
  ########################Intersect with reducedHumanPas Granges
  ###################PositiveStrand
  pcov = coverage(p1, weight=p1$score)
  v1 = Views(pcov, posReducedHumanPas)
  s1 = lapply(v1, function(n) sum(n))
  pc = as.vector(unlist(s1))
  ###################NegativeStrand
  ncov = coverage(n1, weight=n1$score)
  v1 = Views(ncov, negReducedHumanPas)
  s1 = lapply(v1, function(n) sum(n))
  nc = as.vector(unlist(s1))
  ##############################################
  countTable$new1 = c(pc,nc)
  names(countTable)[which(names(countTable) == "new1")] = filestem
}

setwd(TWD)
setwd(TWD2)
coldCountTableRaw = countTable
#save(coldCountTableRaw, file="coldCountTableRaw.Rda")
  
