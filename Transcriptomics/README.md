## TRANSCRIPTOMICS ANALYSIS
### BamToBigwigs
RNA samples from nuclear and cytoplasmic fractions of AC16 and U2OS cells exposed to different temperature conditions were generated into libraries using the QuantSeq 3’ mRNA-Seq Library Prep Kit for Ion Torrent (Lexogen). Libraries were sequenced using the Ion Proton Sequencer and aligned to the hg19 human genome build or a combined dm6 (*D. melanogaster*) / hg19 genome build (for *D. melanogaster* cell spiked in samples) using the TMAP aligner with alignment settings (-tmap mapall stage1 map4) to generate bam files. Bam file reads were narrowed to their 3' nucleotide and converted to bigwig files using R scripts `bamfile_to_narrowed_bigwigs.R` and `bamfile_spikedin_to_narrowed_bigwigs.R`.

### Count Table Generation and Processing
Human polyA site (PAS) annotations were obtained from the Tian lab PolyA_DB 3 (http://exon.umdnj.edu/polya_db/v3) (Wang et al, 2018). Each PAS was extended 20 nt 3’ and 200 nt 5’ from the site of cleavage and those that overlapped on the same strand after extension were combined into a single PAS annotation: see GRanges object `reducedHumanPas.Rda`. 
D. melanogaster PASs were obtained from the Tian lab (Liu et al, 2017) and were extended, in the same way as for human PASs: see GRanges object `reducedFlyPas.Rda`.
Counts per PAS for each RNA-seq sample were obtained using bigwig files (from GEO (http://https://www.ncbi.nlm.nih.gov/geo/) under accession code GSE137003) and R scripts `coldCountTableGeneration.R` and `countTableGenerationFly.R` to give `coldCountTableRaw.Rda` and `FlyCountTableRaw1.Rda`, respectively.

Raw count tables were processed to give:
1) Counts per Genes in RefSeq data base (counts from all PASs associated with each gene are summed) `coldRefSeqGeneCountTable.Rda`.
2) Reads per million PAS-associated reads (RPM)-normalized counts per PAS `coldPasForRefSeqGeneRPMCountTable.Rda`.
3) Counts in reads per million per Genes in RefSeq data base (RPM-normalized counts from all PASs associated with each gene are summed) `coldRefSeqGeneRPMCountTable.Rda`.

### Differential Gene Expression
Differential gene expression between selected temperature conditions (e.g. 18d24h37d2h = 18°C for 24h, then back to 37°C for 2h) was carried out using the DESeq2 package in R (Love et al., 2014) using R script `DESeqDifferentialGeneExpressionRNAseq.R` to give the csv files stored in DESeq/DifferentialGeneExpressionGeneLists/. C or N refers to cytoplasmic or nuclear samples.

### Additional R scripts
Remaining scripts generate charts and statistical analyses presented in Fig. 2-3, Fig. EV2-3 and the source data file of the paper. 
