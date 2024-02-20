library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa") 
a
library(msa)

setwd("/Users/Sydne/OneDrive/Desktop/Github_Bio/Bioinformatic/Midterm")
seqs <- system.file("/Users/Sydne/OneDrive/Desktop/Github_Bio/Bioinformatic/Midterm/sequences.fasta")

sequences <- c("sequences.fasta")

seqaln <- msaClustalW(sequences)

readDNAStringSet(dna,sequences.fasta)


fas <- ("/Users/Sydne/OneDrive/Desktop/Github_Bio/Bioinformatic/Midterm")
fas <- system.file("extdata","sequences.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

AA <- AlignTranslation(dna, type="AAStringSet") # align the translation
BrowseSeqs(AA, highlight=1) # view the alignment
DNA <- AlignSeqs(dna)



