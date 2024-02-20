#Installing the necessary packages 
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msa") 
`force = TRUE`
library(msa)

#set-up a working directory which will keep all my files in the same place and easy located on my computer

setwd("/Users/Sydne/OneDrive/Desktop/Github_Bio/Bioinformatic/Midterm")

#command readDNAStringset creates a string for the fasta file. folder and name of the fasta file is provided. Transfers the data from the fasta file to a variable named 'seq'

seq <- readDNAStringSet("Seq/sequences.fasta")

#command msaClustalW aligns the DNA sequences of the fasta file using a program? ClustalW. Has a new variable name of 'seqaln'
#command nchar provides the amount of bp nucleotides in the sequence 
# alFreq provides a matrix comparing amount of A,T,C,G nucleotides to each individual sequence
seqaln <- msaClustalW(seq)
nchar(seqaln)
alFreq <- alphabetFrequency(seqaln)

#Translates the DNA fasta file to an Amino acid file. New variable of 'seq_AA'
#print 'seq_AA' will paste the the amino acids of all 20 sequences into the console 
#Command seq_AA[["Homo_sapiens_6"]] just views the amino acids from the 6th sequence 
seq_AA <- Biostrings::translate(seq)
print(seq_AA)
seq_AA[["Homo_sapiens_6"]]

#Creates a new Fasta file with the aligned DNA of all 20 sequences
SeqAln_phyDat <- msaConvert(seqaln, type="phangorn::phyDat")
write.phyDat(SeqAln_phyDat, "Seq/Homo_Sapien6_alignment.fasta", format = "fasta")

#Creates a new Fasta file with the translated Amino acid only of the 6th sequence 
write.fasta(sequences = seq_AA, names = attr(seq_AA,"Homo_sapiens_6"), file.out = "seq_AA_6.fasta")


