#Installing the necessary packages 
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)
# comment out package install commands so you don't accidentally re-run them
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa") 
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
# write.fasta(sequences = seq_AA, names = attr(seq_AA,"Homo_sapiens_6"), file.out = "seq_AA_6.fasta")
# you were almost there! The dollar sign allows you to subset items from a larger set
write.fasta(sequences = seq_AA$Homo_sapiens_6, names = "Homo_sapiens_6", file.out = "seq_AA_6.fasta")




#3. The gene is most likely HBB hemoglobin subunit beta for Homo sapiens. Found using BLAST. 
#LC121775.1 is the accession number. It is a protein coding gene. 

# Homo_Sapien_6 sequence is the most different individual which has 641 bp compared to 642 bp of all the other sequences. 
# It codes for Homo Sapiens mutant hemoglobin beta chain (HBB)
#AY356351.1 is the accession number 

#5. The protein that best matches Homo_Sapien_6 is hemoglobin subunit beta (Homo sapiens)
# KAI2558340.1 is the accession number 

#6. possible diseases of the protein structure include sickle cell and beta thalassemia. 
#Homo_Sapien_6 sequence does not have the diseases 


