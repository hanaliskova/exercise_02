rm(list=ls())

# path
setwd('C:/Users/hanic/Desktop/PRG/exercise_02')

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install()

# BiocManager::available()

# BiocManager::install("Biostrings")

# libraries
library(Biostrings)
# install.packages('seqinr')
library(seqinr)
library(pwalign)

## TASK 1
# load data - Biostrings, seqinr
seq <- readDNAStringSet("fishes.fna.gz")
dnafile <- system.file("fishes.fna.gz", package = "seqinr")


## TASK 2
length(seq)
width(seq[1])
names(seq)
names(seq[1])
seq1 <- seq[1]
seq1_sequence <- seq[[1]]
seq1_string <- toString(seq[1])
help(XStringSet)


## TASK 3
#install.packages("pwalign")
data("BLOSUM62")

seq2 <- seq[2]
seq2_sequence <- seq[[2]]
seq2_string <- toString(seq[2])

alignseqs <- pairwiseAlignment(seq1_string, seq2_string, substitutionMatrix = 'BLOSUM62', gapOpening = -1, gapExtension = 1)
alignseqs
help(pairwiseAlignment)


## TASK 4
names_list <- c("anna", "jana", "kamil", "norbert", "pavel",
                "petr", "stanislav", "zuzana")
grep("jana", names_list, perl = TRUE)
grep("n+", names_list, perl = TRUE)
grep("n{2}", names_list, perl = TRUE)
grep("^n", names_list, perl = TRUE)
grep("Anna|Jana", names_list, perl = TRUE)
grep("^z.*a$", names_list, perl = TRUE)


## TASK 5
forwardMID <- "ACGAGTGCGT"
reverseMID <- as.character(reverseComplement(DNAString(forwardMID)))

seq_pattern <- paste0("^",forwardMID,".*",reverseMID,"$")
sample_seq <- grep(seq_pattern, seq, perl = TRUE)
num_sample_seq <- length(sample_seq)
num_sample_seq


## TASK 6
MIDs <- read.csv("fishes_MIDs.csv", sep=";")
path <- ('C:/Users/hanic/Desktop/PRG/exercise_02/fishes.fna.gz')
forwardMIDs <- MIDs$FBarcodeSequence
reverseMIDs <- MIDs$RBarcodeSequence
sample_labels <- MIDs$SampleID


Demultiplexer <- function(path, forwardMIDs, reverseMIDs, sample_labels){
  seqs <- readDNAStringSet(path)
  
  # table for report
  report <- data.frame(Sample = sample_labels, Count = 0, stringsAsFactors = FALSE)
  
  # list for saving sequences for each sample
  sample_seqs <- vector("list", length(sample_labels))
  names(sample_seqs) <- sample_labels
  

  for (sample in seq_along(forwardMIDs)){
    print(sample)
    forwardMID <- forwardMIDs[sample]
    reverseMID <- as.character(reverseComplement(DNAString(reverseMIDs[sample])))
    seq_pattern <- paste0("^", forwardMID, ".*", reverseMID, "$")
    
    for (seq in seq_along(seqs)){
      seq_char <- as.character(seqs[[seq]])
      if (grepl(seq_pattern, seq_char, perl = TRUE)){
        # deleting MIDs from beginning and end
        trimmed <- substr(seq_char,
                          nchar(forwardMID) + 1,
                          nchar(seq_char) - nchar(reverseMID))
        
        # add seq into list
        sample_seqs[[sample]] <- c(sample_seqs[[sample]], trimmed)
      }
    }
    
    # actualization of counter in report
    report$Count[sample] <- length(sample_seqs[[sample]])
  }
  
  # save FASTA file for each sample
  for (i in seq_along(sample_labels)){
    if (length(sample_seqs[[i]]) > 0){
      seqset <- DNAStringSet(sample_seqs[[i]])
      names(seqset) <- paste0(sample_labels[i], "_seq", seq_along(seqset))
      writeXStringSet(seqset, filepath = paste0(sample_labels[i], ".fasta"))
    }
  }
  
  # save report
  write.table(report, file = "report.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(report)
}

report <- Demultiplexer(path, forwardMIDs, reverseMIDs, sample_labels)