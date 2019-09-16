### FIRST STEPS ###

## Installtion Miniconda
## Shell 
## wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
## sh Miniconda3-latest-Linux-x86_64.sh
## 
## conda install -c bioconda bioconductor-biostrings 

setwd("/home/alicia/add_volume/uni/project")

#source("https://bioconductor.org/biocLite.R")
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("Biostrings")
#.libPaths()
#install.packages("reutils", lib = "/home/alicia/R/x86_64-pc-linux-gnu-library/3.5")
#installed.packages()
#install.packages("XML")
#biocLite("Biostrings")
#biocLite("msa")

library(reutils)
library(Biostrings)
library(msa)

##### first steps #######

# Tetracycline class C [Escherichia coli]
# TetA tetracycline resistance protein, class C (plasmid) [Escherichia coli]
# GenBank: QAX88993.1
tmp = tempfile()
idTest <- c(classA="ACQ42041",classB="AOM73298")
res <- efetch(idTest, db="protein", retmode="text", rettype="fasta", outfile=tmp)
aaTest <- readAAStringSet(tmp)
writeXStringSet(aaTest,"NCBI_Tetra_test.fasta")

## EInlesen FASTA Datei
aaSetA <- readAAStringSet("tetA_classA_EC.fasta")
aaSetB <- readAAStringSet("tetA_classB_EC.fasta")
aaSet <- c(aaSetA,aaSetB)

#Alignment
alignAB <- msa(aaTest)

# Challange: Berechnen der relativen Häufigkeiten konservierter AS
cons <- msaConsensusSequence(alignAB)
consVector <- strsplit(cons, split="")
absolut <- 0
for (x in consVector[[1]]){
  for (l in c(LETTERS)){ #läuft über gesamtes Alphabet (nicht nur über AS), ist aber auch egal, oder?
    if (x==l){
      absolut =  absolut+1
    }
  }
}
(absolut/nchar(cons))*100