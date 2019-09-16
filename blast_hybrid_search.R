### BLAST AND HYBRID SEARCH ###

#create fasta file of db seq.
# got to parent directory
#BASH: cat $(find -name "*.faa") > all_seq.faa
#
#makeblastdb
#delete duplicated seq.IDs
faaFile <- "~/add_volume/uni/db/all_seq.faa"
aaSetDB <- readAAStringSet(faaFile)
system2(command = "makeblastdb",
        args = c("-in", faaFile,
                 "-title blastDB",
                 "-out tetraDB"))

# https://blast.ncbi.nlm.nih.gov/Blast.cgi
## BLAST Bash
#script aller Proteine einer Klasse
#"6 qseqid sseqid pident qcovs qlen length mismatch gapopen evalue bitscore"

#Hit table einlesem
resultsDF <- read.table("~/add_volume/uni/project/TetA_N_HitTable.txt")
resultD <- read.table("~/add_volume/uni/project/TetC_C_HitTable.txt")
table(resultsDF[,1])
table(resultsDF[,2])
dim(resultsDF)
#in table filtern (identity + coverage?!)
filteredDF <- subset(resultsDF,resultsDF[,3] > 30 & resultsDF[,4] > 40)
filteredD <- subset(resultD,resultD[,3] > 30 & resultD[,4] > 40)
dim(filteredD)
#Listen vergleichen
intersect(filteredDF$V2,filteredD$V2)

#To-Do:
#funktion schreiben, die Listen vergleicht, Suchstränge als zu übergebende Parameter
#blast-ausgabe Spaltennamen
#Namen der Plots (TetA, TetC usw.)