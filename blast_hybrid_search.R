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

# function for hybrid search
hybrid_search <- function(NTablePath, CTablePath, ident, coverage){
  NTable <- read.table(NTablePath)
  CTable <- read.table(CTablePath)
  filteredDataN <- subset(NTable, NTable[,3]>ident & NTable[,4]>coverage)
  filteredDataC <- subset(CTable, CTable[,3]>ident & CTable[,4]>coverage)
  result <- intersect(filteredDataN$V2,filteredDataC$V2) #besser als File speichern und zurückgeben??
  return(result)
}
#searching for hybrids in BLAST Hit table using 30% identity and 40% coverage
hybrid_search("~/add_volume/uni/project/TetA_N_HitTable.txt","~/add_volume/uni/project/TetC_C_HitTable.txt",30,40)
