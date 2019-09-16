### SEQUENCE DOWNLOAD AND FIRST ANALYSIS ###
library(reutils)
library(Biostrings)
library(msa)

## Marylin roberts
# https://faculty.washington.edu/marilynr/
# tet(A), tet(B), tet(C), tet(D), tet(E)
# tet(G), tet(H), tet(J), tet(V), tet(Y) 
# tet(Z), tet(30), tet(31), tet(33), tet(57), tet(59)

## Download bacterial proteins (diverse species)
## Load all tetracycline proteins
# ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete
# wget ftp://ftp.ncbi.nih.gov/refseq/release/bacteria/*.faa.gz

#TetR -> repressor
# Class A # CAD11598.1 [Escherichia coli] 
# Class B # AAF68935.1 [Salmonella typhi]
# Class C # AAG43219.1 [Salmonella typhimurium]
# Class D # CAA46707.1 [Salmonella ordonez]
# Class E # ABO92307.1 [Aeromonas salmonicida]
# Class G # AAC98495.1 [Salmonella enterica]
# Class H # CAC08219.1 [Pasteurella aerogenes]
# Class J # AAD12754.1 [Proteus mirabilis]
# Class V # 
# Class Y # ALV74598.1
# Class Z # AAD25064.1 (plasmid) [Corynebacterium glutamicum]
# Class 30 # AAD09859.1 [Agrobacterium fabrum str. C58]
# Class 31 # CAC80726.1 (plasmid) [Aeromonas salmonicida subsp. salmonicida]
# Class 33 # CAD12226.1 (plasmid) [Corynebacterium glutamicum]
# Class 57 # AJO67550.1
# Class 59 # AMP42441.1[uncultured bacterium IN-13]

idRep <- c("CAD11598","AAF68935","AAG43219","CAA46707","ABO92307","AAC98495","CAC08219","AAD12754","ALV74598","AAD25064","AAD09859","CAC80726","CAD12226","AJO67550","AMP42441")

tmp = tempfile()
rep <- efetch(idRep, db="protein", retmode="text", rettype="fasta", outfile=tmp)
# AAString
aaSetRep <- readAAStringSet(tmp)
names(aaSetRep)
??gsub
names(aaSetRep) <- gsub(" (t|r).*", "", names(aaSetRep), ignore.case = TRUE)
#create C- and N-terminal subseq. 
aaSetRepN <- subseq(aaSetRep, 1, 50)
aaSetRepC <- subseq(aaSetRep, width(aaSetRep)-49,)
#write fasta
writeXStringSet(aaSetRepC,"NCBI_TC_Repressor_C-terminal.fasta")
writeXStringSet(aaSetRepN,"NCBI_TC_Repressor_N-terminal.fasta")
writeXStringSet(aaSetRep,"NCBI_TC_Repressor.fasta")

alignRep <- msa(aaSetRep)
str(alignRep)
names(alignRep)

## How to compute the number of matches / and mismatches from align result ?
msaPrettyPrint(alignRep, output="pdf",file="MSA_TC_Repressor.pdf", showNames="left", showLogo="none",consensusColor="ColdHot", showLegend=FALSE,askForOverwrite=FALSE)
#nmatch(alignRes)

#install.packages("seqinr")
library(seqinr)
##distance matrix
alignRep <- msaConvert(alignRep, type = "seqinr::alignment")
distRep <- dist.alignment(alignRep, "similarity")
as.matrix(distRep)
##phylogenetic tree
#install.packages("ape")
library(ape)

#agglomerative Clustering
aggloTreeRep <- hclust(distRep)
#Label renaming???????
#aggloTreeRep$labels[aggloTreeRep$labels=="CAD11598.1 repressor (plasmid) [Escherichia coli]"] <- "Class A"
plot(aggloTreeRep, main = "Phylogenetic Tree of Tetracycline Repressor Proteins") #schriftgröße noch anpassen

#NJ
njTreeRep <- nj(distRep)
plot(njTreeRep,hang=-1, main = "Phylogenetic Tree of Tetracycline Repressor Proteins")

## heatmap
#install.packages("pheatmap")
#install.packages("RColorBrewer")
library(pheatmap)
library(RColorBrewer)
#par(mar=c(7,8,8,8),oma=c(7,8,8,8)) <-Wofür?????
mycols <- colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
#label names????????
pheatmap(distRep,main= "Tretracycline repressor genes",color=mycols,show_rownames=TRUE,show_colnames=TRUE,fontsize_row=14,fontsize_col=14, filename = "heatmap_tet_repressor.pdf")
