### This is an R script for cleaning whole proteome CDS files
# This step is necessary before using the AA fasta files to run Orthofinder.
# The script accomplished the following tasks:
# If the CDS files contains multiple isoforms per gene, this scirpt "cleans" CDS files by creating a new CDS file that only contains the longest isoform for each genes
# To ensure perfect agreement between CDS files and AA files (which will save a lot of heartache downstream), this script translates CDS in AA and outputs an AA file.

#Set the working directory to a directory containing the "uncleaned" CDS files.
setwd("./R_test_dir/")

# Load and library the package, seqinr
install.packages("seqinr")
library(seqinr)
#install.packages("biostrings")
#library(biostrings)

# Load input file
# IMPORTANT: Be sure to only uncomment the needed file and to make sure it is consistent with the read.fasta command below.
#fa <- read.fasta("Cs_genes_v2.cds.fa")
#fa <- read.fasta("Alyr_CDS.fa")
#fa <- read.fasta("Athaliana_167_TAIR10.cds.fa")
#fa <- read.fasta("Bstricta_278_v1.2.cds.fa")
#fa <- read.fasta("Cgrandiflora_266_v1.1.cds.fa")
fa <- read.fasta("carhr38.cds.fa")

# NOTE: C. rubella and E. salsugienum CDS files only have 1 isoform so I didn't need to do this first part on them. 


# Get sequences and gene names from fa object
genes <- sapply(strsplit(names(fa), "\\."), function(v) {return(v[1])})

#Convert the sequence, which is stored as a vector of characters, to a single string
sequences <- sapply(fa, c2s)

#Extract longest transcript (or first transcript if they are the same length)
filtered_seq <- tapply(sequences, genes, function(v) {if (length(v)==2 & nchar(v[1])!=nchar(v[2])) 
  {return(v[which(nchar(v)==max(nchar(v)))])} else if (length(v)==3 & nchar(v[1])!=nchar(v[2]) & nchar(v[2])!=nchar(v[3]) & nchar(v[1])!=nchar(v[3])) 
    {return(v[which(nchar(v)==max(nchar(v)))])} else {return(v[1])}})

###Ignore this command
#Extracting longest transcript (or first transcript if they are the same length)
#filtered_seq <- tapply(sequences, genes, function(v) {if (length(v)==2 & nchar(v[1])!=nchar(v[2])) {return(v[which(nchar(v)==max(nchar(v)))])} else {return(v[1])}})


#Make a function for creating an object for write.fasta
s2c_funct<-function(seq) {s2c(seq[[1]])}

# Use the function to create an object suitable for the write.fasta function
# This might throw an error but it should still work
obj <- tapply(filtered_seq, 1:length(filtered_seq), s2c_funct)

# Writing output file
write.fasta(obj, names(filtered_seq), file="Ch_CDS_1iso.fa")

### read new trimmed (cleaned) fasta
# Load the cleaned file
fa2 <- read.fasta("Ch_CDS_1iso.fa")

# Get sequences and gene names from fa object
genes2 <- sapply(strsplit(names(fa2), "\\."), function(v) {return(v[1])})

#Again, Convert the sequence to a single string
sequences2 <- sapply(fa2, c2s)

#Create function for translating
#This will only attempt to translate sequences that are not "na"
translate_func<-function(seq2) {if (seq2[1]!="na") {translate(s2c(seq2[1]))}}

#Run the function
prot_seqs<-tapply(sequences2, 1:length(sequences2), translate_func)

#Write the protein sequences to a fasta file! 
write.fasta(prot_seqs, names(sequences2), file="Ch_PROT_1iso.fa")
