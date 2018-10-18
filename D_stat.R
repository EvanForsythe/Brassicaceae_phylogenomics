#Script for calculating D statistics from gene alignments

setwd("./R_test_dir/")
library(ape)
library(dplyr)

#The directory inwhich the alignment file reside
#Alignments available upon request
aln_dir<-"<alignment_directory>"

#Find the names of all alingment files 
all_files<-list.files(path = aln_dir)

#Initialize the site pattern tickers
AABB<-0
ABBA<-0
BABA<-0
F_denom<-0

#Initialize an empty matrix for results from each alignment
results_mat<-matrix(, nrow = length(all_files), ncol = 8)

#Begin loop to cycle through all alignments
for(x in 1:length(all_files)){
 
#Read the alignment
dna_bin<-read.dna(file = paste(aln_dir,all_files[x], sep = ""), format = "sequential")
aln<-as.alignment(dna_bin)

#initialize a matrix
aln_mat <- matrix(, length(aln$nam), nchar(aln$seq[[1]]))

#Populate the matrix
for(i in 1:length(aln$nam)){
  aln_mat[i, ] <- toupper(unlist(strsplit(aln$seq[[i]], ""))) #toupper is to make sure sequence is uppercase
}

#count number of taxa
Ntaxa<-length(aln$nam)

At_tax<-grep("At_", aln$nam)
Al_tax<-grep("Al_", aln$nam)
Bs_tax<-grep("Bs_", aln$nam)
Cg_tax<-grep("Cg_", aln$nam)
Cr_tax<-grep("Cr_", aln$nam)
Ch_tax<-grep("Ch_", aln$nam)
Es_tax<-grep("Es_", aln$nam)

Cs_tax_temp<-grep("Cs_", aln$nam)
if(length(Cs_tax_temp)==1){
  Cs_tax1<-grep("Cs_", aln$nam)
}else if(length(Cs_tax_temp)==2){
  Cs_tax1<-grep("Cs_", aln$nam)[1]
  Cs_tax2<-grep("Cs_", aln$nam)[2]
}else if(length(Cs_tax_temp)==3){
  Cs_tax1<-grep("Cs_", aln$nam)[1]
  Cs_tax2<-grep("Cs_", aln$nam)[2]
  Cs_tax3<-grep("Cs_", aln$nam)[3]
}

#Initialize the site pattern tickers
AABB<-0
ABBA<-0
BABA<-0
F_denom<-0

#Loop through all sites in the alignment
for(i in 1:ncol(aln_mat)){

#D-STAT
  #Find the informative sites in the alignment
  if(length(unique(aln_mat[, i])) == 2 && #Is site biallelic?
    !any(aln_mat[, i]=="-") && #Does site contain gaps?
    aln_mat[Es_tax, i]==aln_mat[Ch_tax, i] && #Do two outgroups have the same allele?
    aln_mat[At_tax, i]==aln_mat[Al_tax, i] && #Is the A clade monoallelic?
    aln_mat[Cr_tax, i]==aln_mat[Cg_tax, i] && #Is the C clade monoallelic?
    length(unique(c(aln_mat[At_tax, i], aln_mat[Cr_tax, i], aln_mat[Bs_tax, i]))) == 2){ #Is there polymorphism between ingroup clades?
      #Find the site pattern
      if(aln_mat[Es_tax, i] == aln_mat[At_tax, i] && aln_mat[Bs_tax, i] == aln_mat[Cr_tax, i]){
        AABB<-AABB+1
      }else if(aln_mat[Es_tax, i] == aln_mat[Bs_tax, i] && aln_mat[At_tax, i] == aln_mat[Cr_tax, i]){
        ABBA<-ABBA+1
      }else if(aln_mat[Es_tax, i] == aln_mat[Cr_tax, i] && aln_mat[At_tax, i] == aln_mat[Bs_tax, i]){
        BABA<-BABA+1
    }
}#1st if 
  
  #F_stat (See Zheng and Janke, 2018, BMC bioinformatics)
  #Find the informative sites in the alignment
  if(length(unique(aln_mat[, i])) == 2 && #Is site biallelic?
     !any(aln_mat[, i]=="-") && #Does site contain gaps?
     aln_mat[Es_tax, i]==aln_mat[Ch_tax, i] && #Do two outgroups have the same allele?
     aln_mat[At_tax, i]==aln_mat[Al_tax, i] && #Is the A clade monoallelic?
     #aln_mat[Cr_tax, i]==aln_mat[Cg_tax, i] && #Is the C clade monoallelic?
     length(unique(c(aln_mat[At_tax, i], aln_mat[Bs_tax, i]))) == 2){ #Is there polymorphism between clades B and C?
    #Find the site pattern
    if(aln_mat[Es_tax, i] == aln_mat[Bs_tax, i]){
      F_denom<-F_denom+1
    }
  }#1st if 
  
  }#for

results_mat[x,]<-c(all_files[x], aln$nam[At_tax], aln$nam[Cr_tax], nchar(aln$seq[[1]]), AABB, ABBA, BABA, F_denom)

}#End Alignment cycling loop 

#Clean up results
colnames(results_mat)<-c("File", "At_taxon", "Cr_taxon", "Aln_length", "AABB", "ABBA", "BABA", "F_denom")
results_df<-as.data.frame(results_mat)

#Read topology and gene coordinate data
Tops<-read.csv(file = "Topology_results.csv")

#Combine Dstat results and topology results
combined_df<-right_join(Tops, results_df, by = "At_taxon")

###SUBSETTING TO LOOK AT CSS OR INDIVIDUAL CHROMS 

##TOGGLE BELOW##
#Subset for conservatively single-copy genes (toggle this)
#combined_df<-subset(combined_df, single_copy_status=="all_three")

#Subset by chromosome
chrom_keep<-8
combined_df<-subset(combined_df, Crub_chrom==chrom_keep)

#Subset by topology
#combined_df<-subset(combined_df, Topology=="AB_topology")

#Sum the total counts of site patterns
AABB_total_ob<-sum(as.numeric(paste(combined_df$AABB)))
ABBA_total_ob<-sum(as.numeric(paste(combined_df$ABBA)))
BABA_total_ob<-sum(as.numeric(paste(combined_df$BABA)))
F_denom_total<-sum(as.numeric(paste(combined_df$F_denom)))

#Calculate D (observed full dataset)
D<-(ABBA_total_ob-BABA_total_ob)/(ABBA_total_ob+BABA_total_ob)
F<-(ABBA_total_ob-BABA_total_ob)/F_denom_total

#Create loop for resampling
nreps<-10000 #number of bootstrap replicates
D_sample<-numeric()
F_sample<-numeric()

#Loop
for(y in 1:nreps){

resample_df<-results_df[sample(1:nrow(combined_df), replace = TRUE),]

#Sum the total counts of site patterns
AABB_total<-sum(as.numeric(paste(resample_df$AABB)))
ABBA_total<-sum(as.numeric(paste(resample_df$ABBA)))
BABA_total<-sum(as.numeric(paste(resample_df$BABA)))
F_denom_total<-sum(as.numeric(paste(resample_df$F_denom)))

#Calculate D (observed full dataset)
D_sample[y]<-(ABBA_total-BABA_total)/(ABBA_total+BABA_total)
F_sample[y]<-(ABBA_total-BABA_total)/F_denom_total
}

#Plot distribution of D
plot(density(D_sample), xlim=c(-0.05, 0.15))
abline(v=0)
abline(v=D)

#Z-test for D
a<-0
s<-sd(D_sample)
n<-length(D_sample)
xbar<-mean(D_sample)
z_d<-(xbar-a)/s
p_d<-2*pnorm(-abs(z_d))


#Plot distribution of F
plot(density(F_sample), xlim=c(-0.02, 0.04))
abline(v=0)
abline(v=F)

#Z-test for F
a<-0
s<-sd(F_sample)
n<-length(F_sample)
xbar<-mean(F_sample)
z_f<-(xbar-a)/s
p_f<-2*pnorm(-abs(z_f))

fortable<-c(chrom_keep, AABB_total_ob, ABBA_total_ob, BABA_total_ob, D, z_d, paste(p_d), F, z_f, paste(p_f))



