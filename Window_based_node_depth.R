#Script for scanning each gene alignment in 200nt windows and calculating node depths from each window
setwd("./R_test_dir/")

#Library some packages
library(ape)
library(dplyr)

#Assign the directory inwhich the alignment files reside
#Alignments available upon request
aln_dir<-"<alignments_directory>"

#Find the names of all alingment files 
all_files<-list.files(path = aln_dir)

#Initialize an empty matrix for results from each alignment
results_mat<-matrix(, nrow = length(all_files), ncol = 7)

#Initialize a list
aln_lengths<-list()

#Initialize a two lists for storing node depths for all BC and AC windows
BC_ND_allwin<-vector("list", length = 100000)
AC_ND_allwin<-vector("list", length = 100000)

#Set the window length
win_length<-200

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

#identify the species names
At_tax<-grep("At_", aln$nam)
Al_tax<-grep("Al_", aln$nam)
Bs_tax<-grep("Bs_", aln$nam)
Cg_tax<-grep("Cg_", aln$nam)
Cr_tax<-grep("Cr_", aln$nam)
Cs_tax<-grep("Cs_", aln$nam)
Ch_tax<-grep("Ch_", aln$nam)
Es_tax<-grep("Es_", aln$nam)

###Create a loop to look through windows###
###########################################

#define the starting sites of each window
starts <- seq(from =1, to =length(aln_mat[1,])-win_length, by = win_length)

#Initialize some temporary lists
n_wins<-length(starts)
win_top<-list()
win_ND<-list()

#Run loop
for(y in 1:n_wins){
  #Get the window
  window<-aln_mat[1:Ntaxa, starts[y]:(starts[y]+(win_length-1))]
  
  #Get the distance
  dist_mat<-dist.dna(as.DNAbin(window), as.matrix = TRUE)
  
  if(!any(is.na(dist_mat)) && !any(is.infinite(dist_mat))){
  
  #Make a tree
  win_tree<-nj(dist_mat)
  
  #Check topology of tree
  if(
    is.monophyletic(win_tree, c(Cr_tax, Cg_tax, Cs_tax, Bs_tax))
  ){
    win_top[y]<-"BC_top"
    win_ND[y]<-mean(c(dist_mat[Cr_tax, Bs_tax], dist_mat[Cg_tax, Bs_tax]))
  }else if(
    is.monophyletic(win_tree, c(Cr_tax, Cg_tax, Cs_tax, At_tax, Al_tax))
  ){
    win_top[y]<-"AC_top"
    win_ND[y]<-mean(c(dist_mat[Cr_tax, At_tax], dist_mat[Cr_tax, Al_tax], dist_mat[Cg_tax, At_tax], dist_mat[Cg_tax, Al_tax]))
  }else{
    win_top[y]<-"Other"
    win_ND[y]<-"NA"
  }
  }#end NA loop
} #end window loop

#Get the average BC and AC node depth for the gene
BC_ND_av<-mean(as.numeric(paste(win_ND[which(win_top=="BC_top")])))
AC_ND_av<-mean(as.numeric(paste(win_ND[which(win_top=="AC_top")])))

#Get the number of windows displaying each topology
BC_win_count<-length(win_ND[which(win_top=="BC_top")])
AC_win_count<-length(win_ND[which(win_top=="AC_top")])

#Store the name of the A.thaliana sequence for the gene
At_seq_name<-aln$nam[grep("At_", aln$nam)]

#Get the length of the alignments
aln_len<-length(aln_mat[1,])

#Store results
results_mat[x,]<-c(At_seq_name, BC_win_count, BC_ND_av, AC_win_count, AC_ND_av, aln_len, n_wins)

#Add to a list of node depths for all windows accross all alignments that display AC or BC
#BC windows
if(length(win_ND[which(win_top=="BC_top")])>0){
BC_ND_allwin[(length(which(BC_ND_allwin!="NULL"))+1):
    ((length(which(BC_ND_allwin!="NULL"))+1)+length(win_ND[which(win_top=="BC_top")])-1)
  ]<-c(as.numeric(paste(win_ND[which(win_top=="BC_top")])))
}#End if loop to see if there are BC windows

#AC windows
if(length(win_ND[which(win_top=="AC_top")])>0){
  AC_ND_allwin[(length(which(AC_ND_allwin!="NULL"))+1):
                 ((length(which(AC_ND_allwin!="NULL"))+1)+length(win_ND[which(win_top=="AC_top")])-1)
               ]<-c(as.numeric(paste(win_ND[which(win_top=="AC_top")])))
}#End if loop to see if there are AC windows

}#End of alignment loop

#Convert to dataframe and add column names
results_df<-as.data.frame(results_mat)
names(results_df)<-c("At_taxon", "BC_win_count", "BC_ND_av", "AC_win_count", "AC_ND_av", "Aln_len", "n_wins")

#Subset by the majority topology for each gene
BC_top_genes_df<-subset(results_df, as.numeric(paste(results_df$BC_win_count)) > as.numeric(paste(results_df$AC_win_count)))
AC_top_genes_df<-subset(results_df, as.numeric(paste(results_df$BC_win_count)) < as.numeric(paste(results_df$AC_win_count)))

#Plot boxblot to compare BC and AC node depth
boxplot(as.numeric(paste(BC_top_genes_df$BC_ND_av)), as.numeric(paste(AC_top_genes_df$AC_ND_av)), outline = FALSE)

#Perform wilcoxon test of distribution differences
wilcox.test(as.numeric(paste(BC_top_genes_df$BC_ND_av)), as.numeric(paste(AC_top_genes_df$AC_ND_av)))

#Read topology and gene coordinate data
Tops<-read.csv(file = "topology_results170509.csv")

#Combine Dstat results and topology results
combined_df<-right_join(Tops, results_df, by = "At_taxon")

#Subset by topology
BC_ml_top_df<-subset(combined_df, combined_df$Topology == "BC_topology")
AC_ml_top_df<-subset(combined_df, combined_df$Topology == "AC_topology")

#Make a boxplot using topology assignments from ML analysis
boxplot(as.numeric(paste(BC_ml_top_df$BC_ND_av)), as.numeric(paste(AC_ml_top_df$AC_ND_av)), outline = FALSE)

#Perform wilcoxon test of distribution differences
wilcox.test(as.numeric(paste(BC_ml_top_df$BC_ND_av)), as.numeric(paste(AC_ml_top_df$AC_ND_av)))

#Ask how many BC windows are in AC trees
contam_ratio_BCtrees<-as.numeric(paste(BC_ml_top_df$AC_win_count))/as.numeric(paste(BC_ml_top_df$n_wins))
contam_ratio_ACtrees<-as.numeric(paste(AC_ml_top_df$BC_win_count))/as.numeric(paste(AC_ml_top_df$n_wins))

#Ask how many BC windows are in AC trees
not_contam_ratio_BCtrees<-as.numeric(paste(BC_ml_top_df$BC_win_count))/as.numeric(paste(BC_ml_top_df$n_wins))
not_contam_ratio_ACtrees<-as.numeric(paste(AC_ml_top_df$AC_win_count))/as.numeric(paste(AC_ml_top_df$n_wins))

#Histogram of discordant windows
hist(contam_ratio_BCtrees)
hist(contam_ratio_ACtrees)

#Histogram of concordant windows
hist(not_contam_ratio_BCtrees)
hist(not_contam_ratio_ACtrees)



