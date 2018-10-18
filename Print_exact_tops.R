setwd("./R_test_dir/")

#Load packages
install.packages("ape")
library(ape)
library(geiger)
library(lattice)
library(gtable)
install.packages("gridExtra")
library(gridExtra)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)

#Load trees and store the total number of trees
trees<-read.tree("cattrees_allNuc_170624")
Ntrees<-length(trees)


##################
#### FUNCTION ####
##################
print_newick2<-function(tree){
  
  #root if unrooted
  if (is.rooted(tree)){root_tree=tree} else {
    tips<-tree$tip.label
    Esal_tip<-grep("Es_", tips)
    root_tree<-root(tree, Esal_tip, resolve.root=TRUE, edgelabel=TRUE)
  }
  
  #Assign simple names to taxa
  
  tips2<-root_tree$tip.label
  root_tree$tip.label[grep("Cr_", tips2)]<-"Cr"
  root_tree$tip.label[grep("Cg_", tips2)]<-"Cg"
  root_tree$tip.label[grep("At_", tips2)]<-"At"
  root_tree$tip.label[grep("Al_", tips2)]<-"Al"
  root_tree$tip.label[grep("Bs_", tips2)]<-"Bs"
  root_tree$tip.label[grep("Ch_", tips2)]<-"Ch"
  root_tree$tip.label[grep("Es_", tips2)]<-"Es"
  
  #assign simple names to Csat taxa
  if (length(grep("Cs_", tips2)) == 1){
    root_tree$tip.label[grep("Cs_", tips2)]<-"Cs"
  } else if (length(grep("Cs_", tips2)) == 2){
    root_tree$tip.label[grep("Cs_", tips2)][1]<-"Cs"
    root_tree$tip.label[grep("Cs_", tips2)][2]<-"Cs"
  } else if (length(grep("Cs_", tips2)) == 3){
    root_tree$tip.label[grep("Cs_", tips2)][1]<-"Cs"
    root_tree$tip.label[grep("Cs_", tips2)][2]<-"Cs"
    root_tree$tip.label[grep("Cs_", tips2)][3]<-"Cs"
  }
  
  #Remove branch lengths
  root_tree$edge.length<-NULL
  
  #Remove bootstraps
  root_tree$node.label<-NULL
  
  #ladderize tree
  root_tree<-ladderize(root_tree, right = TRUE)

  
### Rotate branches so that all trees have consistent arrangements
  
  #rotate branches of At/Al
  if (isTRUE(grep("Al,At", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["At", "Al"]) 
  }
  
  #rotate branches of Cr/Cg
  if (isTRUE(grep("Cg,Cr", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["Cr", "Cg"])
  }
  
  #rotate branches of Cr/Cg/Cs
  if (isTRUE(grep("(Cs,Cs),(Cr,Cg)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["Cr", "Cs"])
  }
  
  #rotate branches of Cr/Cg/Cs
  if (isTRUE(grep("(Cr,Cg),(At,Al)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["Cr", "At"])
  }
  
  #rotate branches of Cr/Cg/Cs
  if (isTRUE(grep("(Bs,Cs),(At,Al)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["Bs", "At"])
  }
  
  #rotate branches of Cr/Cg/Cs
  if (isTRUE(grep("(Bs,Cs),(Cr,Cg)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["Bs", "Cr"])
  }
  
  #rotate branches of Cr/Cg/Cs
  if (isTRUE(grep("(Bs,Cs),(Cr,Cg)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["Bs", "Cr"])
  }
 
   #rotate branches of Cr/Cg/Cs
  if (isTRUE(grep("Cs,Cs),Cs),((Cr,Cg),Bs", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-rotateNodes(root_tree, mrca.phylo(root_tree) ["Cs", "Cr"])
  }
  
  #if (grep("Es", root_tree$tip.label)>0) {}
  
  #reroot by Es
  if (isTRUE(grep("(Es,Ch)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-root(root_tree, "Es", resolve.root=TRUE)
  }
  
  #reroot by Es
  if (isTRUE(grep("(Ch,Es)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
    root_tree<-root(root_tree, "Es", resolve.root=TRUE)
  }
  
  #reroot by Es
  if (isTRUE(grep("Es),Ch)", write.tree(root_tree, file = ""), fixed = TRUE) >0)){
   root_tree<-root(root_tree, "Es", resolve.root=TRUE)
   }
  
  #print the simple newick string to standard out
  write.tree(root_tree, file = "")
}
####################
### END FUNCTION ###
####################

newick_out2<-lapply(trees, print_newick2)

newick_out_df2<-data.frame(matrix(unlist(newick_out2), nrow=length(trees), byrow=TRUE))

summary(newick_out_df2)





##BELOW THIS LINE IS OLD CODE. NOT USED TO GENERATE TOPOLOGY TABLE





other_tops<-function(tree){

if (is.rooted(tree)){root_tree=tree} else {
  tips<-tree$tip.label
  Esal_tip<-grep("Es_", tips)
  root_tree<-root(tree, Esal_tip, resolve.root=TRUE, edgelabel=TRUE)
}

#store tip names of all relevent species for the tree
tips2<-root_tree$tip.label
Csat_tip<-grep("Cs_", tips2)
Crub_tip<-grep("Cr_", tips2)
Cgrand_tip<-grep("Cg_", tips2)
Athal_tip<-grep("At_", tips2)
Alyr_tip<-grep("Al_", tips2)
Bstri_tip<-grep("Bs_", tips2)
Chir_tip<-grep("Ch_", tips2)
Esal_tip<-grep("Es_", tips2)

###########################################################################################
############                    A group monophyly test                         ############
###########################################################################################

#test if A group seqs are monophyletic and print results
Agroup_mono<-is.monophyletic(phy=root_tree, c(Athal_tip, Alyr_tip))
#if(Agroup_mono) {Agroupmono = "A group monophyletic"} else {Agroupmono = "A group non-monophyletic"}

###########################################################################################
############                    Crub-Cgra group monophyly test                         ############
###########################################################################################

#test if Crub-Cgrand seqs are monophyletic and print results
CrubCgrand_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip))
#if(CrubCgrand_mono) {Capgroupmono = "Crub-Cgrand monophyletic"} else {Capgroupmono = "Crub-Cgrand non-monophyletic"}

###########################################################################################
############                     C.sativa paralog monophyly test                         ############
###########################################################################################

#test if Csativa seqs are monophyletic and print results
Csat_mono<-is.monophyletic(phy=root_tree, c(Csat_tip))
#if(Csat_mono) {Csatmono= "Csativa_monophyletic"} else {Csatmono= "Csativa_non_monophyletic"}


###########################################################################################
############                    Full C group monophyly test                         ############
###########################################################################################

#test if all C group seqs are monophyletic and print results
Cgroup_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip))
#if(Cgroup_mono) {Cgroupmono = "C group monophyletic"} else {Cgroupmono = "C group non-monophyletic"}

###########################################################################################
############                    Full INGROUP group monophyly test             ############
###########################################################################################

#test if all In-group seqs are monophyletic and print results
IN_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip, Bstri_tip))
#if(IN_mono) {INgroup_mono = "IN group monophyletic"} else {INgroup_mono = "C group non-monophyletic"}



#### Find what's sister to Chir

if(length(Chir_tip)>0 & length(Chir_tip)<2 & 
   length(Athal_tip)>0 & length(Athal_tip)<2 & 
   length(Alyr_tip)>0 & length(Alyr_tip)<2 & 
   length(Crub_tip)>0 & length(Crub_tip)<2 & 
   length(Cgrand_tip)>0 & length(Cgrand_tip)<2 & 
   length(Bstri_tip)>0 & length(Bstri_tip)<2 & 
   length(Csat_tip)>0 &
   getSisters(root_tree, Chir_tip, mode="number")>length(root_tree$tip.label)){


Ch_sis<-extract.clade(root_tree, getSisters(root_tree, Chir_tip, mode="number"))$tip.label

ChAt<-any(grepl(Athal_tip, Ch_sis))
ChAl<-any(grepl(Alyr_tip, Ch_sis))
ChCr<-any(grepl(Crub_tip, Ch_sis))
ChCg<-any(grepl(Cgrand_tip, Ch_sis))
ChCs<-any(grepl(Csat_tip[1], Ch_sis))
ChBs<-any(grepl(Bstri_tip, Ch_sis))
}else{
ChAt<-"NA"
ChAl<-"NA"
ChCr<-"NA"
ChCg<-"NA"
ChCs<-"NA"
ChBs<-"NA"
}

return(c(IN_mono, Agroup_mono, Cgroup_mono, CrubCgrand_mono, Csat_mono, ChAt, ChAl, ChCr, ChCg, ChCs, ChBs))

}

output<-lapply(trees, other_tops)

#convert the output from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=Ntrees, byrow=TRUE))
names(output_df) <- c("In_mono", "A_mono", "C_mono", "CrCg_mono", "Cs_mono", "ChAt", "ChAl", "ChCr", "ChCg", "ChCs", "ChBs")

out_other<-subset(output_df, output_df$In_mono==FALSE | output_df$A_mono==FALSE | output_df$C_mono==FALSE | output_df$CrCg_mono==FALSE | output_df$Cs_mono==FALSE)

