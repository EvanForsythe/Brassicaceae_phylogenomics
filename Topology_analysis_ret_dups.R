#set the working directory
setwd("./R_test_dir/")

#load some packages
install.packages("ape")
library(ape)
library(geiger)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)

#The below package may not be available for every version of R studio
#install.packages("ggtree")
#library(ggtree)

#read ret dup trees 
trees<-read.tree("catfile_retdups_170215_new")

#split trees into subtrees
out_trees_string<-sapply(trees, Split_trees)

Split_trees<-function(tree){
  
  #midpoint root the tree
  root_tree<-midpoint(tree, node.labels="support")
  
  #root_tree$node.label=tree$node.label
  #normalize branch lengths
  #this is required to avoid the "incorrect number of dimensions" error from treeSlice
  root_tree<-compute.brlen(root_tree, 100)
  
  #Slice tree at the root (actually very very close to the root)
  #This outputs a multiphylo of the two sub trees
  sliced_trees<-treeSlice(root_tree, 0.01, trivial=FALSE, prompt=FALSE)
  
  #This if/else loop is to verify that the split worked and yeilded two sub trees
  #if it didn't yeild two subtrees, that means that one of the accessions came out sister to all others.  These trees should be discarded. 
  if (length(sliced_trees) == 2){
    alltips<-root_tree$tip.label
    tips1<-sliced_trees[[1]]$tip.label
    tips2<-sliced_trees[[2]]$tip.label
    
    subtree1<-drop.tip(root_tree, setdiff(alltips, tips2))
    subtree2<-drop.tip(root_tree, setdiff(alltips, tips1))
    
    #print trees
    write.tree(c(subtree1, subtree2), append=TRUE)
  }
  #^^^ end of if loop ^^^
  
  #vvv end of split function vvv
}

#Show trees in which the split did not work.  These should be discarded.
bad_split<-grep("NULL", out_trees_string)

#append subtrees into a multiphylo
out_trees<-read.tree(text=c(out_trees_string))

#count the number of trees
Ntrees<-length(out_trees)


# Make a function for analyzing topologies

TopAnalFunc<-function(root_tree){
  
  ###########################################################################################
  ############                    Root the tree by the Esal tip                  ############
  ###########################################################################################
  
  #normalize the branch lengths because short branches are problematic
  #tree<-compute.brlen(tree, 100)
  
  #tips<-tree$tip.label
  #Esal_tip<-grep("Es_", tips)
  #root_tree<-root(tree, Esal_tip, resolve.root=TRUE, edgelabel=TRUE)
  
  #After rooting, the BS scores are all messed up. Below, I take the scores from the unrooted tree and put them on the rooted tree.  I need to double-check to make sure this is working correctly! 
  
  #root_tree$node.label=tree$node.label
  #root_tree$node.label
  
  
  #store tip names of all relevent species for the droptree
  tips2<-root_tree$tip.label
  Csat_tip<-grep("Cs_", tips2)
  Crub_tip<-grep("Cr_", tips2)
  Cgrand_tip<-grep("Cg_", tips2)
  Athal_tip<-grep("At_", tips2)
  Alyr_tip<-grep("Al_", tips2)
  Bstri_tip<-grep("Bs_", tips2)
  Chir_top<-grep("Ch_", tips2)
  Esal_tip<-grep("Es_", tips2)
  
  
  ###########################################################################################
  ############                    A group monophyly test                         ############
  ###########################################################################################
  
  
  #test if A group seqs are monophyletic and print results
  Agroup_mono<-is.monophyletic(phy=root_tree, c(Athal_tip, Alyr_tip))
  if(Agroup_mono) {Agroupmono = "A group monophyletic"} else {Agroupmono = "A group non-monophyletic"}
  
  
  ###########################################################################################
  ############                    Crub-Cgra group monophyly test                         ############
  ###########################################################################################
  
  
  #test if Crub-Cgrand seqs are monophyletic and print results
  CrubCgrand_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip))
  if(CrubCgrand_mono) {Capgroupmono = "Crub-Cgrand monophyletic"} else {Capgroupmono = "Crub-Cgrand non-monophyletic"}
  
  ###########################################################################################
  ############                     C.sativa paralog monophyly test                         ############
  ###########################################################################################
  
  #test if Csativa seqs are monophyletic and print results
  Csat_mono<-is.monophyletic(phy=root_tree, c(Csat_tip))
  if(Csat_mono) {Csatmono= "Csativa_monophyletic"} else {Csatmono= "Csativa_non_monophyletic"}
  
  
  ###########################################################################################
  ############                    Full C group monophyly test                         ############
  ###########################################################################################
  
  
  #test if all C group seqs are monophyletic and print results
  Cgroup_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip))
  if(Cgroup_mono) {Cgroupmono = "C group monophyletic"} else {Cgroupmono = "C group non-monophyletic"}
  
  ###########################################################################################
  ############                    Full INGROUP group monophyly test                         ############
  ###########################################################################################
  
  
  #test if all In-group seqs are monophyletic and print results
  IN_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip, Bstri_tip))
  if(IN_mono) {INgroup_mono = "IN group monophyletic"} else {INgroup_mono = "C group non-monophyletic"}
  
  
  
  ###########################################################################################
  ############                   Topology analysis of keeper trees               ############
  ###########################################################################################
  
  #Check which clade is monophyletic
  
  BC_clade<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Bstri_tip))
  
  AC_clade<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip))
  
  AB_clade<-is.monophyletic(phy=root_tree, c(Athal_tip, Alyr_tip, Bstri_tip))
  
  #Store the correct topology (STRICT: requires that A, B, and/or C clade are perfectly monogomous)
  if(BC_clade & Cgroup_mono) {final_topology = "BC_topology"} else if(AC_clade & Cgroup_mono & Agroup_mono) {final_topology = "AC_topology"} else if (AB_clade & Agroup_mono) {final_topology = "AB_topology"} else {final_topology = "Other_topology"}
  
  #Store the correct topology (LOOSE: Topology analysis with less-strict monophyly requirements)
  if(BC_clade) {final_topology_loose = "BC_topology"} else if(AC_clade) {final_topology_loose = "AC_topology"} else if (AB_clade) {final_topology_loose = "AB_topology"} else {final_topology_loose = "Other_topology"}
  
  #Store the node representing the MRCA of each potential clade
  #This is the node label at which the crucial BS score resides
  BC_MRCA<-getMRCA(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Bstri_tip))
  AC_MRCA<-getMRCA(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip))
  AB_MRCA<-getMRCA(phy=root_tree, c(Athal_tip, Alyr_tip, Bstri_tip))
  
  #plot.phylo(root_tree, show.node.label=TRUE)
  
  #retrieve the supporting BS score
  if(final_topology_loose == "BC_topology") {BS_score = (root_tree$node.label[(BC_MRCA - length(root_tree$tip.label))])} else if(final_topology_loose == "AC_topology") {BS_score = (root_tree$node.label[(AC_MRCA - length(root_tree$tip.label))])} else if(final_topology_loose == "AB_topology") {BS_score = (root_tree$node.label[(AB_MRCA - length(root_tree$tip.label))])} else {BS_score = "BS_scoreNA"} 
  
  ###Investigating the non-monophyletic topologies
  #What is sister to A. thaliana?
  #Athal_sisters<-c(tips(root_tree, getSisters(root_tree, Athal_tip, mode="number")))
  #Athal_sis_count<-length(Athal_sisters)
  #if(Athal_sis_count==1) 
  #{Athal_sis_paste = Athal_sisters} else if(Athal_sis_count==2) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2])} else if(Athal_sis_count==3) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3])} else if(Athal_sis_count==4) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4])} else if(Athal_sis_count==5) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4], Athal_sisters[5])} else {Athal_sis_paste = "Many_sisters"}
  
  #What is sister to the C.grand - C. rubella clade?
  #CrubCgra_MRCA<-getMRCA(phy=root_tree, c(Crub_tip, Cgrand_tip))
  #CrubCgra_sisters<-c(tips(root_tree, getSisters(root_tree, CrubCgra_MRCA, mode="number")))
  #CrubCgra_sis_count<-length(CrubCgra_sisters)
  #if(CrubCgra_sis_count==1) 
  #{CrubCgra_sis_paste = CrubCgra_sisters} else if(CrubCgra_sis_count==2) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2])} else if(CrubCgra_sis_count==3) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3])} else if(CrubCgra_sis_count==4) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3], CrubCgra_sisters[4])} else if(CrubCgra_sis_count==5) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3], CrubCgra_sisters[4], CrubCgra_sisters[5])} else {CrubCgra_sis_paste = "Many_sisters"}
  
  #plot.phylo(root_tree, show.node.label=TRUE)
  
  
  #find the name of the Athal tip(s)
  Athal_tip2<-grep("At_", tips2, value=TRUE)
  if (length(Athal_tip2) == 0)
  {Athal_tip_name = "No_Athal_tips"} else if (length(Athal_tip2)==1)
  {Athal_tip_name = Athal_tip2} else if (length(Athal_tip2)==2)
  {Athal_tip_name = paste(Athal_tip2[1], Athal_tip2[2])} else if (length(Athal_tip2)==3)
  {Athal_tip_name = paste(Athal_tip2[1], Athal_tip2[2], Athal_tip2[3])}
  
  #find the name of the Crub tip(s)
  Crub_tip2<-grep("Cr_", tips2, value=TRUE)
  if (length(Crub_tip2) == 0)
  {Crub_tip_name = "No_Crub_tips"} else if (length(Crub_tip2)==1)
  {Crub_tip_name = Crub_tip2} else if (length(Crub_tip2)==2)
  {Crub_tip_name = paste(Crub_tip2[1], Crub_tip2[2])} else if (length(Crub_tip2)==3)
  {Crub_tip_name = paste(Crub_tip2[1], Crub_tip2[2], Crub_tip2[3])}
  
  
  return(c(Athal_tip_name, Crub_tip_name, IN_mono, Agroup_mono, CrubCgrand_mono, Csat_mono, Cgroup_mono, final_topology, final_topology_loose, BS_score))
  
}

### END OF FUNCTION ###

#Apply the function to all the trees
output<-lapply(out_trees, TopAnalFunc)


#convert the output from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=Ntrees, byrow=TRUE))
names(output_df) <- c("Athl_seq", "Crub_seq", "In_group_mono", "Agroup_monophyly", "Crub_Cgrand_monophyly", "Csat_monophyly", "C_group_monophyly", "Topology", "Topology_loose", "Bootstrap_Support")


#make a pie chart of the ingroup monophyly
labels<-names(summary(output_df$In_group_mono))
labels<-paste(labels, summary(output_df$In_group_mono))
pie(summary(output_df$In_group_mono), labels=labels, main="Ingroup monophyletic?")


#make a pie chart of topologies (LOOSE)
labels_loose<-names(summary(output_df$Topology_loose))
labels_loose<-paste(labels_loose, summary(output_df$Topology_loose))
pie(summary(output_df$Topology_loose), labels=labels_loose, main="Nuclear Topologies (loose)")


#Subset to include only monophyletic ingroups and non-other topoloies
subsettedIN<-subset(output_df, output_df$In_group_mono=="TRUE")
subsetted_noOther<-subset(subsettedIN, subsettedIN$Topology_loose!="Other_topology")

#piecharts of monophyletic treees only
labels2<-names(summary(subsetted_noOther$Topology_loose))
labels2<-paste(labels2, summary(subsetted_noOther$Topology_loose))
pie(summary(subsetted_noOther$Topology_loose), labels=labels2, main="Topologies")


###Histogram of bootstrap scores
#filter for the ones that aren't NA
boots_only<-subset(subsetted_noOther, subsetted_noOther$Bootstrap_Support != "BS_scoreNA")

#plot  the histogram
full_hist<-hist(as.numeric(as.vector(boots_only$Bootstrap_Support)), breaks=10)

#Split the topologies
boots_AC<-subset(boots_only, boots_only$Topology_loose == "AC_topology")
AC_hist<-hist(as.numeric(as.vector(boots_AC$Bootstrap_Support)), breaks=10)

boots_BC<-subset(boots_only, boots_only$Topology_loose == "BC_topology")
BC_hist<-hist(as.numeric(as.vector(boots_BC$Bootstrap_Support)), breaks=10)

boots_AB<-subset(boots_only, boots_only$Topology_loose == "AB_topology")
AB_hist<-hist(as.numeric(as.vector(boots_AB$Bootstrap_Support)), breaks=10)

plot(full_hist, main="distibution of BS scores", xlab= "BS score", border="black", xlim=c(0,100))
plot(AC_hist, main="distibution of BS scores (AC topology)", xlab= "BS score", border="red", xlim=c(0,100))
plot(BC_hist, main="distibution of BS scores (BC topology)", xlab= "BS score", border="blue", xlim=c(0,100))
plot(AB_hist, main="distibution of BS scores (AB topology)", xlab= "BS score", border="green", xlim=c(0,100))


#export a CSV file
write.csv(output_df, file = "topology_results_retdups_170214.csv")
