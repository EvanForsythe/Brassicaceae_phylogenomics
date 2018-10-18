#Set working directory
setwd("./R_test_dir/")

#load some packages
install.packages("ape")
library(ape)
library(geiger)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)
library(plyr)
install.packages("ggplot2")
library(ggplot2)

#Import concatinated trees
trees<-read.tree("cat_trees_nuc170213.txt")

#Assign one tree for testing
tree<-trees[[500]]

#take a look at the tree
plot.phylo(tree)

#Make function
lambda_set<-100
lambda<-lambda_set
chromo_nodedepth<-function(tree){
#Root unrooted trees
if (is.rooted(tree)){root_tree=tree} else {
  tips<-tree$tip.label
  Esal_tip<-grep("Es_", tips)
  root_tree<-root(tree, Esal_tip, resolve.root=TRUE, edgelabel=TRUE)
}

#take a look at the new rooted tree
#plot.phylo(root_tree)

#Turn tree into chronogram
#Note that lambda=0 here. I may need to explore this parameter
chrono_tree<-chronopl(root_tree, lambda = lambda_set)

#Take a look at the tree
#plot.phylo(chrono_tree)

#store tip names of all relevent species for the tree
tips2<-chrono_tree$tip.label
Csat_tip<-grep("Cs_", tips2)
Crub_tip<-grep("Cr_", tips2)
Cgrand_tip<-grep("Cg_", tips2)
Athal_tip<-grep("At_", tips2)
Alyr_tip<-grep("Al_", tips2)
Bstri_tip<-grep("Bs_", tips2)
Chir_top<-grep("Ch_", tips2)
Esal_tip<-grep("Es_", tips2)

#store Athal tip name as text
Athal_tip_name<-tree$tip.label[Athal_tip]

#Get the node depths for the the chrono tree
all_depths<-node.depth.edgelength(chrono_tree)

#Retun the Most Recent Common Ancestor node for each clade
A_MRCA<-getMRCA(phy=chrono_tree, c(Athal_tip, Alyr_tip))
C_MRCA<-getMRCA(phy=chrono_tree, c(Crub_tip, Cgrand_tip, Csat_tip))
AC_MRCA<-getMRCA(phy=chrono_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip))
BC_MRCA<-getMRCA(phy=chrono_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Bstri_tip))
AB_MRCA<-getMRCA(phy=chrono_tree, c(Athal_tip, Alyr_tip, Bstri_tip))

#Return the node depth (distance from the tips)
A_Ndepth<-1-all_depths[A_MRCA]
C_Ndepth<-1-all_depths[C_MRCA]
AC_Ndepth<-1-all_depths[AC_MRCA]
BC_Ndepth<-1-all_depths[BC_MRCA]
AB_Ndepth<-1-all_depths[AB_MRCA]

#Print node depths 
return(c(Athal_tip_name, as.numeric(A_Ndepth), as.numeric(C_Ndepth), as.numeric(AC_Ndepth), as.numeric(BC_Ndepth), as.numeric(AB_Ndepth)))
}

#############################
###### END FUNCTION  ########
#############################

#run the script
output<-lapply(trees, chromo_nodedepth) 

#convert the output from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=length(trees), byrow=TRUE), stringsAsFactors = FALSE)
names(output_df) <- c("Name", "A_depth", "C_depth", "AC_depth", "BC_depth", "AB_depth")
#output_df$Name<-as.factor(output_df$Name)

#read in a table with topology results
tops<-read.csv("topology_results170509.csv", header = TRUE)

#Join the two data frames by the Athal seq name
Node_df<-join(output_df, tops, type = "inner")

#Subset by topology
BC_sub<-subset(Node_df, Node_df$Topology_loose=="BC_topology")
AC_sub<-subset(Node_df, Node_df$Topology_loose=="AC_topology")
AB_sub<-subset(Node_df, Node_df$Topology_loose=="AB_topology")

#Subset chloro and mito
CP_sub<-subset(Node_df, Node_df$Genome=="cp")
MT_sub<-subset(Node_df, Node_df$Genome=="mt")

#Stack the two histograms
layout(matrix(c(1,2), byrow=TRUE))
 par(mar=c(2,4,1,1))
#Density plots for T2
plot(density(as.numeric(BC_sub$BC_depth)), col="orange", ylim=c(1,12), xlim=c(0,1), main =paste("lambda=", lambda), xlab = "")
lines(density(as.numeric(AC_sub$AC_depth)), col="green")
lines(density(as.numeric(AB_sub$AB_depth)), col="purple")
par(mar=c(2,4,1,1))
#Density plots for T1
plot(density(as.numeric(BC_sub$AC_depth)), col="orange", ylim=c(1,12), xlim=c(0,1), main = "", xlab = "node depth")
lines(density(as.numeric(AC_sub$BC_depth)), col="green")
lines(density(as.numeric(AB_sub$AC_depth)), col="purple")

#set par back to normal
par(mfrow=c(1,1))

###Use the Shapiro-Wilk test to see if distriubtions are normal
#This test is only valid for sample sizes of 5000 or less.
#Further, this test may not be valid ever. See the link below:
#https://stats.stackexchange.com/questions/2492/is-normality-testing-essentially-useless
shapiro.test(sample(as.numeric(BC_sub$BC_depth), 5000))
shapiro.test(as.numeric(AC_sub$AC_depth))
shapiro.test(as.numeric(AB_sub$AB_depth))



TOPOLOGY<-c("1_BC_T1", "2_AC_T1", "3_AB_T1", "4_BC_T2", "5_AC_T2", "6_AB_T2")
  
MEAN<-c(
        mean(as.numeric(BC_sub$AC_depth)),
        mean(as.numeric(AC_sub$BC_depth)),
        mean(as.numeric(AB_sub$AC_depth)),
        mean(as.numeric(BC_sub$BC_depth)),
        mean(as.numeric(AC_sub$AC_depth)),
        mean(as.numeric(AB_sub$AB_depth))
)

YMAX<-c(
  mean(as.numeric(BC_sub$AC_depth))+sd(as.numeric(BC_sub$AC_depth)),
  mean(as.numeric(AC_sub$BC_depth))+sd(as.numeric(AC_sub$BC_depth)),
  mean(as.numeric(AB_sub$AC_depth))+sd(as.numeric(AB_sub$AC_depth)),
  mean(as.numeric(BC_sub$BC_depth))+sd(as.numeric(BC_sub$BC_depth)),
  mean(as.numeric(AC_sub$AC_depth))+sd(as.numeric(AC_sub$AC_depth)),
  mean(as.numeric(AB_sub$AB_depth))+sd(as.numeric(AB_sub$AB_depth))
)

YMIN<-c(
  mean(as.numeric(BC_sub$AC_depth))-sd(as.numeric(BC_sub$AC_depth)),
  mean(as.numeric(AC_sub$BC_depth))-sd(as.numeric(AC_sub$BC_depth)),
  mean(as.numeric(AB_sub$AC_depth))-sd(as.numeric(AB_sub$AC_depth)),
  mean(as.numeric(BC_sub$BC_depth))-sd(as.numeric(BC_sub$BC_depth)),
  mean(as.numeric(AC_sub$AC_depth))-sd(as.numeric(AC_sub$AC_depth)),
  mean(as.numeric(AB_sub$AB_depth))-sd(as.numeric(AB_sub$AB_depth))
)


T1_T2_df<-data.frame(TOPOLOGY, MEAN, YMAX, YMIN)



#Plot means +/= one standard deviation
ggplot(T1_T2_df, aes(x=T1_T2_df$TOPOLOGY, y=MEAN)) + geom_point() + geom_errorbar(width=.1, aes(ymin=YMIN, ymax=YMAX)) + 
  geom_point(shape=21, size=3, fill="white") + 
  ylim(0.2, 1) + 
  xlab("Topology") + 
  ylab("Node depth") +
  geom_text(aes(x=2, y=0.2, label= paste("lambda=", lambda))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#T1 BC vs AC
t.test(as.numeric(BC_sub$AC_depth), as.numeric(AC_sub$BC_depth))

#T1 BC vs AB
t.test(as.numeric(BC_sub$AC_depth), as.numeric(AB_sub$AC_depth))

#T1 AC vs AB
t.test(as.numeric(AC_sub$BC_depth), as.numeric(AB_sub$AC_depth))

#T2 BC vs AC
t.test(as.numeric(BC_sub$BC_depth), as.numeric(AC_sub$AC_depth))

#T2 BC vs AB
t.test(as.numeric(BC_sub$BC_depth), as.numeric(AB_sub$AB_depth))

#T2 AC vs AB
t.test(as.numeric(AC_sub$AC_depth), as.numeric(AB_sub$AB_depth))

#Boxplot with all but no AC comparison
boxplot(
  #T1
    as.numeric(BC_sub$AC_depth), as.numeric(AC_sub$BC_depth), as.numeric(AB_sub$BC_depth), 
    #T2
      as.numeric(BC_sub$BC_depth), as.numeric(AC_sub$AC_depth), as.numeric(AB_sub$AB_depth), outline = FALSE)

#rank sums tests
###T1
#BC vs AC
wilcox.test(as.numeric(BC_sub$AC_depth), as.numeric(AC_sub$BC_depth))

#BC vs AB
wilcox.test(as.numeric(BC_sub$AC_depth), as.numeric(AB_sub$BC_depth))

#AC vs AB
wilcox.test(as.numeric(AC_sub$BC_depth), as.numeric(AB_sub$BC_depth))

###T2
#BC vs AC
wilcox.test(as.numeric(BC_sub$BC_depth), as.numeric(AC_sub$AC_depth))

#BC vs AB
wilcox.test(as.numeric(BC_sub$BC_depth), as.numeric(AB_sub$AB_depth))

#AC vs AB
wilcox.test(as.numeric(AC_sub$AC_depth), as.numeric(AB_sub$AB_depth))
