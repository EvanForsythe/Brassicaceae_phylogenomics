setwd("./R_test_dir/")

#load some packages
install.packages("ape")
library(ape)
install.packages("phytools")
library(phytools)
#load packages for TICR test
install.packages("phylolm")
library(phylolm)

#Read the phlyonet population tree (0 reticulations)
phylonet_tree<-read.tree("phylonet171218.tre")

#unroot tree
phylonet_tree<-unroot(phylonet_tree)

#Plot tree
plot(phylonet_tree, use.edge.length=T)

plot(phylonet_tree, use.edge.length=F)
edgelabels(round(phylonet_tree$edge.length,3), frame="none", adj=c(0.5,0))

#load quartet table
quartet_csv<-read.csv("tableCF_171012.csv")
dim(quartet_csv)

#### TICR ###
#Prelim

prelim = test.tree.preparation(quartet_csv, phylonet_tree)

Ntaxa = length(phylonet_tree$tip.label)

internal.edges = which(phylonet_tree$edge[,2] > Ntaxa)

internal.edges # indices of internal edges: those that need a length

#Run actual test 
res <- test.one.species.tree(quartet_csv,phylonet_tree,prelim,edge.keep=internal.edges)

res[1:6]

#output

#$alpha
#[1] 281.7089
#$minus.pll
#[1] -250.9742
#$X2
#[1] 17.4127
#$chisq.pval
#[1] 0.0005812105
#$chisq.conclusion
#[1] "The chi-square test is significant at level .05,\nbut there is a deficit of outlier quartets (with outlier p-value<=0.01).\nThis pattern does not have a simple evolutionary explanation.\n"

#$outlier.table
#.01 .05  .10 large
#observed 0.00 6.0 0.00  29.0
#expected 0.35 1.4 1.75  31.5

