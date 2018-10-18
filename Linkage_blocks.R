#This is an R script for asking if linear clustering is observed among genes with the same topology
#we use two different meaures of clustering, 
	#(1) genes of like-topology within 10kb of eachother
	#(2) genes of like-topology neighboring eachother
	
#Observed values are compared against distributions calculated from chromosome maps with topologies randomized


#set working directory
setwd("./R_test_dir/")

install.packages("data.table")
library(data.table)
install.packages("dplyr")
library(dplyr)


#Running entire script below here, but use comment to choose correct topology
loci<-read.csv("TopologyResults4linkage170518.csv", header = TRUE)

### CALCULATE GENES WITHIN 10KB OF LIKE-TOPOLOGY ###
n<-nrow(loci)

###IMPORTANT: remember to toggle this to the desired topology
#keeper_top<-"BC_topology"
#keeper_top<-"AC_topology"
keeper_top<-"AB_topology"


#Create function for randomizing 
neighbor_distance <- function(X=loci){
  x.rando<-data.frame(X$Crub_chrom, X$Crub_mid)
  x.rando$Topology_loose<-sample(X$Topology_loose, n, replace = FALSE)
  sub_rando<-subset(x.rando, x.rando$Topology_loose==keeper_top)
  setDT(sub_rando)
  sub_rando[, lag_Crub_mid := c(NA, X.Crub_mid[-.N]), by = X.Crub_chrom]
  sub_rando[, diff := X.Crub_mid - lag_Crub_mid]
  nrow(subset(sub_rando, sub_rando$diff<10001))
}

#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
rando_replicates <- replicate(N, neighbor_distance())

#Calculate the number of "keeper top" genes 10000 nt or closer to another AC gene
sub_loci<-subset(loci, loci$Topology_loose==keeper_top)
setDT(sub_loci)
sub_loci[, lag_Crub_mid := c(NA, Crub_mid[-.N]), by = Crub_chrom]
sub_loci[, diff := Crub_mid - lag_Crub_mid]
obs_close<-nrow(subset(sub_loci, sub_loci$diff<10001))

plot(density(rando_replicates)#, xlim=c(4700,5000)
     )
abline(v=obs_close)

a<-obs_close
s<-sd(rando_replicates)
n<-length(rando_replicates)
xbar<-mean(rando_replicates)
z<-(xbar-a)/s
2*pnorm(-abs(z))


###  CALCULATE NEIGHBORING GENES WITH SAME TOPOLOGY  ###
library(dplyr)
library(magrittr)

loci$backback_top<-loci$Topology_loose
setDT(loci)
equal_prev <- function(x) {
  x %>% equals(lag(x, default = x[1])) %>% not %>% as.numeric
}

top_lag<-loci %>%
  group_by(Crub_chrom) %>%
  mutate_each(funs(equal_prev),starts_with("Topology_loo"))

obs_full<-nrow(top_lag)-sum(top_match$Topology_loose)

sub_keeper<-subset(top_lag, top_lag$backback_top==keeper_top)

obs_keeper<-nrow(sub_keeper)-sum(sub_keeper$Topology_loose)

### RANDOMIZED REPLICATIONS ###
n<-nrow(loci)

#Create function for randomizing 
neighbor_same <- function(X=loci){
  x.rando<-data.frame(Crub_chrom=X$Crub_chrom)
  x.rando$Topology_loose<-sample(X$Topology_loose, n, replace = FALSE)
  setDT(x.rando)

  top_lag_temp<-x.rando %>%
    group_by(Crub_chrom) %>%
    mutate_each(funs(equal_prev),starts_with("Topology_loo"))
  
  nrow(top_lag_temp)-sum(top_lag_temp$Topology_loose)
}

#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
rando_replicates_same_full <- replicate(N, neighbor_same())


### Plot rando distribution
plot(density(rando_replicates_same), xlim=c(7400, 7700))

abline(v=obs_full)

a<-obs_full
s<-sd(rando_replicates_same)
n<-length(rando_replicates_same)
xbar<-mean(rando_replicates_same)
z<-(xbar-a)/s
2*pnorm(-abs(z))

#GENERAL NEXT DOOR NEIGHBORS (FULL)

### Randomized with per-topology scores ###
###########################################

###IMPORTANT: remeber to toggle this to the desired topology
#keeper_top<-"AC_topology"
#keeper_top<-"AB_topology"
#keeper_top<-"BC_topology"

n<-nrow(loci)
#Create function for randomizing 
neighbor_same2 <- function(X=loci){
  x.rando<-data.frame(Crub_chrom=X$Crub_chrom)
  x.rando$Topology_loose<-sample(X$Topology_loose, n, replace = FALSE)
  x.rando$BACKTopology_loose<-x.rando$Topology_loose
  setDT(x.rando)

  top_lag_temp<-x.rando %>%
    group_by(Crub_chrom) %>%
    mutate_each(funs(equal_prev),starts_with("Topology_loo"))
  
  setDF(top_lag_temp)
  sub_temp<-subset(top_lag_temp, top_lag_temp$BACKTopology_loose==keeper_top)
  nrow(sub_temp)-sum(sub_temp$Topology_loose)
}

#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
bytop_reps<- replicate(N, neighbor_same2())

### Plot rando distribution
plot(density(bytop_reps), xlim=c(0, 50)
     )

abline(v=obs_keeper)

a<-obs_keeper
s<-sd(bytop_reps)
n<-length(bytop_reps)
xbar<-mean(bytop_reps)
z<-(xbar-a)/s
2*pnorm(-abs(z))



