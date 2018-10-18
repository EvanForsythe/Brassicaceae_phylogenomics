setwd("./R_test_dir/")

#Install and load packages
install.packages("boot")
library(boot)
install.packages("resample")
library(resample)
install.packages("ggplot2")
library(ggplot2)

#read in dataframe
data<-read.csv("All_nuctrees170217.csv", header=TRUE)

#filter out "other topologies"
data_noOthers<-subset(data, data$Topology_loose != "Other_topology")

####This is to filter for conservatively single-copy gene families only
#MAKE SURE TO COMMENT THIS OUT FOR FULL ANALYSIS
data_noOthers<-subset(data_noOthers, data_noOthers$single_copy_status == "ortho_desmet" | data_noOthers$single_copy_status == "ortho_duarte" | data_noOthers$single_copy_status == "ret_dup" | data_noOthers$single_copy_status == "all_three")

#Find counts of each topology
table<-table(data_noOthers$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-data_noOthers$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_full<-c("All", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#############################################################
###Now do the Dstat calculation for individual chromosomes###
#############################################################

#CHROM 1#
chrom1_data<-subset(data_noOthers, Crub_chrom == 1)

table<-table(chrom1_data$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom1_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom1<-c("Chr1", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 2#
chrom2_data<-subset(data_noOthers, Crub_chrom == 2)

table<-table(chrom2_data$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom2_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom2<-c("Chr2", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 3#
chrom3_data<-subset(data_noOthers, Crub_chrom == 3)

table<-table(chrom3_data$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom3_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom3<-c("Chr3", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 4#
chrom4_data<-subset(data_noOthers, Crub_chrom == 4)

table<-table(chrom4_data$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom4_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom4<-c("Chr4", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 5#
chrom5_data<-subset(data_noOthers, Crub_chrom == 5)

table<-table(chrom5_data$Topology_loose)
AB_count<-unname(table["AB_topology"])
#AB_count<-unname(AB_count)
AC_count<-unname(table["AC_topology"])
#AC_count<-unname(AC_count)
BC_count<-unname(table["BC_topology"])
#BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom5_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom5<-c("Chr5", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 6#
chrom6_data<-subset(data_noOthers, Crub_chrom == 6)

table<-table(chrom6_data$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom6_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom6<-c("Chr6", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 7#
chrom7_data<-subset(data_noOthers, Crub_chrom == 7)

table<-table(chrom7_data$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom7_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom7<-c("Chr7", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 8#
chrom8_data<-subset(data_noOthers, Crub_chrom == 8)

table<-table(chrom8_data$Topology_loose)
AB_count<-table["AB_topology"]
AB_count<-unname(AB_count)
AC_count<-table["AC_topology"]
AC_count<-unname(AC_count)
BC_count<-table["BC_topology"]
BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-chrom8_data$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom8<-c("Chr8", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#Concat the results from all the chromosomes
cat_out<-c(out_full, out_chrom1, out_chrom2, out_chrom3, out_chrom4, out_chrom5, out_chrom6, out_chrom7, out_chrom8)

#convert into dataframe
cat_out_df<- data.frame(matrix(unlist(cat_out), nrow=9, byrow=TRUE))

#append column names and rownames
names(cat_out_df) <- c("Chromosome", "AC_trees", "AB_trees", "Dstat_obs", "SD", "Z_score", "P_value")
#rownames(cat_out_df) <- c("All", "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8")

#########################################################
##########    Plot the D-stat figure   ##################
#########################################################

ymin<-as.numeric(paste(cat_out_df$Dstat_obs))-as.numeric(paste(cat_out_df$SD))
ymax<-as.numeric(paste(cat_out_df$Dstat_obs))+as.numeric(paste(cat_out_df$SD))
Dstat<-as.numeric(paste(cat_out_df$Dstat_obs))

#command for printing pie to pdf
#pdf(file="Jack_BC_170626.pdf",width=6,height=4)

#command for printing pie to pdf TOGGLE FOR CSC
pdf(file="Jack_BC_170626_CSC.pdf",width=6,height=4)

#Plot with All as dot
ggplot(cat_out_df, aes(x=cat_out_df$Chromosome, y=Dstat)) + geom_point() + geom_errorbar(width=.1, aes(ymin=ymin, ymax=ymax)) + 
  geom_hline(aes(yintercept=0), color="red", linetype="dashed") + 
  geom_point(shape=21, size=c(7, 3, 3, 3, 3, 3, 3, 3, 3), fill=c("dark gray", "dark gray", "dark gray", "dark gray", "white", "dark gray", "dark gray", "dark gray", "dark gray")) + 
  ylim(-0.1, 0.7) + 
  xlab("Chromosome") + 
  ylab("Patterson's D-statistic") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#end the printing to PDF
dev.off()

#################################################################
################      Bootstrapping     #########################
#################################################################

#Find counts of each topology
table<-table(data_noOthers$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)
#This is the numerator for calculating the proportion of gene flow
Num4propGF<-(AC_count - AB_count)
Num4propGF<-unname(Num4propGF)
#Num4prop = 229 (17-01-06)


#Create function for bootstrapping
Boot_Dstat <- function(X=data_noOthers$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_full<-c("All", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#############################################################
###Now do the Dstat calculation for individual chromosomes###
#############################################################

#CHROM 1#
chrom1_data<-subset(data_noOthers, Crub_chrom == 1)

table<-table(chrom1_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom1_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom1<-c("Chr1", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 2#
chrom2_data<-subset(data_noOthers, Crub_chrom == 2)

table<-table(chrom2_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom2_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom2<-c("Chr2", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 3#
chrom3_data<-subset(data_noOthers, Crub_chrom == 3)

table<-table(chrom3_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom3_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom3<-c("Chr3", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 4#
chrom4_data<-subset(data_noOthers, Crub_chrom == 4)

table<-table(chrom4_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom4_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom4<-c("Chr4", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 5#
#CHROM 5#
chrom5_data<-subset(data_noOthers, Crub_chrom == 5)

table<-table(chrom5_data$Topology_loose)
AB_count<-unname(table["AB_topology"])
#AB_count<-unname(AB_count)
AC_count<-unname(table["AC_topology"])
#AC_count<-unname(AC_count)
BC_count<-unname(table["BC_topology"])
#BC_count<-unname(BC_count)

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom5_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom5<-c("Chr5", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 6#
chrom6_data<-subset(data_noOthers, Crub_chrom == 6)

table<-table(chrom6_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom6_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom6<-c("Chr6", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 7#
chrom7_data<-subset(data_noOthers, Crub_chrom == 7)

table<-table(chrom7_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom7_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom7<-c("Chr7", AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 8#
chrom8_data<-subset(data_noOthers, Crub_chrom == 8)

table<-table(chrom8_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom8_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom8<-c("Chr8", AC_count, AB_count, Dstat_observed, SD, z, P_value)


#Concat the results from all the chromosomes
cat_out<-c(out_full, out_chrom1, out_chrom2, out_chrom3, out_chrom4, out_chrom5, out_chrom6, out_chrom7, out_chrom8)

#convert into dataframe
cat_out_df<- data.frame(matrix(unlist(cat_out), nrow=9, byrow=TRUE))

#append column names and rownames
names(cat_out_df) <- c("Chromosome", "AC_trees", "AB_trees", "Dstat_obs", "SD", "Z_score", "P_value")
#rownames(cat_out_df) <- c("All", "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8")

ymin<-as.numeric(paste(cat_out_df$Dstat_obs))-as.numeric(paste(cat_out_df$SD))
ymax<-as.numeric(paste(cat_out_df$Dstat_obs))+as.numeric(paste(cat_out_df$SD))
Dstat<-as.numeric(paste(cat_out_df$Dstat_obs))

#command for printing pie to pdf
#pdf(file="Boot_BC_170626.pdf",width=6,height=4)
#command for printing pie to pdf TOGGLE FOR CSC
pdf(file="Boot_BC_170626_CSC.pdf",width=6,height=4)

ggplot(cat_out_df, aes(x=cat_out_df$Chromosome, y=Dstat)) + geom_point() + geom_errorbar(width=.1, aes(ymin=ymin, ymax=ymax)) + 
  geom_hline(aes(yintercept=0), color="red", linetype="dashed") + 
  geom_point(shape=21, size=c(7, 3, 3, 3, 3, 3, 3, 3, 3), fill=c("dark gray", "dark gray", "dark gray", "dark gray", "white", "dark gray", "dark gray", "dark gray", "dark gray")) + 
  ylim(-0.1, 0.7) + 
  xlab("Chromosome") + 
  ylab("Patterson's D-statistic") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#end the printing to PDF
dev.off()

write.csv(cat_out_df, file = "boot_out170806.csv")
