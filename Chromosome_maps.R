setwd("./R_test_dir/")

#Install bioconductor, which is not available from CRAN
source("https://bioconductor.org/biocLite.R")
#Activate sushi package
biocLite('Sushi')

#read in dataset
full_data<-read.csv("Chrom_map_data170621.csv", head=TRUE)

#Round midpoint up to interger
full_data$Crub_mid<-ceiling(full_data$Crub_mid)
full_data$Atha_mid<-ceiling((full_data$Atha_start+full_data$Atha_stop)/2)

#Find the length of each chrome
max1<-max(subset(full_data, full_data$Crub_chrom==1)$Crub_stop)
max2<-max(subset(full_data, full_data$Crub_chrom==2)$Crub_stop)
max3<-max(subset(full_data, full_data$Crub_chrom==3)$Crub_stop)
max4<-max(subset(full_data, full_data$Crub_chrom==4)$Crub_stop)
max5<-max(subset(full_data, full_data$Crub_chrom==5)$Crub_stop)
max6<-max(subset(full_data, full_data$Crub_chrom==6)$Crub_stop)
max7<-max(subset(full_data, full_data$Crub_chrom==7)$Crub_stop)
max8<-max(subset(full_data, full_data$Crub_chrom==8)$Crub_stop)

#Subset for only BC top
BC_data<-subset(full_data, full_data$Topology_loose=="BC_topology")
#Make data file
BCbed<-data.frame(chrom=BC_data$Crub_chrom, start=BC_data$Crub_start, end=BC_data$Crub_stop, value=BC_data$Bootstrap_Support)

#Subset for only AC top
AC_data<-subset(full_data, full_data$Topology_loose=="AC_topology")
#Make data file
ACbed<-data.frame(chrom=AC_data$Crub_chrom, start=AC_data$Crub_start, end=AC_data$Crub_stop, value=AC_data$Bootstrap_Support)

#Subset for only AB top
AB_data<-subset(full_data, full_data$Topology_loose=="AB_topology")
#Make data file
ABbed<-data.frame(chrom=AB_data$Crub_chrom, start=AB_data$Crub_start, end=AB_data$Crub_stop, value=AB_data$Bootstrap_Support)

#PLOT CHR1
#Start plotting
#Use width to scale by the length of the chromosome
pdf("Chrom1.pdf",width=6,height=1)

#Set graphical parameters for multipanel graph
par(mfrow=c(3,1))
par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
Sushi::plotBedgraph(BCbed,1,1,max1, color = "orange", ymax = 1)
#Add AC
Sushi::plotBedgraph(ACbed,1,1,max1, color = "green", ymax = 1, overlay = FALSE)
#Add AB
Sushi::plotBedgraph(ABbed,1,1,max1, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
dev.off()

#PLOT CHR2
#Start plotting
#Use width to scale by the length of the chromosome
pdf("Chrom2.pdf",width=4.3,height=1)

#Set graphical parameters for multipanel graph
par(mfrow=c(3,1))
par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
Sushi::plotBedgraph(BCbed,2,1,max2, color = "orange", ymax = 1)
#Add AC
Sushi::plotBedgraph(ACbed,2,1,max2, color = "green", ymax = 1, overlay = FALSE)
#Add AB
Sushi::plotBedgraph(ABbed,2,1,max2, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
dev.off()



### NOTE: FOR SOME REASON 3 and 4 are buggy. I provide the work-around below
#PLOT CHR3
#Start plotting
#Use width to scale by the length of the chromosome
#pdf("Chrom3.pdf",width=4.6,height=1)

#Set graphical parameters for multipanel graph
#par(mfrow=c(3,1))
#par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
#Sushi::plotBedgraph(BCbed,3,1,max3, color = "orange", ymax = 1)
#Add AC
#Sushi::plotBedgraph(ACbed,3,1,max3, color = "green", ymax = 1, overlay = FALSE)
#Add AB
#Sushi::plotBedgraph(ABbed,3,1,max3, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
#dev.off()

#PLOT CHR4
#Start plotting
#Use width to scale by the length of the chromosome
#pdf("Chrom4.pdf",width=4.6,height=1)

#Set graphical parameters for multipanel graph
#par(mfrow=c(3,1))
#par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
#Sushi::plotBedgraph(BCbed,4,1,max4, color = "orange", ymax = 1)
#Add AC
#Sushi::plotBedgraph(ACbed,4,1,max4, color = "green", ymax = 1, overlay = FALSE)
#Add AB
#Sushi::plotBedgraph(ABbed,4,1,max4, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
#dev.off()

#PLOT CHR5
#Start plotting
#Use width to scale by the length of the chromosome
#pdf("Chrom5.pdf",width=4.2,height=1)

#Set graphical parameters for multipanel graph
#par(mfrow=c(3,1))
#par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
#Sushi::plotBedgraph(BCbed,5,1,max5, color = "orange", ymax = 1)
#Add AC
#Sushi::plotBedgraph(ACbed,5,1,max5, color = "green", ymax = 1, overlay = FALSE)
#Add AB
#Sushi::plotBedgraph(ABbed,5,1,max5, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
#dev.off()

#PLOT CHR6
#Start plotting
#Use width to scale by the length of the chromosome
pdf("Chrom6.pdf",width=5.1,height=1)

#Set graphical parameters for multipanel graph
par(mfrow=c(3,1))
par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
Sushi::plotBedgraph(BCbed,6,1,max6, color = "orange", ymax = 1)
#Add AC
Sushi::plotBedgraph(ACbed,6,1,max6, color = "green", ymax = 1, overlay = FALSE)
#Add AB
Sushi::plotBedgraph(ABbed,6,1,max6, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
dev.off()

#PLOT CHR7
#Start plotting
#Use width to scale by the length of the chromosome
pdf("Chrom7.pdf",width=5.3,height=1)

#Set graphical parameters for multipanel graph
par(mfrow=c(3,1))
par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
Sushi::plotBedgraph(BCbed,7,1,max7, color = "orange", ymax = 1)
#Add AC
Sushi::plotBedgraph(ACbed,7,1,max7, color = "green", ymax = 1, overlay = FALSE)
#Add AB
Sushi::plotBedgraph(ABbed,7,1,max7, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
dev.off()

#PLOT CHR8
#Start plotting
#Use width to scale by the length of the chromosome
pdf("Chrom8.pdf",width=4.1,height=1)

#Set graphical parameters for multipanel graph
par(mfrow=c(3,1))
par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
Sushi::plotBedgraph(BCbed,8,1,max8, color = "orange", ymax = 1)
#Add AC
Sushi::plotBedgraph(ACbed,8,1,max8, color = "green", ymax = 1, overlay = FALSE)
#Add AB
Sushi::plotBedgraph(ABbed,8,1,max8, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
dev.off()



###### REDO CHROM 3 and 4 #########

only3<-subset(full_data, full_data$Crub_chrom==3)
only4<-subset(full_data, full_data$Crub_chrom==4)
only34<-rbind(only3, only4)

BC_only34<-subset(only34, only34$Topology_loose=="BC_topology")
AC_only34<-subset(only34, only34$Topology_loose=="AC_topology")
AB_only34<-subset(only34, only34$Topology_loose=="AB_topology")

BC_only34$NewStart<-(BC_only34$Crub_mid-500)
BC_only34$NewEnd<-(BC_only34$Crub_mid+500)

AC_only34$NewStart<-(AC_only34$Crub_mid-500)
AC_only34$NewEnd<-(AC_only34$Crub_mid+500)

AB_only34$NewStart<-(AB_only34$Crub_mid-500)
AB_only34$NewEnd<-(AB_only34$Crub_mid+500)

#Make data file
BCbed_34only<-data.frame(chrom=BC_only34$Crub_chrom, start=BC_only34$NewStart, end=BC_only34$NewEnd, value=BC_only34$Bootstrap_Support)
#Make data file
ACbed_34only<-data.frame(chrom=AC_only34$Crub_chrom, start=AC_only34$NewStart, end=AC_only34$NewEnd, value=AC_only34$Bootstrap_Support)
#Make data file
ABbed_34only<-data.frame(chrom=AB_only34$Crub_chrom, start=AB_only34$NewStart, end=AB_only34$NewEnd, value=AB_only34$Bootstrap_Support)

#PLOT CHR3
#Start plotting
#Use width to scale by the length of the chromosome
pdf("Chrom3_new.pdf",width=4.6,height=1)

#Set graphical parameters for multipanel graph
par(mfrow=c(3,1))
par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
Sushi::plotBedgraph(BCbed_34only,3,1,max3, color = "orange", ymax = 1)
#Add AC
Sushi::plotBedgraph(ACbed_34only,3,1,max3, color = "green", ymax = 1, overlay = FALSE)
#Add AB
Sushi::plotBedgraph(ABbed_34only,3,1,max3, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
dev.off()

#PLOT CHR4
#Start plotting
#Use width to scale by the length of the chromosome
pdf("Chrom4_new.pdf",width=4.6,height=1)

#Set graphical parameters for multipanel graph
par(mfrow=c(3,1))
par(mar=c(0,0,0,0), oma= c(0.1,0.1,0.1,0.1))
#plot the chromosome plots
Sushi::plotBedgraph(BCbed_34only,4,1,max4, color = "orange", ymax = 1)
#Add AC
Sushi::plotBedgraph(ACbed_34only,4,1,max4, color = "green", ymax = 1, overlay = FALSE)
#Add AB
Sushi::plotBedgraph(ABbed_34only,4,1,max4, color = "purple", ymax = 1, overlay = FALSE)

#End plotting
dev.off()



##################################
###                           ####
###   PLOTTING ALL MARKERS    ####
###                           ####
##################################


# Make bed file 
#Allbed<-data.frame(chrom=full_data$Crub_chrom, start=full_data$Crub_mid - 500, end=full_data$Crub_mid + 500, value=100)

## TOGGLE THIS ON FOR ATHAL PLOTS

# Make bed file 
Allbed<-data.frame(chrom=full_data$Atha_chrom, start=full_data$Atha_mid - 500, end=full_data$Atha_mid + 500, value=100)

# Make bed file 
Allbed<-data.frame(chrom=full_data$Atha_chrom, start=full_data$Atha_start, end=full_data$Atha_stop, value=1)



#Start plotting
#df("AllMarkers_AllChromes.pdf",width=6,height=8)

#Set graphical parameters for multipanel graph
par(mfrow=c(8,1))
par(mar=c(0.1,0.1,0.1,0.1), oma= c(0.1,0.1,0.1,0.1))

#plot the chromosome plots
Sushi::plotBedgraph(Allbed,1,1,max1, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,2,1,max2, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,3,1,max3, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,4,1,max4, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,5,1,max5, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,6,1,max6, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,7,1,max7, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,8,1,max8, color = "black", ymax = 1)

#End plotting
dev.off()





###CSC###
library('plyr')

#IMPORT SINGLE COPY STATUS STUFF
SCS_data<-read.csv("Single_copy_status.csv", header = TRUE)
joined_data<-join(full_data, SCS_data, type= "inner")

### TOGGLE THIS FOR  Consesus single copy
data_CSC<-subset(joined_data, joined_data$single_copy_status == "ortho_desmet" | joined_data$single_copy_status == "ortho_duarte" | joined_data$single_copy_status == "ret_dup" | joined_data$single_copy_status == "all_three")

#Find the length of each chrome
csc_max1<-max(subset(data_CSC, data_CSC$Crub_chrom==1)$Crub_stop)
csc_max2<-max(subset(data_CSC, data_CSC$Crub_chrom==2)$Crub_stop)
csc_max3<-max(subset(data_CSC, data_CSC$Crub_chrom==3)$Crub_stop)
csc_max4<-max(subset(data_CSC, data_CSC$Crub_chrom==4)$Crub_stop)
csc_max5<-max(subset(data_CSC, data_CSC$Crub_chrom==5)$Crub_stop)
csc_max6<-max(subset(data_CSC, data_CSC$Crub_chrom==6)$Crub_stop)
csc_max7<-max(subset(data_CSC, data_CSC$Crub_chrom==7)$Crub_stop)
csc_max8<-max(subset(data_CSC, data_CSC$Crub_chrom==8)$Crub_stop)

# Make bed file 
CSCbed<-data.frame(chrom=data_CSC$Crub_chrom, start=data_CSC$Crub_mid - 500, end=data_CSC$Crub_mid + 500, value=100)


#Start PDF Consesus single copy
pdf("CSCMarkers_AllChromes.pdf",width=6,height=8)

#Set graphical parameters for multipanel graph
par(mfrow=c(8,1))
par(mar=c(0.1,0.1,0.1,0.1), oma= c(0.1,0.1,0.1,0.1))

#plot the chromosome plots
Sushi::plotBedgraph(CSCbed,1,1,csc_max1, color = "black", y_max = 1)

Sushi::plotBedgraph(CSCbed,2,1,csc_max2, color = "black", y_max = 1)

Sushi::plotBedgraph(CSCbed,3,1,csc_max3, color = "black", y_max = 1)

Sushi::plotBedgraph(CSCbed,4,1,csc_max4, color = "black", y_max = 1)

Sushi::plotBedgraph(CSCbed,5,1,csc_max5, color = "black", y_max = 1)

Sushi::plotBedgraph(CSCbed,6,1,csc_max6, color = "black", y_max = 1)

Sushi::plotBedgraph(CSCbed,7,1,csc_max7, color = "black", y_max = 1)

Sushi::plotBedgraph(CSCbed,8,1,csc_max8, color = "black", y_max = 1)

#End plotting
dev.off()



### ATHAL
data1<-subset(full_data, full_data$Atha_chrom==1)
bed1<-data.frame(chrom=data1$Atha_chrom, start=data1$Atha_start, end=data1$Atha_stop, value=100)

# Make bed file 
Allbed<-data.frame(chrom=full_data$Atha_chrom, start=full_data$Atha_mid - 600, end=full_data$Atha_mid + 500, value=100)

# Make bed file 
#Allbed<-data.frame(chrom=full_data$Atha_chrom, start=full_data$Atha_start, end=full_data$Atha_stop, value=1)


#Start plotting
#Find the length of each chrome
ATmax1<-max(subset(full_data, full_data$Atha_chrom==1)$Atha_stop)
ATmax2<-max(subset(full_data, full_data$Atha_chrom==2)$Atha_stop)
ATmax3<-max(subset(full_data, full_data$Atha_chrom==3)$Atha_stop)
ATmax4<-max(subset(full_data, full_data$Atha_chrom==4)$Atha_stop)
ATmax5<-max(subset(full_data, full_data$Atha_chrom==5)$Atha_stop)


#TRY A WEIRD SUBSET STEP AS A WORK AROUND FOR ANNOYING BUG
bed1<-subset(Allbed, Allbed$chrom==1)


pdf("ATHA_AllMarkers_AllChromes.pdf",width=6,height=8)

par(mfrow=c(5,1))

par(mar=c(0.1,0.1,0.1,0.1), oma= c(0.1,0.1,0.1,0.1))

#plot the chromosome plots
par(mfrow=c(1,1))
Sushi::plotBedgraph(bed1,1,2,ATmax1, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,2,2,ATmax2, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,3,2,ATmax3, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,4,2,ATmax4, color = "black", ymax = 1)

Sushi::plotBedgraph(Allbed,5,1,ATmax5, color = "black", ymax = 1)

#End plotting
dev.off()



