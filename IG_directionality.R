setwd("./R_test_dir/")
library(plyr)

#Read in dS results
x<-read.csv("dS_results_worganelle.csv", header = TRUE)

#read in topology results
y<-read.csv("topology_results170509.csv", header = TRUE)

#Join dS values and topologies
z<-join(x, y, type = "inner")

#calculate average dS between clades
z$dS_BC<-((z$CrBs+z$CgBs)/2)
z$dS_AC<-((z$AtCr+z$AtCg+z$AlCr+z$AlCg)/4)
z$dS_AB<-((z$AtBs+z$AlBs)/2)

#calculate average "T1" distance with dS
#This is the distance from any one species to the other two species
z$dS_A_BC<-((z$AtCr+z$AtCg+z$AlCr+z$AlCg+z$AtBs+z$AlBs)/6)
z$dS_B_AC<-((z$CrBs+z$CgBs+z$AtBs+z$AlBs)/4)
z$dS_C_AB<-((z$AtCr+z$AtCg+z$AlCr+z$AlCg+z$CrBs+z$CgBs)/6)


#subset data to exclude rows inwhich average dS is greater than or equal to 0
z_sub<-subset(z, z$dS_BC<=1 & z$dS_AC<=1 & z$dS_AB<=1 & z$dS_BC>=0 & z$dS_AC>=0 & z$dS_AB>=0)


#subset data to exclude rows inwhich BS is greater less than 70
#IMPORTANT: comment this out to analyze all topologies


#Subset by topology
BC_sub<-subset(z_sub, z_sub$Topology_loose=="BC_topology")
AC_sub<-subset(z_sub, z_sub$Topology_loose=="AC_topology")
AB_sub<-subset(z_sub, z_sub$Topology_loose=="AB_topology")

#Subset chloro and mito
CP_sub<-subset(z_sub, z_sub$Genome=="cp")
MT_sub<-subset(z_sub, z_sub$Genome=="mt")


#Plot dS for subsetted data
plot(density(na.omit(BC_sub$dS_BC), adjust = 0.5), col="red")
lines(density(na.omit(AC_sub$dS_AC), adjust = 0.5), col="green")
lines(density(na.omit(AB_sub$dS_AB), adjust = 0.5), col="blue")

#plot same dS comparing topolgies
plot(density(na.omit(BC_sub$dS_AC), adjust = 0.5), col="dark green")
lines(density(na.omit(AC_sub$dS_AC), adjust = 0.5), col="green")

plot(density(na.omit(BC_sub$dS_BC), adjust = 0.5), col="red")
lines(density(na.omit(AC_sub$dS_BC), adjust = 0.5), col="dark red")

plot(density(na.omit(BC_sub$dS_AB), adjust = 0.5), col="dark blue")
lines(density(na.omit(AB_sub$dS_AB), adjust = 0.5), col="blue")

#Plot distance to outgroup by topology
plot(density(na.omit(BC_sub$dS_out), adjust = 0.5), col="red")
lines(density(na.omit(AC_sub$dS_out), adjust = 0.5), col="green")
lines(density(na.omit(AB_sub$dS_out), adjust = 0.5), col="blue")

#Plot dS for subsetted data
plot(density(na.omit(BC_sub$dS_AC), adjust = 0.5), col="red")
lines(density(na.omit(AC_sub$dS_AC), adjust = 0.5), col="green")
lines(density(na.omit(AB_sub$dS_AC), adjust = 0.5), col="blue")
lines(density(na.omit(CP_sub$dS_AC), adjust = 0.5), col="black")

#Plot dS for subsetted data
plot(density(na.omit(BC_sub$dS_AC), adjust = 0.5), col="red")
lines(density(na.omit(AC_sub$dS_AC), adjust = 0.5), col="green")
lines(density(na.omit(AB_sub$dS_AC), adjust = 0.5), col="blue")
lines(density(na.omit(CP_sub$dS_AC), adjust = 0.5), col="black")
lines(density(na.omit(MT_sub$dS_AC), adjust = 0.5), col="gray")


####################################
#### Direction of introgression ####
####################################

#dS
#Plot dS AC for subsetted data
plot(density(na.omit(BC_sub$dS_AC)), col="red")
lines(density(na.omit(AC_sub$dS_AC)), col="dark green")
abline(v=median(na.omit(BC_sub$dS_AC)),col="red")
abline(v=median(na.omit(AC_sub$dS_AC)),col="dark green")

#lines(density(na.omit(AB_sub$dS_AC), adjust = 0.5), col="blue")
#lines(density(na.omit(CP_sub$dS_AC), adjust = 0.5), col="black")
#lines(density(na.omit(MT_sub$dS_AC), adjust = 0.5), col="gray")

#T test
t.test(na.omit(BC_sub$dS_AC), na.omit(AC_sub$dS_AC))
#wilcox
wilcox.test(na.omit(BC_sub$dS_AC), na.omit(AC_sub$dS_AC))

#Plot dS BC for subsetted data
plot(density(na.omit(BC_sub$dS_BC)), col="red")
lines(density(na.omit(AC_sub$dS_BC)), col="dark green")
abline(v=median(na.omit(BC_sub$dS_BC)),col="red")
abline(v=median(na.omit(AC_sub$dS_BC)),col="dark green")


#lines(density(na.omit(AB_sub$dS_BC), adjust = 0.5), col="blue")
#lines(density(na.omit(CP_sub$dS_BC), adjust = 0.5), col="black")
#lines(density(na.omit(MT_sub$dS_BC), adjust = 0.5), col="gray")

#T test
t.test(na.omit(BC_sub$dS_BC), na.omit(AC_sub$dS_BC))
#wilcox
wilcox.test(na.omit(BC_sub$dS_BC), na.omit(AC_sub$dS_BC))

#Plot dS AB for subsetted data
plot(density(na.omit(BC_sub$dS_AB)), col="red")
lines(density(na.omit(AC_sub$dS_AB)), col="dark green")
abline(v=median(na.omit(BC_sub$dS_AB)),col="red")
abline(v=median(na.omit(AC_sub$dS_AB)),col="dark green")

#lines(density(na.omit(AB_sub$dS_AB), adjust = 0.5), col="blue")
#lines(density(na.omit(CP_sub$dS_AB), adjust = 0.5), col="black")
#lines(density(na.omit(MT_sub$dS_AB), adjust = 0.5), col="gray")

#T test
t.test(na.omit(BC_sub$dS_AB), na.omit(AC_sub$dS_AB))
#wilcox
wilcox.test(na.omit(BC_sub$dS_AB), na.omit(AC_sub$dS_AB))


#dS boxplots
boxplot(na.omit(BC_sub$dS_AC), na.omit(AC_sub$dS_AC), outline = FALSE,
        main="dS(A-C)", 
        xlab="Topology", 
        ylab="dS distance A-C",
        ylim=c(0,0.5),
        names=c("BC", "AC"))

boxplot(na.omit(BC_sub$dS_BC), na.omit(AC_sub$dS_BC), outline = FALSE,
        main="dS(B-C)", 
        xlab="Topology", 
        ylab="dS distance B-C",
        ylim=c(0,0.5),
        names=c("BC", "AC"))

boxplot(na.omit(BC_sub$dS_AB), na.omit(AC_sub$dS_AB), outline = FALSE,
        main="dS(A-B)", 
        xlab="Topology", 
        ylab="dS distance A-B",
        ylim=c(0,0.5),
        names=c("BC", "AC"))

#Boxplot with all
boxplot(na.omit(BC_sub$dS_AC), na.omit(AC_sub$dS_AC),
        na.omit(BC_sub$dS_BC), na.omit(AC_sub$dS_BC),
        na.omit(BC_sub$dS_AB), na.omit(AC_sub$dS_AB), outline = FALSE)
        outline = FALSE,
        xlab="Topology", 
        ylab="dS distance A-C",
        ylim=c(0,0.5),
        col=c()
        
#### THIS IS THE PLOT USED FOR FIGURE #### 170914      
        #Boxplot with all
        #S1
        boxplot(na.omit(AC_sub$dS_BC), na.omit(BC_sub$dS_BC),
                #S2
                na.omit(AC_sub$dS_AC), na.omit(BC_sub$dS_AC),
                #S3
                na.omit(AC_sub$dS_AB), na.omit(BC_sub$dS_AB), outline = FALSE)
        outline = FALSE,
        xlab="Topology", 
        ylab="dS distance",
        ylim=c(0,0.5),
        col=c()
        

#Wilcox test for signifcance. 17-07-31

#S1-sp vs S1-ig
        wilcox.test(na.omit(AC_sub$dS_BC), na.omit(BC_sub$dS_BC))
        #p-value < 2.2e-16
       
        #S2-sp vs S2-ig
         wilcox.test(na.omit(AC_sub$dS_AC), na.omit(BC_sub$dS_AC))
        #p-value = 2.365e-12
        
        #S3-sp vs S3-ig 
        wilcox.test(na.omit(AC_sub$dS_AB), na.omit(BC_sub$dS_AB))
        #p-value = 0.1056
        
#Density plots for supplemental figure
        #S1
        plot(density(na.omit(AC_sub$dS_BC)), col="dark green", ylim=c(0,10))
        lines(density(na.omit(BC_sub$dS_BC)), col="orange")
        #S2
        plot(density(na.omit(AC_sub$dS_AC)), col="dark green", ylim=c(0,10))
        lines(density(na.omit(BC_sub$dS_AC)), col="orange")
        #S3
        plot(density(na.omit(AC_sub$dS_AB)), col="dark green", ylim=c(0,10))
        lines(density(na.omit(BC_sub$dS_AB)), col="orange")


library(ggplot2)
qplot(dS_BC, dS_AC, data = z_sub, colour = Topology_loose)


plot(density(na.omit(BC_sub$RND_B_C), adjust = 0.5), col="red")
lines(density(na.omit(AC_sub$RND_A_C), adjust = 0.5), col="green")
lines(density(na.omit(CP_sub$RND_A_C), adjust = 0.5), col="black")

plot(density(na.omit(BC_sub$dS_BC), adjust = 0.5), col="red")
lines(density(na.omit(AC_sub$dS_AC), adjust = 0.5), col="green")
lines(density(na.omit(CP_sub$dS_AC), adjust = 0.5), col="black")


#### T1 and T2 ####

#dS
boxplot(
#T1  
  na.omit(BC_sub$dS_A_BC), na.omit(AC_sub$dS_B_AC), na.omit(AB_sub$dS_C_AB),
#T2        
        na.omit(BC_sub$dS_BC), na.omit(AC_sub$dS_AC), na.omit(AB_sub$dS_AB),
        outline = FALSE)





