setwd("./R_test_dir/")
library(plyr)
install.packages("dplyr")
library(dplyr)

#Import topology results
my_data<-read.csv("mydata.csv", header = TRUE)

#Import TAIR data
tair_GO<-read.csv("TAIR_GO_data.csv", header = TRUE)

#Localization aspect only (Cellular component = C)
LocGO<-subset(tair_GO, tair_GO$Aspect=="C")

#Subset for chloroplast localized
LocCP<-subset(LocGO, LocGO$GO_slim=="chloroplast")

#remove duplicated gene names
LocCP_u<-LocCP[!duplicated(LocCP$Name), ]

#Append a "Localization" column onto this file.
LocCP_u$Localization<-"CP_loc"

#Join my data with CP (GO) localization data
My_CP<-join(my_data, LocCP_u, type ="left")

#Count AC top genes
AC_CP<-nrow(subset(My_CP, My_CP$Topology_loose == "AC_topology" & My_CP$Localization == "CP_loc"))
AC_noCP<-(nrow(subset(My_CP, My_CP$Topology_loose == "AC_topology")))-(nrow(subset(My_CP, My_CP$Topology_loose == "AC_topology" & My_CP$Localization == "CP_loc")))

#Count BC top genes
BC_CP<-nrow(subset(My_CP, My_CP$Topology_loose == "BC_topology" & My_CP$Localization == "CP_loc"))
BC_noCP<-(nrow(subset(My_CP, My_CP$Topology_loose == "BC_topology")))-(nrow(subset(My_CP, My_CP$Topology_loose == "BC_topology" & My_CP$Localization == "CP_loc")))

#create a twoway table 
twoway_matrix_fullCP<-matrix(c(AC_CP, AC_noCP, BC_CP, BC_noCP), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_fullCP) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_fullCP) <- c("Chloro_localized", "Non_Chloro_localized")
twoway_table_fullCP<-as.table(twoway_matrix_fullCP)

#Percent AC
propAC_CP<-AC_CP/(AC_CP+AC_noCP)

#Percent BC
propBC_CP<-BC_CP/(BC_CP+BC_noCP)

#Enrichment score
enrich_CP<-(propAC_CP-propBC_CP)/(propAC_CP+propBC_CP)

#Significance
fisher.test(twoway_matrix_fullCP, alternative = "two.sided")

fisher.test(twoway_matrix_fullCP, alternative = "greater")

### CP PATHWAYS ###
#list for ensuring that all of these are CP localized
LocCP_NAME<-data.frame(Name=LocCP_u$Name)

###PHOTOSYNTHESIS
light1<-tair_GO[grep("photosystem", tair_GO$String), ]
light2<-tair_GO[grep("chloroplast ATP synthase", tair_GO$String), ]
light3<-tair_GO[grep("photosynthesis  light", tair_GO$String), ]
All_light<-rbind(light1, light2, light3)
dark1<-tair_GO[grep("reductive pentose-phosphate cycle", tair_GO$String), ]
dark2<-tair_GO[grep("chloroplast ribulose", tair_GO$String), ]
dark3<-tair_GO[grep("carbon fixation", tair_GO$String), ]
dark4<-tair_GO[grep("photosynthesis  dark", tair_GO$String), ]
All_dark<-rbind(dark1, dark2, dark3, dark4)
phot_resp1<-tair_GO[grep("photoresp", tair_GO$String), ]
phot_resp2<-tair_GO[grep("oxidative photosynthetic", tair_GO$String), ]
All_photoresp<-rbind(phot_resp1, phot_resp2)
ALLphotosynth<-tair_GO[grep("photosynthesis", tair_GO$String), ]
photosynth<-rbind(All_light, All_dark, All_photoresp, ALLphotosynth)

#remove duplicated gene names
photosynth_u<-photosynth[!duplicated(photosynth$Name), ]
#Make sure genes are CP localized
photosynth_u<-join(photosynth_u, LocCP_NAME, type="inner")
#Append a "Localization" column onto this file.
photosynth_u$Localization<-"photosynth"
#Join my data with CP (GO) localization data
My_photosynth<-join(my_data, photosynth_u, type ="left")

#Count AC top genes
AC_photosynth<-nrow(subset(My_photosynth, My_photosynth$Topology_loose == "AC_topology" & My_photosynth$Localization == "photosynth"))
AC_nophotosynth<-(nrow(subset(My_photosynth, My_photosynth$Topology_loose == "AC_topology")))-(nrow(subset(My_photosynth, My_photosynth$Topology_loose == "AC_topology" & My_photosynth$Localization == "photosynth")))
#Count BC top genes
BC_photosynth<-nrow(subset(My_photosynth, My_photosynth$Topology_loose == "BC_topology" & My_photosynth$Localization == "photosynth"))
BC_nophotosynth<-(nrow(subset(My_photosynth, My_photosynth$Topology_loose == "BC_topology")))-(nrow(subset(My_photosynth, My_photosynth$Topology_loose == "BC_topology" & My_photosynth$Localization == "photosynth")))

#Table
twoway_matrix_photosynth<-matrix(c(AC_photosynth, AC_nophotosynth, BC_photosynth, BC_nophotosynth), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_photosynth) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_photosynth) <- c("photosynthesis", "Non_photosynthesis")
twoway_table_photosynth<-as.table(twoway_matrix_photosynth)

#Percent AC
propAC_photosynth<-AC_photosynth/(AC_photosynth+AC_nophotosynth)
#Percent BC
propBC_photosynth<-BC_photosynth/(BC_photosynth+BC_nophotosynth)

#Enrichment score
enrich_photosynth<-(propAC_photosynth-propBC_photosynth)/(propAC_photosynth+propBC_photosynth)

#Significance
fisher.test(twoway_matrix_photosynth, alternative = "two.sided")
fisher.test(twoway_matrix_photosynth, alternative = "greater")

###LIGHT REACTIONS
light1<-tair_GO[grep("photosystem", tair_GO$String), ]
light2<-tair_GO[grep("chloroplast ATP synthase", tair_GO$String), ]
light3<-tair_GO[grep("photosynthesis  light", tair_GO$String), ]
All_light<-rbind(light1, light2, light3)
#remove duplicated gene names
All_light_u<-All_light[!duplicated(All_light$Name), ]
#Make sure genes are CP localized
All_light_u<-join(All_light_u, LocCP_NAME, type="inner")
#Append a "Localization" column onto this file.
All_light_u$Localization<-"light"
#Join my data with light reaction localization data
My_light<-join(my_data, All_light_u, type ="left")
#Count AC top genes
AC_light<-nrow(subset(My_light, My_light$Topology_loose == "AC_topology" & My_light$Localization == "light"))
AC_nolight<-(nrow(subset(My_light, My_light$Topology_loose == "AC_topology")))-(nrow(subset(My_light, My_light$Topology_loose == "AC_topology" & My_light$Localization == "light")))
#Count BC top genes
BC_light<-nrow(subset(My_light, My_light$Topology_loose == "BC_topology" & My_light$Localization == "light"))
BC_nolight<-(nrow(subset(My_light, My_light$Topology_loose == "BC_topology")))-(nrow(subset(My_light, My_light$Topology_loose == "BC_topology" & My_light$Localization == "light")))
#Table
twoway_matrix_light<-matrix(c(AC_light, AC_nolight, BC_light, BC_nolight), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_light) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_light) <- c("light", "Non_light")
twoway_table_light<-as.table(twoway_matrix_light)

#Percent AC
propAC_light<-AC_light/(AC_light+AC_nolight)

#Percent BC
propBC_light<-BC_light/(BC_light+BC_nolight)

#Enrichment score
enrich_light<-(propAC_light-propBC_light)/(propAC_light+propBC_light)

#Significance
fisher.test(twoway_matrix_light, alternative = "two.sided")
fisher.test(twoway_matrix_light, alternative = "greater")


#DARK REACTIONS
dark1<-tair_GO[grep("reductive pentose-phosphate cycle", tair_GO$String), ]
dark2<-tair_GO[grep("chloroplast ribulose", tair_GO$String), ]
dark3<-tair_GO[grep("carbon fixation", tair_GO$String), ]
dark4<-tair_GO[grep("photosynthesis  dark", tair_GO$String), ]
All_dark<-rbind(dark1, dark2, dark3, dark4)
#remove duplicated gene names
All_dark_u<-All_dark[!duplicated(All_dark$Name), ]
#Make sure genes are CP localized
All_dark_u<-join(All_dark_u, LocCP_NAME, type="inner")
#Append a "Localization" column onto this file.
All_dark_u$Localization<-"dark"
#Join my data with light reaction localization data
My_dark<-join(my_data, All_dark_u, type ="left")

#Count AC top genes
AC_dark<-nrow(subset(My_dark, My_dark$Topology_loose == "AC_topology" & My_dark$Localization == "dark"))
AC_nodark<-(nrow(subset(My_dark, My_dark$Topology_loose == "AC_topology")))-(nrow(subset(My_dark, My_dark$Topology_loose == "AC_topology" & My_dark$Localization == "dark")))
#Count BC top genes
BC_dark<-nrow(subset(My_dark, My_dark$Topology_loose == "BC_topology" & My_dark$Localization == "dark"))
BC_nodark<-(nrow(subset(My_dark, My_dark$Topology_loose == "BC_topology")))-(nrow(subset(My_dark, My_dark$Topology_loose == "BC_topology" & My_dark$Localization == "dark")))
#Table
twoway_matrix_dark<-matrix(c(AC_dark, AC_nodark, BC_dark, BC_nodark), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_dark) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_dark) <- c("dark", "Non_dark")
twoway_table_dark<-as.table(twoway_matrix_dark)

#Percent AC
propAC_dark<-AC_dark/(AC_dark+AC_nodark)
#Percent BC
propBC_dark<-BC_dark/(BC_dark+BC_nodark)

#Enrichment score
enrich_dark<-(propAC_dark-propBC_dark)/(propAC_dark+propBC_dark)

#Significance
fisher.test(twoway_matrix_dark, alternative = "two.sided")
fisher.test(twoway_matrix_dark, alternative = "greater")

#PHOTORESPIRATION
phot_resp1<-tair_GO[grep("photoresp", tair_GO$String), ]
phot_resp2<-tair_GO[grep("oxidative photosynthetic", tair_GO$String), ]
All_photoresp<-rbind(phot_resp1, phot_resp2)
#remove duplicated gene names
All_photoresp_u<-All_photoresp[!duplicated(All_photoresp$Name), ]
#Make sure genes are CP localized
All_photoresp_u<-join(All_photoresp_u, LocCP_NAME, type="inner")
#Append a "Localization" column onto this file.
All_photoresp_u$Localization<-"photoresp"
#Join my data with light reaction localization data
My_photoresp<-join(my_data, All_photoresp_u, type ="left")

#Count AC top genes
AC_photoresp<-nrow(subset(My_photoresp, My_photoresp$Topology_loose == "AC_topology" & My_photoresp$Localization == "photoresp"))
AC_nophotoresp<-(nrow(subset(My_photoresp, My_photoresp$Topology_loose == "AC_topology")))-(nrow(subset(My_photoresp, My_photoresp$Topology_loose == "AC_topology" & My_photoresp$Localization == "photoresp")))
#Count BC top genes
BC_photoresp<-nrow(subset(My_photoresp, My_photoresp$Topology_loose == "BC_topology" & My_photoresp$Localization == "photoresp"))
BC_nophotoresp<-(nrow(subset(My_photoresp, My_photoresp$Topology_loose == "BC_topology")))-(nrow(subset(My_photoresp, My_photoresp$Topology_loose == "BC_topology" & My_photoresp$Localization == "photoresp")))
#Table
twoway_matrix_photoresp<-matrix(c(AC_photoresp, AC_nophotoresp, BC_photoresp, BC_nophotoresp), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_photoresp) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_photoresp) <- c("photoresp", "Non_photoresp")
twoway_table_photoresp<-as.table(twoway_matrix_photoresp)

#Percent AC
propAC_photoresp<-AC_photoresp/(AC_photoresp+AC_nophotoresp)
#Percent BC
propBC_photoresp<-BC_photoresp/(BC_photoresp+BC_nophotoresp)

#Enrichment score
enrich_photoresp<-(propAC_photoresp-propBC_photoresp)/(propAC_photoresp+propBC_photoresp)

#Significance
fisher.test(twoway_matrix_photoresp, alternative = "two.sided")
fisher.test(twoway_matrix_photoresp, alternative = "greater")


################################################
############   MITOCHONDRIA     ################
################################################

LocMT<-subset(LocGO, LocGO$GO_slim=="mitochondria")
#remove duplicated gene names
LocMT_u<-LocMT[!duplicated(LocMT$Name), ]
#Append a "Localization" column onto this file.
LocMT_u$Localization<-"MT_loc"
#Join my data with MT (GO) localization data
My_MT<-join(my_data, LocMT_u, type ="left")

#Count AC top genes
AC_MT<-nrow(subset(My_MT, My_MT$Topology_loose == "AC_topology" & My_MT$Localization == "MT_loc"))
AC_noMT<-(nrow(subset(My_MT, My_MT$Topology_loose == "AC_topology")))-(nrow(subset(My_MT, My_MT$Topology_loose == "AC_topology" & My_MT$Localization == "MT_loc")))

#Count BC top genes
BC_MT<-nrow(subset(My_MT, My_MT$Topology_loose == "BC_topology" & My_MT$Localization == "MT_loc"))
BC_noMT<-(nrow(subset(My_MT, My_MT$Topology_loose == "BC_topology")))-(nrow(subset(My_MT, My_MT$Topology_loose == "BC_topology" & My_MT$Localization == "MT_loc")))

twoway_matrix_fullMT<-matrix(c(AC_MT, AC_noMT, BC_MT, BC_noMT), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_fullMT) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_fullMT) <- c("Mito_localized", "Non_Mito_localized")
twoway_table_fullMT<-as.table(twoway_matrix_full)

#Percent AC
propAC_MT<-AC_MT/(AC_MT+AC_noMT)
#Percent BC
propBC_MT<-BC_MT/(BC_MT+BC_noMT)

#Enrichment score
enrich_MT<-(propAC_MT-propBC_MT)/(propAC_MT+propBC_MT)

#Significance
fisher.test(twoway_matrix_fullMT, alternative = "two.sided")
fisher.test(twoway_matrix_fullMT, alternative = "greater")

## MT FUNCTIONAL CATEGORIES
#list for ensuring that all of these are MT localized
LocMT_NAME<-data.frame(Name=LocMT_u$Name)

resp1<-tair_GO[grep("cellular respiration", tair_GO$String), ]
cyt1<-tair_GO[grep("cytochrome c", tair_GO$String), ]
glyc1<-tair_GO[grep("glycolytic process", tair_GO$String), ]
Etrans1<-tair_GO[grep("respiratory electron transport", tair_GO$String), ]
Etrans2<-tair_GO[grep("aerobic electron transport", tair_GO$String), ]
Etrans3<-tair_GO[grep("mitochondrial electron transport", tair_GO$String), ]

#Combine all
All_resp<-rbind(resp1, cyt1, glyc1, Etrans1, Etrans2, Etrans3)

#remove duplicated gene names
All_resp_u<-All_resp[!duplicated(All_resp$Name), ]
#Make sure genes are CP localized
All_resp_u<-join(All_resp_u, LocMT_NAME, type="inner")
#Append a "Localization" column onto this file.
All_resp_u$Localization<-"resp"
#Join my data with MT (GO) localization data
My_resp<-join(my_data, All_resp_u, type ="left")

#Count AC top genes
AC_resp<-nrow(subset(My_resp, My_resp$Topology_loose == "AC_topology" & My_resp$Localization == "resp"))
AC_noresp<-(nrow(subset(My_resp, My_resp$Topology_loose == "AC_topology")))-(nrow(subset(My_resp, My_resp$Topology_loose == "AC_topology" & My_resp$Localization == "resp")))

#Count BC top genes
BC_resp<-nrow(subset(My_resp, My_resp$Topology_loose == "BC_topology" & My_resp$Localization == "resp"))
BC_noresp<-(nrow(subset(My_resp, My_resp$Topology_loose == "BC_topology")))-(nrow(subset(My_resp, My_resp$Topology_loose == "BC_topology" & My_resp$Localization == "resp")))

twoway_matrix_resp<-matrix(c(AC_resp, AC_noresp, BC_resp, BC_noresp), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_resp) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_resp) <- c("resp", "Non_resp")
twoway_table_resp<-as.table(twoway_matrix_resp)

#Percent AC
propAC_resp<-AC_resp/(AC_resp+AC_noresp)
#Percent BC
propBC_resp<-BC_resp/(BC_resp+BC_noresp)

#Enrichment score
enrich_resp<-(propAC_resp-propBC_resp)/(propAC_resp+propBC_resp)

#Significance
fisher.test(twoway_matrix_resp, alternative = "two.sided")
fisher.test(twoway_matrix_resp, alternative = "greater")


###  CYTOCHROME C  ###
#remove duplicated gene names
cyt1_u<-cyt1[!duplicated(cyt1$Name), ]
#Make sure genes are CP localized
cyt1_u<-join(cyt1_u, LocMT_NAME, type="inner")
#Append a "Localization" column onto this file.
cyt1_u$Localization<-"cyt"
#Join my data with MT (GO) localization data
My_cyt<-join(my_data, cyt1_u, type ="left")

#Count AC top genes
AC_cyt<-nrow(subset(My_cyt, My_cyt$Topology_loose == "AC_topology" & My_cyt$Localization == "cyt"))
AC_nocyt<-(nrow(subset(My_cyt, My_cyt$Topology_loose == "AC_topology")))-(nrow(subset(My_cyt, My_cyt$Topology_loose == "AC_topology" & My_cyt$Localization == "cyt")))

#Count BC top genes
BC_cyt<-nrow(subset(My_cyt, My_cyt$Topology_loose == "BC_topology" & My_cyt$Localization == "cyt"))
BC_nocyt<-(nrow(subset(My_cyt, My_cyt$Topology_loose == "BC_topology")))-(nrow(subset(My_cyt, My_cyt$Topology_loose == "BC_topology" & My_cyt$Localization == "cyt")))

twoway_matrix_cyt<-matrix(c(AC_cyt, AC_nocyt, BC_cyt, BC_nocyt), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_cyt) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_cyt) <- c("cyt", "Non_cyt")
twoway_table_cyt<-as.table(twoway_matrix_cyt)

#Percent AC
propAC_cyt<-AC_cyt/(AC_cyt+AC_nocyt)
#Percent BC
propBC_cyt<-BC_cyt/(BC_cyt+BC_nocyt)

#Enrichment score
enrich_cyt<-(propAC_cyt-propBC_cyt)/(propAC_cyt+propBC_cyt)

#Significance
fisher.test(twoway_matrix_cyt, alternative = "two.sided")
fisher.test(twoway_matrix_cyt, alternative = "less")

###  Glycolitic process  ###
#remove duplicated gene names
glyc1_u<-glyc1[!duplicated(glyc1$Name), ]
#Make sure genes are CP localized
glyc1_u<-join(glyc1_u, LocMT_NAME, type="inner")
#Append a "Localization" column onto this file.
glyc1_u$Localization<-"glyc"
#Join my data with MT (GO) localization data
My_glyc<-join(my_data, glyc1_u, type ="left")

#Count AC top genes
AC_glyc<-nrow(subset(My_glyc, My_glyc$Topology_loose == "AC_topology" & My_glyc$Localization == "glyc"))
AC_noglyc<-(nrow(subset(My_glyc, My_glyc$Topology_loose == "AC_topology")))-(nrow(subset(My_glyc, My_glyc$Topology_loose == "AC_topology" & My_glyc$Localization == "glyc")))

#Count BC top genes
BC_glyc<-nrow(subset(My_glyc, My_glyc$Topology_loose == "BC_topology" & My_glyc$Localization == "glyc"))
BC_noglyc<-(nrow(subset(My_glyc, My_glyc$Topology_loose == "BC_topology")))-(nrow(subset(My_glyc, My_glyc$Topology_loose == "BC_topology" & My_glyc$Localization == "glyc")))

twoway_matrix_glyc<-matrix(c(AC_glyc, AC_noglyc, BC_glyc, BC_noglyc), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_glyc) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_glyc) <- c("glyc", "Non_glyc")
twoway_table_glyc<-as.table(twoway_matrix_glyc)

#Percent AC
propAC_glyc<-AC_glyc/(AC_glyc+AC_noglyc)
#Percent BC
propBC_glyc<-BC_glyc/(BC_glyc+BC_noglyc)

#Enrichment score
enrich_glyc<-(propAC_glyc-propBC_glyc)/(propAC_glyc+propBC_glyc)

#Significance
fisher.test(twoway_matrix_glyc, alternative = "two.sided")
fisher.test(twoway_matrix_glyc, alternative = "greater")

### Electron transport ###
All_Etrans<-rbind(Etrans1, Etrans2, Etrans3)
#remove duplicated gene names
All_Etrans_u<-All_Etrans[!duplicated(All_Etrans$Name), ]
#Make sure genes are CP localized
All_Etrans_u<-join(All_Etrans_u, LocMT_NAME, type="inner")
#Append a "Localization" column onto this file.
All_Etrans_u$Localization<-"etrans"
#Join my data with MT (GO) localization data
My_etrans<-join(my_data, All_Etrans_u, type ="left")

#Count AC top genes
AC_etrans<-nrow(subset(My_etrans, My_etrans$Topology_loose == "AC_topology" & My_etrans$Localization == "etrans"))
AC_noetrans<-(nrow(subset(My_etrans, My_etrans$Topology_loose == "AC_topology")))-(nrow(subset(My_etrans, My_etrans$Topology_loose == "AC_topology" & My_etrans$Localization == "etrans")))

#Count BC top genes
BC_etrans<-nrow(subset(My_etrans, My_etrans$Topology_loose == "BC_topology" & My_etrans$Localization == "etrans"))
BC_noetrans<-(nrow(subset(My_etrans, My_etrans$Topology_loose == "BC_topology")))-(nrow(subset(My_etrans, My_etrans$Topology_loose == "BC_topology" & My_etrans$Localization == "etrans")))

#create twoway table
twoway_matrix_etrans<-matrix(c(AC_etrans, AC_noetrans, BC_etrans, BC_noetrans), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_etrans) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_etrans) <- c("etrans", "Non_etrans")
twoway_table_etrans<-as.table(twoway_matrix_etrans)

#Percent AC
propAC_etrans<-AC_etrans/(AC_etrans+AC_noetrans)
#Percent BC
propBC_etrans<-BC_etrans/(BC_etrans+BC_noetrans)

#Enrichment score
enrich_etrans<-(propAC_etrans-propBC_etrans)/(propAC_etrans+propBC_etrans)

#Significance
fisher.test(twoway_matrix_etrans, alternative = "two.sided")
fisher.test(twoway_matrix_etrans, alternative = "greater")


##############################################
#####      PPR genes          ################
##############################################
#List of PPR genes was obtained from:
  #(Lurin et al., 2004) Genome-wide analysis of Arabidopsis pentatricopeptide repeat proteins reveals their essential role in organelle biogenesis

#read in PPR gene list
PPR_genes<-read.csv("PPR_genes", header = TRUE)

PPR_genes$Localization<-"PPR"

#remove duplicates
PPR_genes_u<-PPR_genes[!duplicated(PPR_genes$Name), ]

My_PPR<-join(my_data, PPR_genes_u, type="left")

#Count AC top genes
AC_PPR<-nrow(subset(My_PPR, My_PPR$Topology_loose == "AC_topology" & My_PPR$Localization == "PPR"))
AC_noPPR<-(nrow(subset(My_PPR, My_PPR$Topology_loose == "AC_topology")))-(nrow(subset(My_PPR, My_PPR$Topology_loose == "AC_topology" & My_PPR$Localization == "PPR")))

#Count BC top genes
BC_PPR<-nrow(subset(My_PPR, My_PPR$Topology_loose == "BC_topology" & My_PPR$Localization == "PPR"))
BC_noPPR<-(nrow(subset(My_PPR, My_PPR$Topology_loose == "BC_topology")))-(nrow(subset(My_PPR, My_PPR$Topology_loose == "BC_topology" & My_PPR$Localization == "PPR")))

#create twoway table
twoway_matrix_PPR<-matrix(c(AC_PPR, AC_noPPR, BC_PPR, BC_noPPR), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_PPR) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_PPR) <- c("PPR", "Non_PPR")
twoway_table_PPR<-as.table(twoway_matrix_PPR)

#Percent AC
propAC_PPR<-AC_PPR/(AC_PPR+AC_noPPR)
#Percent BC
propBC_PPR<-BC_PPR/(BC_PPR+BC_noPPR)

#Enrichment score
enrich_PPR<-(propAC_PPR-propBC_PPR)/(propAC_PPR+propBC_PPR)

#Significance
fisher.test(twoway_matrix_PPR, alternative = "two.sided")
fisher.test(twoway_matrix_PPR, alternative = "greater")


##############   NON CYTOPLASMIC GENE   ###################
#Localization aspect only (Cellular component = C)

LocNUC<-subset(LocGO, LocGO$GO_slim=="nucleus")

#remove duplicated gene names
LocNUC_u<-LocNUC[!duplicated(LocNUC$Name), ]
#Append a "Localization" column onto this file.
LocNUC_u$Localization<-"NUC_loc"
#Join my data with NUC (GO) localization data
My_NUC<-join(my_data, LocNUC_u, type ="left")

#Count AC top genes
AC_NUC<-nrow(subset(My_NUC, My_NUC$Topology_loose == "AC_topology" & My_NUC$Localization == "NUC_loc"))
AC_noNUC<-(nrow(subset(My_NUC, My_NUC$Topology_loose == "AC_topology")))-(nrow(subset(My_NUC, My_NUC$Topology_loose == "AC_topology" & My_NUC$Localization == "NUC_loc")))
#Count BC top genes
BC_NUC<-nrow(subset(My_NUC, My_NUC$Topology_loose == "BC_topology" & My_NUC$Localization == "NUC_loc"))
BC_noNUC<-(nrow(subset(My_NUC, My_NUC$Topology_loose == "BC_topology")))-(nrow(subset(My_NUC, My_NUC$Topology_loose == "BC_topology" & My_NUC$Localization == "NUC_loc")))

twoway_matrix_fullNUC<-matrix(c(AC_NUC, AC_noNUC, BC_NUC, BC_noNUC), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_fullNUC) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_fullNUC) <- c("Nuc_localized", "Non_Nuc_localized")
twoway_table_fullNUC<-as.table(twoway_matrix_fullNUC)

#Percent AC
propAC_NUC<-AC_NUC/(AC_NUC+AC_noNUC)
#Percent BC
propBC_NUC<-BC_NUC/(BC_NUC+BC_noNUC)

#Enrichment score
enrich_NUC<-(propAC_NUC-propBC_NUC)/(propAC_NUC+propBC_NUC)

#Significance
fisher.test(twoway_matrix_fullNUC, alternative = "two.sided")
fisher.test(twoway_matrix_fullNUC, alternative = "less")

#NUCLEAR LOCALIZED COMPLEXES
#list for ensuring that all of these are MT localized
LocNUC_NAME<-data.frame(Name=LocNUC_u$Name)

# ALL RNA POLYMERASES
RNApol<-tair_GO[grep("RNA polymerase", tair_GO$String), ]
#remove duplicated gene names
RNApol_u<-RNApol[!duplicated(RNApol$Name), ]
#Make sure genes are Nuc localized
RNApol_u<-join(RNApol_u, LocNUC_NAME, type="inner")
#Append a "Localization" column onto this file.
RNApol_u$Localization<-"RNApol"
#Join my data with NUC (GO) localization data
My_RNApol<-join(my_data, RNApol_u, type ="left")

#Count AC top genes
AC_RNApol<-nrow(subset(My_RNApol, My_RNApol$Topology_loose == "AC_topology" & My_RNApol$Localization == "RNApol"))
AC_noRNApol<-(nrow(subset(My_RNApol, My_RNApol$Topology_loose == "AC_topology")))-(nrow(subset(My_RNApol, My_RNApol$Topology_loose == "AC_topology" & My_RNApol$Localization == "RNApol")))
#Count BC top genes
BC_RNApol<-nrow(subset(My_RNApol, My_RNApol$Topology_loose == "BC_topology" & My_RNApol$Localization == "RNApol"))
BC_noRNApol<-(nrow(subset(My_RNApol, My_RNApol$Topology_loose == "BC_topology")))-(nrow(subset(My_RNApol, My_RNApol$Topology_loose == "BC_topology" & My_RNApol$Localization == "RNApol")))

#Create twoway table
twoway_matrix_fullRNApol<-matrix(c(AC_RNApol, AC_noRNApol, BC_RNApol, BC_noRNApol), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_fullRNApol) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_fullRNApol) <- c("RNA_pol", "Non_RNA_pol")
twoway_table_fullRNApol<-as.table(twoway_matrix_fullRNApol)

#Percent AC
propAC_RNApol<-AC_RNApol/(AC_RNApol+AC_noRNApol)
#Percent BC
propBC_RNApol<-BC_RNApol/(BC_RNApol+BC_noRNApol)

#Enrichment score
enrich_RNApol<-(propAC_RNApol-propBC_RNApol)/(propAC_RNApol+propBC_RNApol)

#Significance
fisher.test(twoway_matrix_fullRNApol, alternative = "two.sided")
fisher.test(twoway_matrix_fullRNApol, alternative = "less")

# RNA POLYMERASE I
RNApolI<-tair_GO[grep("RNA polymerase I ", tair_GO$String), ]
#remove duplicated gene names
RNApolI_u<-RNApolI[!duplicated(RNApolI$Name), ]
#Make sure genes are Nuc localized
RNApolI_u<-join(RNApolI_u, LocNUC_NAME, type="inner")
#Append a "Localization" column onto this file.
RNApolI_u$Localization<-"RNApolI"

#Join my data with NUC (GO) localization data
My_RNApolI<-join(my_data, RNApolI_u, type ="left")

#Count AC top genes
AC_RNApolI<-nrow(subset(My_RNApolI, My_RNApolI$Topology_loose == "AC_topology" & My_RNApolI$Localization == "RNApolI"))
AC_noRNApolI<-(nrow(subset(My_RNApolI, My_RNApolI$Topology_loose == "AC_topology")))-(nrow(subset(My_RNApolI, My_RNApolI$Topology_loose == "AC_topology" & My_RNApolI$Localization == "RNApolI")))

#Count BC top genes
BC_RNApolI<-nrow(subset(My_RNApolI, My_RNApolI$Topology_loose == "BC_topology" & My_RNApolI$Localization == "RNApolI"))
BC_noRNApolI<-(nrow(subset(My_RNApolI, My_RNApolI$Topology_loose == "BC_topology")))-(nrow(subset(My_RNApolI, My_RNApolI$Topology_loose == "BC_topology" & My_RNApolI$Localization == "RNApolI")))

twoway_matrix_RNApolI<-matrix(c(AC_RNApolI, AC_noRNApolI, BC_RNApolI, BC_noRNApolI), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_RNApolI) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_RNApolI) <- c("RNA_polI", "Non_RNA_polI")
twoway_table_RNApolI<-as.table(twoway_matrix_RNApolI)

#Percent AC
propAC_RNApolI<-AC_RNApolI/(AC_RNApolI+AC_noRNApolI)
#Percent BC
propBC_RNApolI<-BC_RNApolI/(BC_RNApolI+BC_noRNApolI)

#Enrichment score
enrich_RNApolI<-(propAC_RNApolI-propBC_RNApolI)/(propAC_RNApolI+propBC_RNApolI)

#Significance
fisher.test(twoway_matrix_RNApolI, alternative = "two.sided")
fisher.test(twoway_matrix_RNApolI, alternative = "less")

# RNA POLYMERASE II
RNApolII<-tair_GO[grep("RNA polymerase II ", tair_GO$String), ]
#remove duplicated gene names
RNApolII_u<-RNApolII[!duplicated(RNApolII$Name), ]
#Make sure genes are Nuc localized
RNApolII_u<-join(RNApolII_u, LocNUC_NAME, type="inner")
#Append a "Localization" column onto this file.
RNApolII_u$Localization<-"RNApolII"

#Join my data with NUC (GO) localization data
My_RNApolII<-join(my_data, RNApolII_u, type ="left")

#Count AC top genes
AC_RNApolII<-nrow(subset(My_RNApolII, My_RNApolII$Topology_loose == "AC_topology" & My_RNApolII$Localization == "RNApolII"))
AC_noRNApolII<-(nrow(subset(My_RNApolII, My_RNApolII$Topology_loose == "AC_topology")))-(nrow(subset(My_RNApolII, My_RNApolII$Topology_loose == "AC_topology" & My_RNApolII$Localization == "RNApolII")))

#Count BC top genes
BC_RNApolII<-nrow(subset(My_RNApolII, My_RNApolII$Topology_loose == "BC_topology" & My_RNApolII$Localization == "RNApolII"))
BC_noRNApolII<-(nrow(subset(My_RNApolII, My_RNApolII$Topology_loose == "BC_topology")))-(nrow(subset(My_RNApolII, My_RNApolII$Topology_loose == "BC_topology" & My_RNApolII$Localization == "RNApolII")))

twoway_matrix_RNApolII<-matrix(c(AC_RNApolII, AC_noRNApolII, BC_RNApolII, BC_noRNApolII), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_RNApolII) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_RNApolII) <- c("RNApolII", "Non_RNApolII")
twoway_table_RNApolII<-as.table(twoway_matrix_RNApolII)

#Percent AC
propAC_RNApolII<-AC_RNApolII/(AC_RNApolII+AC_noRNApolII)
#Percent BC
propBC_RNApolII<-BC_RNApolII/(BC_RNApolII+BC_noRNApolII)

#Enrichment score
enrich_RNApolII<-(propAC_RNApolII-propBC_RNApolII)/(propAC_RNApolII+propBC_RNApolII)

#Significance
fisher.test(twoway_matrix_RNApolII, alternative = "two.sided")
fisher.test(twoway_matrix_RNApolII, alternative = "less")

# RNA POLYMERASE III
RNApolIII<-tair_GO[grep("RNA polymerase III ", tair_GO$String), ]
#remove duplicated gene names
RNApolIII_u<-RNApolIII[!duplicated(RNApolIII$Name), ]
#Make sure genes are Nuc localized
RNApolIII_u<-join(RNApolIII_u, LocNUC_NAME, type="inner")
#Append a "Localization" column onto this file.
RNApolIII_u$Localization<-"RNApolIII"
#Join my data with NUC (GO) localization data
My_RNApolIII<-join(my_data, RNApolIII_u, type ="left")

#Count AC top genes
AC_RNApolIII<-nrow(subset(My_RNApolIII, My_RNApolIII$Topology_loose == "AC_topology" & My_RNApolIII$Localization == "RNApolIII"))
AC_noRNApolIII<-(nrow(subset(My_RNApolIII, My_RNApolIII$Topology_loose == "AC_topology")))-(nrow(subset(My_RNApolIII, My_RNApolIII$Topology_loose == "AC_topology" & My_RNApolIII$Localization == "RNApolIII")))

#Count BC top genes
BC_RNApolIII<-nrow(subset(My_RNApolIII, My_RNApolIII$Topology_loose == "BC_topology" & My_RNApolIII$Localization == "RNApolIII"))
BC_noRNApolIII<-(nrow(subset(My_RNApolIII, My_RNApolIII$Topology_loose == "BC_topology")))-(nrow(subset(My_RNApolIII, My_RNApolIII$Topology_loose == "BC_topology" & My_RNApolIII$Localization == "RNApolIII")))


twoway_matrix_RNApolIII<-matrix(c(AC_RNApolIII, AC_noRNApolIII, BC_RNApolIII, BC_noRNApolIII), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_RNApolIII) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_RNApolIII) <- c("RNApolIII", "Non_RNApolIII")
twoway_table_RNApolIII<-as.table(twoway_matrix_RNApolIII)

#Percent AC
propAC_RNApolIII<-AC_RNApolIII/(AC_RNApolIII+AC_noRNApolIII)
#Percent BC
propBC_RNApolIII<-BC_RNApolIII/(BC_RNApolIII+BC_noRNApolIII)

#Enrichment score
enrich_RNApolIII<-(propAC_RNApolIII-propBC_RNApolIII)/(propAC_RNApolIII+propBC_RNApolIII)

#Significance
fisher.test(twoway_matrix_RNApolIII, alternative = "two.sided")
fisher.test(twoway_matrix_RNApolIII, alternative = "greater")


# RNA POLYMERASE IV
RNApolIV<-tair_GO[grep("RNA polymerase IV ", tair_GO$String), ]
#remove duplicated gene names
RNApolIV_u<-RNApolIV[!duplicated(RNApolIV$Name), ]
#Make sure genes are Nuc localized
RNApolIV_u<-join(RNApolIV_u, LocNUC_NAME, type="inner")
#Append a "Localization" column onto this file.
RNApolIV_u$Localization<-"RNApolIV"
#Join my data with NUC (GO) localization data
My_RNApolIV<-join(my_data, RNApolIV_u, type ="left")

#Count AC top genes
AC_RNApolIV<-nrow(subset(My_RNApolIV, My_RNApolIV$Topology_loose == "AC_topology" & My_RNApolIV$Localization == "RNApolIV"))
AC_noRNApolIV<-(nrow(subset(My_RNApolIV, My_RNApolIV$Topology_loose == "AC_topology")))-(nrow(subset(My_RNApolIV, My_RNApolIV$Topology_loose == "AC_topology" & My_RNApolIV$Localization == "RNApolIV")))

#Count BC top genes
BC_RNApolIV<-nrow(subset(My_RNApolIV, My_RNApolIV$Topology_loose == "BC_topology" & My_RNApolIV$Localization == "RNApolIV"))
BC_noRNApolIV<-(nrow(subset(My_RNApolIV, My_RNApolIV$Topology_loose == "BC_topology")))-(nrow(subset(My_RNApolIV, My_RNApolIV$Topology_loose == "BC_topology" & My_RNApolIV$Localization == "RNApolIV")))

twoway_matrix_RNApolIV<-matrix(c(AC_RNApolIV, AC_noRNApolIV, BC_RNApolIV, BC_noRNApolIV), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_RNApolIV) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_RNApolIV) <- c("RNApolIV", "Non_RNApolIV")
twoway_table_RNApolIV<-as.table(twoway_matrix_RNApolIV)

#Percent AC
propAC_RNApolIV<-AC_RNApolIV/(AC_RNApolIV+AC_noRNApolIV)
#Percent BC
propBC_RNApolIV<-BC_RNApolIV/(BC_RNApolIV+BC_noRNApolIV)

#Enrichment score
enrich_RNApolIV<-(propAC_RNApolIV-propBC_RNApolIV)/(propAC_RNApolIV+propBC_RNApolIV)

#Significance
fisher.test(twoway_matrix_RNApolIV, alternative = "two.sided")
fisher.test(twoway_matrix_RNApolIV, alternative = "less")


# RNA POLYMERASE V
RNApolV<-tair_GO[grep("RNA polymerase V ", tair_GO$String), ]
#remove duplicated gene names
RNApolV_u<-RNApolV[!duplicated(RNApolV$Name), ]
#Make sure genes are Nuc localized
RNApolV_u<-join(RNApolV_u, LocNUC_NAME, type="inner")
#Append a "Localization" column onto this file.
RNApolV_u$Localization<-"RNApolV"

#Join my data with NUC (GO) localization data
My_RNApolV<-join(my_data, RNApolV_u, type ="left")

#Count AC top genes
AC_RNApolV<-nrow(subset(My_RNApolV, My_RNApolV$Topology_loose == "AC_topology" & My_RNApolV$Localization == "RNApolV"))
AC_noRNApolV<-(nrow(subset(My_RNApolV, My_RNApolV$Topology_loose == "AC_topology")))-(nrow(subset(My_RNApolV, My_RNApolV$Topology_loose == "AC_topology" & My_RNApolV$Localization == "RNApolV")))

#Count BC top genes
BC_RNApolV<-nrow(subset(My_RNApolV, My_RNApolV$Topology_loose == "BC_topology" & My_RNApolV$Localization == "RNApolV"))
BC_noRNApolV<-(nrow(subset(My_RNApolV, My_RNApolV$Topology_loose == "BC_topology")))-(nrow(subset(My_RNApolV, My_RNApolV$Topology_loose == "BC_topology" & My_RNApolV$Localization == "RNApolV")))

twoway_matrix_RNApolV<-matrix(c(AC_RNApolV, AC_noRNApolV, BC_RNApolV, BC_noRNApolV), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_RNApolV) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_RNApolV) <- c("RNApolV", "Non_RNApolV")
twoway_table_RNApolV<-as.table(twoway_matrix_RNApolV)

#Percent AC
propAC_RNApolV<-AC_RNApolV/(AC_RNApolV+AC_noRNApolV)
#Percent BC
propBC_RNApolV<-BC_RNApolV/(BC_RNApolV+BC_noRNApolV)

#Enrichment score
enrich_RNApolV<-(propAC_RNApolV-propBC_RNApolV)/(propAC_RNApolV+propBC_RNApolV)

#Significance
fisher.test(twoway_matrix_RNApolV, alternative = "two.sided")
fisher.test(twoway_matrix_RNApolV, alternative = "less")


# TELOMERE COMPLEX
telo<-tair_GO[grep("telomere", tair_GO$String), ]
#remove duplicated gene names
telo_u<-telo[!duplicated(telo$Name), ]
#Make sure genes are Nuc localized
telo_u<-join(telo_u, LocNUC_NAME, type="inner")
#Append a "Localization" column onto this file.
telo_u$Localization<-"telo"
#Join my data with NUC (GO) localization data
My_telo<-join(my_data, telo_u, type ="left")

#Count AC top genes
AC_telo<-nrow(subset(My_telo, My_telo$Topology_loose == "AC_topology" & My_telo$Localization == "telo"))
AC_notelo<-(nrow(subset(My_telo, My_telo$Topology_loose == "AC_topology")))-(nrow(subset(My_telo, My_telo$Topology_loose == "AC_topology" & My_telo$Localization == "telo")))

#Count BC top genes
BC_telo<-nrow(subset(My_telo, My_telo$Topology_loose == "BC_topology" & My_telo$Localization == "telo"))
BC_notelo<-(nrow(subset(My_telo, My_telo$Topology_loose == "BC_topology")))-(nrow(subset(My_telo, My_telo$Topology_loose == "BC_topology" & My_telo$Localization == "telo")))

twoway_matrix_telo<-matrix(c(AC_telo, AC_notelo, BC_telo, BC_notelo), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_telo) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_telo) <- c("Nuc_localized", "Non_Nuc_localized")
twoway_table_telo<-as.table(twoway_matrix_telo)

#Percent AC
propAC_telo<-AC_telo/(AC_telo+AC_notelo)
#Percent BC
propBC_telo<-BC_telo/(BC_telo+BC_notelo)

#Enrichment score
enrich_telo<-(propAC_telo-propBC_telo)/(propAC_telo+propBC_telo)

#Significance
fisher.test(twoway_matrix_telo, alternative = "two.sided")
fisher.test(twoway_matrix_telo, alternative = "greater")

tair_GO$String[grep("cyclin-dependent", tair_GO$String, ignore.case = TRUE)]

# CYCLIN-DEPENDANT KINASE COMPLEX
cyclin<-tair_GO[grep("cyclin-dependent", tair_GO$String), ]
#remove duplicated gene names
cyclin_u<-cyclin[!duplicated(cyclin$Name), ]

#Append a "Localization" column onto this file.
cyclin_u$Localization<-"cyclin"

#Join my data with NUC (GO) localization data
My_cyclin<-join(my_data, cyclin_u, type ="left")

#Count AC top genes
AC_cyclin<-nrow(subset(My_cyclin, My_cyclin$Topology_loose == "AC_topology" & My_cyclin$Localization == "cyclin"))
AC_nocyclin<-(nrow(subset(My_cyclin, My_cyclin$Topology_loose == "AC_topology")))-(nrow(subset(My_cyclin, My_cyclin$Topology_loose == "AC_topology" & My_cyclin$Localization == "cyclin")))

#Count BC top genes
BC_cyclin<-nrow(subset(My_cyclin, My_cyclin$Topology_loose == "BC_topology" & My_cyclin$Localization == "cyclin"))
BC_nocyclin<-(nrow(subset(My_cyclin, My_cyclin$Topology_loose == "BC_topology")))-(nrow(subset(My_cyclin, My_cyclin$Topology_loose == "BC_topology" & My_cyclin$Localization == "cyclin")))

twoway_matrix_cyclin<-matrix(c(AC_cyclin, AC_nocyclin, BC_cyclin, BC_nocyclin), ncol = 2, byrow = FALSE)
colnames(twoway_matrix_cyclin) <- c("AC_topology", "BC_topology")
rownames(twoway_matrix_cyclin) <- c("CDK", "Non_CDK")
twoway_table_cyclin<-as.table(twoway_matrix_cyclin)

#Percent AC
propAC_cyclin<-AC_cyclin/(AC_cyclin+AC_nocyclin)
#Percent BC
propBC_cyclin<-BC_cyclin/(BC_cyclin+BC_nocyclin)

#Enrichment score
enrich_cyclin<-(propAC_cyclin-propBC_cyclin)/(propAC_cyclin+propBC_cyclin)

#Significance
fisher.test(twoway_matrix_cyclin, alternative = "two.sided")
fisher.test(twoway_matrix_cyclin, alternative = "greater")


### MAKE TABLE OF ALL THE RESULTS  ###
#results_matrix<-matrix(c(enrich_CP,enrich_photosynth,enrich_light,
#         enrich_dark,enrich_photoresp,enrich_MT,
#         enrich_resp,enrich_cyt,enrich_glyc,enrich_etrans,
#         enrich_PPR,enrich_NUC,enrich_RNApol,enrich_RNApolI,
#         enrich_RNApolII,enrich_RNApolIII,enrich_RNApolIV,
#         enrich_RNApolV,enrich_telo,enrich_cyclin), ncol=1, byrow = FALSE)


#enrich_df<-data.frame(Enrich=results_matrix)

all_results<-matrix(c(AC_CP,propAC_CP,BC_CP,propBC_CP,enrich_CP,0.007756,0.004429,
                      AC_photosynth,propAC_photosynth,BC_photosynth,propBC_photosynth,enrich_photosynth,0.01419,0.01184,
                      AC_light,propAC_light,BC_light,propBC_light,enrich_light,0.005326,0.005326,
                      AC_dark,propAC_dark,BC_dark,propBC_dark,enrich_dark,0.04469,0.04469,
                      AC_photoresp,propAC_photoresp,BC_photoresp,propBC_photoresp,enrich_photoresp,0.4572,0.4572,
                      AC_MT,propAC_MT,BC_MT,propBC_MT,enrich_MT,0.004364,0.002504,
                      AC_resp,propAC_resp,BC_resp,propBC_resp,enrich_resp,0.4475,0.3246,
                      AC_cyt,propAC_cyt,BC_cyt,propBC_cyt,enrich_cyt,1,0.7054,
                      AC_glyc,propAC_glyc,BC_glyc,propBC_glyc,enrich_glyc,0.2329,0.2329,
                      AC_etrans,propAC_etrans,BC_etrans,propBC_etrans,enrich_etrans,0.5441,0.5441,
                      AC_PPR,propAC_PPR,BC_PPR,propBC_PPR,enrich_PPR,0.02507,0.01682,
                      AC_NUC,propAC_NUC,BC_NUC,propBC_NUC,enrich_NUC,0.01786,0.009361,
                      AC_RNApol,propAC_RNApol,BC_RNApol,propBC_RNApol,enrich_RNApol,0.1478,0.08138,
                      AC_RNApolI,propAC_RNApolI,BC_RNApolI,propBC_RNApolI,enrich_RNApolI,1,0.5924,
                      AC_RNApolII,propAC_RNApolII,BC_RNApolII,propBC_RNApolII,enrich_RNApolII,0.1274,0.05929,
                      AC_RNApolIII,propAC_RNApolIII,BC_RNApolIII,propBC_RNApolIII,enrich_RNApolIII,0.5823,0.5823,
                      AC_RNApolIV,propAC_RNApolIV,BC_RNApolIV,propBC_RNApolIV,enrich_RNApolIV,1,0.7054,
                      AC_RNApolV,propAC_RNApolV,BC_RNApolV,propBC_RNApolV,enrich_RNApolV,1,0.7054,
                      AC_telo,propAC_telo,BC_telo,propBC_telo,enrich_telo,0.1061,0.1061,
                      AC_cyclin,propAC_cyclin,BC_cyclin,propBC_cyclin,enrich_cyclin,0.2322,0.1928), 
                      ncol=7, byrow = TRUE)

results_df<-data.frame(AC_count=all_results[,1], AC_prop=all_results[,2], BC_count=all_results[,3], BC_prop=all_results[,4], Enrich_score=all_results[,5], two_sided=all_results[,6], one_sided=all_results[,7])

