#Set the working directory
setwd("./R_test_dir/")

#Check out Blanc and Wolfe, 2004

#Load igraph
library(igraph)

#load seqinr
#install.packages("seqinr")
#library(seqinr)

#Read in the datatable
data_wtops<-read.csv("interactome_mydata170421.csv", header=TRUE)

#Get only the partner names
partners_df<-data.frame(PartA=data_wtops$Partner_A, PartB=data_wtops$Partner_B)

#make a network graph
network_graph<-graph.data.frame(partners_df, directed = FALSE)

#plot the graph
plot.igraph(network_graph, vertex.label=NA, vertex.size=1)

#make an adjecency matrix (working?)
#network_adj_mat<-as_adjacency_matrix(network_graph)

V(network_graph)$name

#import the sample_attributes
a<-read.csv("node_attributes170422.csv")

a<-unique(a)

#assign Topology attribute
V(network_graph)$Top=as.character(a$Topology[match(V(network_graph)$name,a$Partner)]) # This code says to create a vertex attribute called "Sex" by extracting the value of the column "Sex" in the attributes file when the Bird ID number matches the vertex name.
V(network_graph)$Top # This will print the new vertex attribute, "Top"

#assign BS attribute
V(network_graph)$BS=as.numeric(a$BS[match(V(network_graph)$name,a$Partner)]) # This code says to create a vertex attribute called "Sex" by extracting the value of the column "Sex" in the attributes file when the Bird ID number matches the vertex name.
V(network_graph)$BS # This will print the new vertex attribute, "BS"

#Assign colors to the topologies
V(network_graph)$color=V(network_graph)$Top
V(network_graph)$color=gsub("BC", "red", V(network_graph)$color)
V(network_graph)$color=gsub("AC", "green", V(network_graph)$color)
V(network_graph)$color=gsub("AB", "blue", V(network_graph)$color)

#Assign sizes to the BS scores
V(network_graph)$size=(as.numeric(V(network_graph)$BS)/100)*5
V(network_graph)$size

#Plot the network
plot.igraph(network_graph,
            #vertex.label.cex=0.4
            vertex.label=NA)

#Make an interactive viewer
install.packages('visNetwork')
library('visNetwork') 

#Use visigraph to generate an interactive map
visIgraph(network_graph)

####ASSORTATIVITY ######
#Assign number code for the topologies
#These codes are needed to run assortativity 
V(network_graph)$Top_code=V(network_graph)$Top
V(network_graph)$Top_code=gsub("BC", "1", V(network_graph)$Top_code)
V(network_graph)$Top_code=gsub("AC", "2", V(network_graph)$Top_code)
V(network_graph)$Top_code=gsub("AB", "3", V(network_graph)$Top_code)

#Nominal assortativity coefficient calculation
obs_assort<-assortativity.nominal(network_graph, types = V(network_graph)$Top_code, directed = FALSE)

###create randomized attributes on the same network###

#make function
rando_assort<-function(){
#Assign number code for the topologies
#These codes are needed to run assortativity 
V(network_graph)$Top_code_rando=sample(V(network_graph)$Top, replace = FALSE)
V(network_graph)$Top_code_rando=gsub("BC", "1", V(network_graph)$Top_code_rando)
V(network_graph)$Top_code_rando=gsub("AC", "2", V(network_graph)$Top_code_rando)
V(network_graph)$Top_code_rando=gsub("AB", "3", V(network_graph)$Top_code_rando)
#Nominal assortativity coefficient calculation
assortativity.nominal(network_graph, types = V(network_graph)$Top_code_rando, directed = FALSE)
}
#run the function
assort_reps<-replicate(10000, rando_assort())

#plot the curve of the null distribution
plot(density(assort_reps))
#plot a line where the observed Assortativity is
abline(v=obs_assort)

#perform a Z-test of observed against the null
a<-obs_assort
s<-sd(assort_reps)
n<-length(assort_reps)
xbar<-mean(assort_reps)
z<-(xbar-a)/s

#Two tailed p-value
2*pnorm(-abs(z))

