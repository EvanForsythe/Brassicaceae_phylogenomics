# Brassicaceae_phylogenomics

The repository contains scripts used to infer gene trees from genes spanning nuclear and organellar genomes. Also included are scripts for performing various downstream analyses of gene tree topologies. 

Examples input files are located in R_test_dir/

Scripts are written in perl, shell, and R.

List of scripts included:

Clean_CDSandAAseqs.R
-An R script that cleans CDS files to contain only the longest isoform for each gene. This script also translates CDS to AA and output an AA fasta file to ensure that CDS and AA files are compatible.

batch_raxml.pl
-A perl wrapper script for running for inferring gene trees from multiple alignments 

Topology_analysis_nuc.R
-Script for assessing the topology of single-copy nuclear gene trees

Topology_analysis_ret_dups.R
-Script for splitting retained duplicate trees and assessing the topology of subtrees

Topology_analysis_CP.R
-Script for assessing the topology of CP trees

Topology_analysis_MT.R
-Script for assessing the topology of MT trees

Print_exact_tops.R
-Prints gene tree topologies with branchlengths and support values removed to make them more human readable        

D_stat.R
-script used to calculate stie-pattern based (standard application) D-statistics.

Dgt_stat_BCvsAB.R
-Script used to calculate the Dgt-statistic comparing the frequencies of A(BC) and C(AB) for nuclear genes

Dgt_stat_ACvsAB.R
-Script used to calculate the Dgt-statistic comparing the frequencies of B(AC) and C(AB) for nuclear genes

Phylonet.nex
-Nexus file containing rooted nuclear gene trees. The last line of the file is the executable line for PhyloNet (this example is to run PhlyoNet with 0 reticulations). Analysis is run from command line with: java -jar /pathtophylonet/PhyloNet_3.6.1.jar Phylonet.nex >output_file

TICR.R 
-Script for performing the TICR test. 

Node_depths.R
-Rate smooth gene trees and calculate node depths

IG_directionality.R
-Uses dS distances to perform a test of IG directionality (i.e. determine donor and recipient lineages)

GO_categories.R
-Extracts GO category data and calculated proportion of A(BC) and B(AC) with each GO category        
       
Prot_Interact_Network.R 
-Script used to perform a network analysis using protein-protein interaction data

Chromosome_maps.R
-Script to generate chromosome map figures

Linkage_blocks.R
-Script used to detect linear clustering of genes with the same topology along the C. rubella chromosome map.

Window_based_node_depth.R
-Script used for window-based node depth analysis

