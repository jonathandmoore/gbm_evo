
	
#### This part parses command line options or sets defaults

if(!require(optparse)){
	install.packages("optparse")
	library(optparse)
}

library("optparse")

option_list = list(
	make_option(c("-p", "--project"), action="store", default=NA, type='character',
              help="project code for location of bs_sequel output"),
	make_option(c("-c", "--context"), action="store", default=NA, type='character',
              help="methylation context (CG, CHG or CHH)"),
	make_option(c("-a", "--action"), action="store", default="all", type='character',
              help="actions to perform (load, plot, call, merge, analyse or all)"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]"),
	make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Make the program not be verbose.")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$verbose) {
    # show the user which options were chosen
	cat("Project: ")
    cat(opt$project)
    cat("\n")
    
    cat("Context: ")
    cat(opt$context)
    cat("\n")

    cat("Action: ")
    cat(opt$action)
    cat("\n")
}

# project_id and meth_context together are used for finding raw data and naming output files
project_id=NULL
if (is.na(opt$project)) {
	# No project defined - set a default
	project_id="test"
	#project_id="SRA035939"  # Schmitz et al, 2011
	#project_id="PRJEB2678"  # Becker et al, 2011
} else {
	project_id=opt$project
}
	
meth_context=NULL
if (is.na(opt$context)) {
	# No context defined - set a default
	meth_context="CG"
	#meth_context="CHG"cd 
	#meth_context="CHH"
} else {
	meth_context=opt$context
}

action=NULL
if (is.na(opt$action)) {
	# No action defined - set a default
	action="all"
	#action="load"
	#action="plot"
	#action="call"
	#action="merge"
	#action="analyse"
} else {
	action=opt$action
}

# Platform-specific stuff:
#  Base of path to project
#  Where to find bioconductor (packages need to be pre-installed in R using singularity if running on the cluster)
on_cluster = FALSE
pathroot = ""
if (.Platform$OS.type=="windows") {
  pathroot="D:"
  source("http://bioconductor.org/biocLite.R")
} else if (.Platform$OS.type=="unix") {
  if(substr(getwd(),1,20)=="/jic/scratch/groups/") {
    # assume we are running on cluster
    pathroot="/jic/scratch/groups/Daniel-Zilberman"
    on_cluster = TRUE
  } else {
    # assume we are running in Virtualbox with shared folder /media/sf_D_DRIVE
    pathroot="/media/sf_D_DRIVE"
    source("http://bioconductor.org/biocLite.R")
  }
}

if (!on_cluster) {
  # Assume that if we are on the cluster, then we are running in a VM with all relevant packages pre-installed
  BiocManager::install(c("GenomicRanges", "BSgenome", "BSgenome.Athaliana.TAIR.TAIR9", "MethylSeekR", "karyoploteR", "ggtree"))
}


if(!require(data.table)){
	install.packages("data.table")
	library(data.table)
}
if(!require(stringr)){
	install.packages("stringr")
	library(stringr)
}
if(!require(ggplot2)){
	install.packages("ggplot2")
	library(ggplot2)
}
if(!require(mixtools)){
	install.packages("mixtools")
	library(mixtools)
}
if(!require(rootSolve)){
	install.packages("rootSolve")
	library(rootSolve)
}
if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr)
}

	
#### This part sets up libraries and contextual information and metadata about the project

# Source directory containing alignments from bs_sequel
source_dir = paste0(pathroot,"/Projects/Ancestral_Methylation_Analysis/")

# Working directory where outputs will be directed to
setwd(dir = paste0(pathroot,"/Projects/Ancestral_Methylation_Analysis/"))

reference_fasta = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta")  # TAIR10 genome assembly This never gets used at the moment
reference_CG_sites = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.CG_sites.tsv")   # We prepared this earlier using a perl script to parse the TAIR10 genome assembly
reference_gff = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff")  # AraPort11 annotation
reference_exons = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-exon.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_introns = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-intron.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_5UTR = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-five_prime_UTR.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_3UTR = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-three_prime_UTR.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
### SNP_loci.txt - we need to import a list of all known SNPs between lines considered for analysis, so SNP loci can inform site exclusion criteria
reference_chromatin_states = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Chromatin_states/Sequeira-Mendes_et_al_2014_tpc124578_SupplementalIDS2.txt")  # genomic ranges assigning each segment of the nuclear genome to one of 9 chromatin states
reference_DHS_loci = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/DHS_scores/TAIR10_DHSs.gff")  # DNase-hypersensitivity sites (indicative of open chromatin)

# Nucleosomes position data from http://plantdhs.org/Download (Jiming Jiang lab) Originally described in Wu Y.F. Zhang W.L. Jiang J.M. Genome-wide nucleosome positioning is orchestrated by genomic regions associated with DNase I hypersensitivity in rice PLoS Genet.  2014 10 e1004378
# For nucleosome positioning data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.RDS")
#reference_nucleosomes = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.RDS")  

# Nucleosomes position data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2807196
# Lyons DB, Zilberman D. DDM1 and Lsh remodelers allow methylation of DNA wrapped in nucleosomes. Elife 2017 Nov 15;6. PMID: 29140247
# GSE96994 wt.rep1.mnase_seq results
# GSM2807196_wt.rep1.nucleosomes.bed is all nucleosome calls
# GSM2807196_wt.group_1.bed is well-positioned nucleosome calls

#reference_nucleosomes = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Nucleosomes/GSM2807196_wt.rep1.nucleosomes.bed")
reference_nucleosomes = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Nucleosomes/GSM2807196_wt.group_1.bed")

# H3K9me2 normalised WT data from Dave Lyons
# For H3K9me2 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/dk9.jacobsen.normalized.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/dk9.jacobsen.normalized.RDS")
reference_H3K9me2 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/dk9.jacobsen.normalized.RDS")  

# H3K9me1 Plant DNase I hypersensitivity data. Downloaded from http://plantdhs.org/Download (Jiming Jiang lab) on 05/12/2017
# For H3K4me1 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.RDS")
reference_H3K4me1 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.RDS")

# H3K9me1 Plant DNase I hypersensitivity data. Downloaded from http://plantdhs.org/Download (Jiming Jiang lab) on 05/12/2017
# For H3K27me3 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.RDS")
reference_H3K27me3 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.RDS")

# H2AW data from Jaemyung
reference_H2AW = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/h2aw.w50.gff")


# get a list of the samples from the subdirectories of the source directory
samples = list.dirs(path = source_dir, full.names = FALSE, recursive = FALSE)
no_samples=length(samples)

# Read the sample metadata in from a text file in the source directory
# Commented out this line as we don't use sample metadata in this workflow right now
#sample_metadata=read.delim(paste0(source_dir,project_id,"_metadata.tsv"), header=TRUE, sep="\t")
#rownames(sample_metadata)=sample_metadata$Identifier

# Read the Schmitz et al, 2011 DMRs data, so this can be used in masking where needed
Schmitz_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_CG_DMRs.txt")
Schmitz_non_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_non-CG_DMRs.txt")


mixedCaseChr <- function(s, strict = FALSE) {
	paste0(toupper(substr(s,1,1)),tolower(substr(s,2,3)),toupper(substr(s,4,4)))
}

# Read in the matrix of meth status values
# genes only:
#tree_samples_meth_status <- readRDS(paste0("4-summary/AncestralAnalysis_CG_gene_CG_genmeth.rds"))

# all sites in genome in 99 samples from literature:
tree_samples_meth_status <- readRDS(paste0("4-summary/AncestralAnalysis_CG_all_samples_meth_status.rds"))

# all sites in genome in runs including root taxa:
tree_samples_meth_status2 <- readRDS(paste0("4-summary/SeqRunGenAll_CG_all_samples_meth_status.rds"))

# Merge the three sets of 'root taxa' data with the other samples
tree_samples_meth_status = cbind(tree_samples_meth_status, tree_samples_meth_status2[,4:6])
rm(tree_samples_meth_status2)

# previously calculated genome-wide set:
#tree_samples_meth_status = valid_samples_meth_status

tree_samples_meth_status$Chromosome = mixedCaseChr(tree_samples_meth_status$Chromosome)

# Read in the segmentation model in case we need it later
segmentation_model.gr = readRDS(paste0("4-summary/SRA035939_CG_segmentation_model_draft3.rds"))
levels(segmentation_model.gr@seqnames@values) = mixedCaseChr(as.character(levels(segmentation_model.gr@seqnames@values)))
segmentation_model.gr@seqinfo@seqnames = levels(segmentation_model.gr@seqnames@values)

# Generate a masked segmentation model
segment_mask_loci.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="Chr2", "start_site"=3239693, "end_site"=3505260))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")		
masked_segmentation_model.gr = segmentation_model.gr[unique(queryHits(findOverlaps(segmentation_model.gr, setdiff(segmentation_model.gr,segment_mask_loci.gr))))]

# Read in the table of genes with their characteristics
gff.genes =  readRDS(paste0("4-summary/SRA035939_CG_gff.genes.rds"))

# Make a set of genomic ranges corresponding to GBM genes
GBM_gene_loci=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5","gene_ID")]
colnames(GBM_gene_loci)=c("Chromosome", "start_site", "end_site", "gene_ID")
GBM_gene_loci$Chromosome=mixedCaseChr(GBM_gene_loci$Chromosome)
GBM_gene_loci.gr=reduce(makeGRangesFromDataFrame(df=GBM_gene_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome"))

# Make a set of genomic ranges corresponding to the sites in the matrix of meth status values
meth_status_sites.gr = makeGRangesFromDataFrame(df=tree_samples_meth_status, start.field="Locus", end.field="Locus", seqnames.field="Chromosome")

# Find overlaps between sites in the matrix of meth status values and GBM genes
#GBM_sites_hits=findOverlaps(GBM_gene_loci.gr, meth_status_sites.gr)

# Extract a data matrix corresponding to the GBM gene loci
tree_samples_GBM_gene_meth_status = tree_samples_meth_status[subjectHits(findOverlaps(GBM_gene_loci.gr, meth_status_sites.gr)),]


# Transform the matrix, and convert character values to numeric values (they represent coordinates in a k dimensional space where k = number of CG sites)
# This version does all genomic loci:
#tree_samples_data_matrix = t(tree_samples_meth_status[,seq(from=4, to=ncol(tree_samples_meth_status))])

# tree_samples_meth_status is finished with now.  It is a large object so delete it for space
#rm(tree_samples_meth_status)
#gc()

# This version does GBM gene loci only:
tree_samples_data_matrix = t(tree_samples_GBM_gene_meth_status[,seq(from=4, to=ncol(tree_samples_GBM_gene_meth_status))])

rm(tree_samples_GBM_gene_meth_status)
gc()

# Basic scheme representing unmethylated as 0, methylated as 1, and all else as missing data  -- later version does it better
#tree_samples_data_matrix[tree_samples_data_matrix == "M"] = 1
#tree_samples_data_matrix[tree_samples_data_matrix == "U"] = 0
#tree_samples_data_matrix[tree_samples_data_matrix %in% c("D", "I", "O", "P")] = NA

 
# Sanity check numbers of sites with each possible balance of M vs U across all samples

# Put the data matrix into a data frame with numerics  NB THIS STEP FAILS IF THERE ARE ONLY "1" VALUES AT A LOCUS - assigns them as 1s instead of 2s
#z1=as.matrix(as.data.frame(lapply(as.data.frame(tree_samples_data_matrix), as.numeric)))
# Convert factor IDs to numbers
#z1[z1==1]=0
#z1[z1==2]=1

# THIS STEP FIXES THE PROBLEM WITH THE APPROACH ABOVE but is a bit slow and uses an enormous amount of RAM for large data sets
#z1=data.frame(lapply(tree_samples_data_matrix, function(x) as.numeric(as.character(x))), check.names=F, row.names = rownames(tree_samples_data_matrix))

# This version does it neatly in one go rather than converting the original matrix first
z1=matrix(,nrow=nrow(tree_samples_data_matrix), ncol=ncol(tree_samples_data_matrix))
 
z1[tree_samples_data_matrix == "M"] = 1
z1[tree_samples_data_matrix == "P"] = 1
z1[tree_samples_data_matrix == "U"] = 0
z1[tree_samples_data_matrix %in% c("D", "I", "O")] = NA
#z1[tree_samples_data_matrix %in% c("D", "I", "O", "P")] = NA

# Alternative scheme to take account of partially methylated sites specifically:
#z1[tree_samples_data_matrix == "M"] = 1
#z1[tree_samples_data_matrix == "U"] = 0
#z1[tree_samples_data_matrix == "P"] = 0.5
#z1[tree_samples_data_matrix %in% c("D", "I", "O")] = NA
 
z1=data.frame(z1, check.names=F, row.names = rownames(tree_samples_data_matrix))
 

# Find the number of lines with methylation, by site (sum up the 1==Methylated columns)
z2=colSums(z1, na.rm=TRUE) # na.rm=TRUE allows to include sites where some lines have NA call.  This reduces the count at these sites somewhat, but otherwise, very few sites would be included

ggplot() + geom_histogram(data=as.data.frame(z2), aes(x=z2), binwidth=1)
# Leave out zero counts 
ggplot() + geom_histogram(data=as.data.frame(z2[z2>0]), aes(x=z2[z2>0]), binwidth=1)


# Calculate a distance matrix between samples 
tree_samples_dist_matrix = dist(z1, method="euclidean")

# Plot a hierarchical clustering of the samples based on the distance matrix
plot(hclust(tree_samples_dist_matrix))

# Extract sorted vector of k nearest neighbours of a given sample:

# Make a variable with no of samples
no_samples = attr(tree_samples_dist_matrix, "Size")
# Make another variable with the number of nearest neighbours we are interested in
knn_samples = no_samples  # set k to be all samples, in case some of the nearer ones aren't the ones we want, then we potentially have all the neighbours in case we need more

# Make empty matrices to hold the nearest neighbours and their distances, for each member of the distance matrix
knn_matrix = matrix(0, ncol = knn_samples, nrow = no_samples)
knd_matrix = knn_matrix

# Go through the samples in their order in the distance matrix, and populate the n x k matrices with one row per sample, one column per heighbour, in order of increasing distance
for(sample_no in 1:no_samples) {
	# Populate the matrices of identifiers and distances from all other samples to this sample, sorted by distance
	# First put the list of sample numbers into the relevant row of the knn_matrix, sorted by distance (the first column is always going to contain the sample itself, which has distance 0)
	knn_matrix[sample_no,] = order(as.matrix(tree_samples_dist_matrix)[sample_no,])[1:knn_samples]
	# Now put the corresponding distances into the other matrix, in the same order, so we have them in case we want to use them
	knd_matrix[sample_no,] = as.matrix(tree_samples_dist_matrix)[sample_no,as.numeric(knn_matrix[sample_no,])]
	# Finally, replace the sample numbers with the sample names
	knn_matrix[sample_no,] = labels(tree_samples_dist_matrix)[as.numeric(knn_matrix[sample_no,])]
}
#colnames(knn_matrix) = labels(tree_samples_dist_matrix)
rownames(knn_matrix) = labels(tree_samples_dist_matrix)
#colnames(knd_matrix) = labels(tree_samples_dist_matrix)
rownames(knd_matrix) = labels(tree_samples_dist_matrix)


# Compress the variation in the distance matrix down to 2 dimensions (for MultiDimensional Scaling plot)
fit = cmdscale( tree_samples_dist_matrix, eig=TRUE, k=2)

# MDS plot with labels for each sample
plot(fit$points[,1], fit$points[,2])
text(fit$points[,1], fit$points[,2], labels=row.names(tree_samples_data_matrix), cex=0.7)


# UMAP package allows for manifold approximation and projection to allow new samples to be added to an existing 2D projection of a collection of samples 
# In other words, it's a bit like a PCA plot, but it lets you add new samples to a pre-existing plot and is more flexible
install.packages("umap")
library(umap)
 
# Make UMAP of data matrix, but need to convert to numeric first. The UMAP algoritm seems to fail when data matrix contains NAs - replace with 0.5?
# Also, runs out of space trying to make a 77GB vector
#this.umap = umap(matrix(as.numeric(unlist(tree_samples_data_matrix)), nrow=nrow(tree_samples_data_matrix)))
#this.umap = umap(matrix(as.numeric(unlist(z1)), nrow=nrow(z1)))
this.umap = umap(matrix(as.numeric(unlist(ifelse(is.na(z1),0.5,z1))), nrow=nrow(z1)))
plot(this.map$layout)
# Not sure if the labels are in the right order for this part - check!
text(this.map$layout[,1], this.map$layout[,2], labels=row.names(z1), cex=0.7)



# Write data matrix in Nexus-like format for tree building in MrBayes
write.table(tree_samples_data_matrix, file="5-Analysis/tree_samples.txt", quote=FALSE, sep="", col.names=FALSE, na="?", row.names=paste0("  ",row.names(tree_samples_data_matrix)," "))
cat(paste0("#NEXUS\nBegin Data;\nDimensions Ntax=",nrow(tree_samples_data_matrix)," Nchar=",ncol(tree_samples_data_matrix),";\nFormat Datatype=Restriction missing=?;\nMatrix\n"))

# Write data matrix in FASTA-like format for tree building in RAxML
#write.table(tree_samples_data_matrix, file="4-Summary/tree_samples.fas", quote=FALSE, sep="", col.names=FALSE, na="?", row.names=paste0(">",row.names(tree_samples_data_matrix),"\n"))
write.table(z1, file="4-Summary/tree_samples_P.fas", quote=FALSE, sep="", col.names=FALSE, na="?", row.names=paste0(">",row.names(z1),"\n"))

# At this point we can use MrBayes or RAxML to make a tree.  See e.g. script 5-RAxML_2018-07-30.bat



# Doing phylogenetic analyses in R:

## install CRAN Task View for phylogenetics
install.packages('ctv')
library('ctv')
install.views('Phylogenetics')
update.views('Phylogenetics')

install.packages("phylotate")


# Make a data matrix suitable for phylogenetic inference in R
#this_data = tree_samples_data_matrix 
this_data = z1
this_data[is.na(this_data)] = "?"
# THIS TAKES AGES!


# let's try ape
#library(ape)

#library(phylotate)
#this_tree = read_annotated("5-Analysis/tree_samples.nex.con.tre", format="nexus")

# some ideas:
#this_tree= read_annotated("5-analysis/RAxML_rootedTree.draft_tree_schmitz_rooted", format="newick")
#states = ancRECON(this_tree,this_data,method='marginal', hrm=TRUE,model="ARD")
#fitER = ace(this_data[,1], this_tree, model="ER", type="discrete")
### This fails first because it doesn't like zero length branches, and the branch from the root to first bifurcation is zero length ( I think)
# Second it fails because it doesn't seem to like something else about our data 

 
# let's try phangorn
library(phangorn)
# Transform this_data matrix into a phangorn object
this_phydat = as.phyDat(this_data, type="USER", levels = c(0,1), ambiguity = "?")
# construct a draft tree using neighbour-joiing
tr = nj(dist.ml(this_phydat))
#tr.r = root(tr, outgroup=13:14)

# Check parsimony of NJ tree and compare to maximum parsimony
parsimony(tr,this_phydat)
# 99 samples 685717
optim.parsimony(tr,this_phydat)
# 99 samples 680210

# Add gamma distribution to HKY model
fit <- update(pml(tr,this_phydat,model="HKY"),k=4)


# set up an evolutionary model for the treem0 = pml(tr, this_phydat)
m0 = pml(tr, this_phydat)
# Try Jukes-Cantor 69 model
m.jc69 = optim.pml(m0, optNni=TRUE)
# Try GTR model
m.gtr = optim.pml(m0, optNni=TRUE, model = "GTR")
# Try GTR + M + I
### THIS STAGE ALWAYS CRASHES R IN WINDOWS AND LINUX
m.gtr.G.I = optim.pml(update(m.gtr, k=4), model="GTR", optNni=TRUE, optGamma=TRUE, optInv=TRUE)
 
 
 
# Read and plot ancestral reconstruction from RAxML
# Start by reading the labelled tree- this gives identifiers for the internal nodes
labelled_tree=read.tree(file="5-analysis/RAxML_nodeLabelledRootedTree.draft_tree_rooted_anc")
plot(labelled_tree, cex=0.5)
nodelabels(labelled_tree$node.label, cex=0.5, frame="none")

# Read in the tree also containing the branch lengths
branch_length_tree=read.tree(file="5-analysis/RAxML_rootedTree.draft_tree_rooted")
plot(branch_length_tree, cex=0.5)
nodelabels(branch_length_tree$node.label, cex=0.5, frame="none")

# Read in list of loci missing from tree reconstruction due to missing data in all samples:

# Read in states of leaves of tree (THIS TAKES MANY MINUTES)
observed_states = read.table(file="4-summary/tree_samples.fas.reduced", header=FALSE, row.names=1, sep=" ", stringsAsFactors = FALSE)

# Read in ancestral states (THIS TAKES MANY MINUTES)
ancestral_states = read.table(file="5-analysis/RAxML_marginalAncestralStates.draft_tree_rooted_anc", header=FALSE, row.names=1, sep=" ", stringsAsFactors = FALSE)
no_sites = nchar(ancestral_states[1,1])

# Merge tip nodes with internal nodes
node_states = rbind.data.frame(observed_states,ancestral_states)
node_states = as.data.frame(node_states[2:nrow(node_states),])
rownames(node_states) = c(rownames(observed_states)[2:length(rownames(observed_states))], rownames(ancestral_states))

# Function to list the differences between two strings, from here: https://fabiomarroni.wordpress.com/2012/01/03/extract-different-characters-between-two-strings-of-equal-length/
list.string.diff <- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE)
{
if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
if(ignore.case)
{
a <- toupper(a)
b <- toupper(b)
}
split_seqs <- strsplit(c(a, b), split = "")
only.diff <- (split_seqs[[1]] != split_seqs[[2]])
only.diff[
(split_seqs[[1]] %in% exclude) |
(split_seqs[[2]] %in% exclude)
] <- NA
diff.info<-data.frame(which(is.na(only.diff)|only.diff),
split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
if(!show.excluded) diff.info<-na.omit(diff.info)
diff.info
}

# Unpack character strings into separate integer columns, one per locus (THIS TAKES MANY MINUTES)
#loci_states = data.frame(matrix(NA_integer_, ncol = no_sites, nrow = nrow(node_states)))
#rownames(loci_states) = rownames(node_states)
#for (locus in 1:no_sites) {
#	loci_states[,locus] = as.integer(substr(node_states[,1],locus,locus))
#} 

# Label each branch with the number of gains and losses of methylation at each end of the branch
# Make a data table to hold the results

# make a data frame to store the gains and losses for each edge
#edge_changes = data.frame(ancestor=character(), derived=character(), losses=integer(), gains=integer(), stringsAsFactors=FALSE, nrow=nrow(labelled_tree$edge))
edge_changes = data.frame(ancestor=character(nrow(labelled_tree$edge)), derived=character(nrow(labelled_tree$edge)), losses=integer(nrow(labelled_tree$edge)), gains=integer(nrow(labelled_tree$edge)), stringsAsFactors=FALSE)

annotated_tree = labelled_tree

# for each edge
for (edge in 1:nrow(labelled_tree$edge)) {
	# For ancestral node, get label from node.label in tree
	locus1 = labelled_tree$node.label[labelled_tree$edge[edge,1] - length(labelled_tree$tip.label)]
	# For derived node, label may come either from node.label or tip.label
	locus2 = ifelse(labelled_tree$edge[edge,2]<=length(labelled_tree$tip.label), labelled_tree$tip.label[labelled_tree$edge[edge,2]], labelled_tree$node.label[labelled_tree$edge[edge,2] - length(labelled_tree$tip.label)])
	z=list.string.diff(as.character(node_states[locus1,1]), as.character(node_states[locus2,1]))
#	ifelse(labelled_tree$edge[edge,2]==labelled_tree$Nnode + length(labelled_tree$tip.label), "ROOT", ifelse(labelled_tree$edge[edge,2]<=length(labelled_tree$tip.label), labelled_tree$tip.label[labelled_tree$edge[edge,2]],as.character(labelled_tree$edge[edge,2])))
#	z=list.string.diff(as.character(node_states[locus1,1]), as.character(node_states[locus2,1]))
			# Count number of 1->0 transitions
			losses = sum(as.numeric(levels(z$poly.seq.a)[z$poly.seq.a]))
			# Count number of 0->1 transitions
			gains = sum(as.numeric(levels(z$poly.seq.b)[z$poly.seq.b]))

			# Store the counts on the tree
			
			# Store the bias towards gaining methylation as a proportion of the change (figure between -1 and +1)
			bias = 2*(gains/(losses+gains) - 0.5)
			cat(paste0(locus1, " ", locus2, " Gains:", gains, " Losses:", losses, " Delta:",gains-losses,"\n"))			
			edge_changes[edge,] = list(locus1, locus2, as.integer(losses), as.integer(gains))
			#annotated_tree$node.label[edge] = 
			rm(z, losses, gains, bias)
	
}

103 CB3 Gains:13872 Losses:18044 Delta:-4172
103 CBA Gains:13650 Losses:17485 Delta:-3835
103 CBB Gains:14270 Losses:17245 Delta:-2975
103 suvh2_Stroud_2013 Gains:16690 Losses:17413 Delta:-723
103 cmt1_Stroud_2013 Gains:14319 Losses:13350 Delta:969
103 drm2_3_Zemach_2013 Gains:4632 Losses:12670 Delta:-8038
103 suvh1_Stroud_2013 Gains:14567 Losses:19206 Delta:-4639
103 dcl4_2_Stroud_2013 Gains:13377 Losses:20728 Delta:-7351
103 cmt2_6_Zemach_2013 Gains:12589 Losses:18764 Delta:-6175
103 suvh6_Stroud_2013 Gains:16803 Losses:18924 Delta:-2121
103 WC14 Gains:17163 Losses:20695 Delta:-3532
103 WCB Gains:14569 Losses:21250 Delta:-6681
103 WCA Gains:14999 Losses:21451 Delta:-6452
103 WC16 Gains:16971 Losses:20959 Delta:-3988
103 WT_Zemach_2013 Gains:140 Losses:3422 Delta:-3282
103 cmt2_5_Zemach_2013 Gains:9688 Losses:20050 Delta:-10362
103 WT_Stroud_2014 Gains:15178 Losses:15935 Delta:-757
103 suvh7_Stroud_2013 Gains:13925 Losses:16530 Delta:-2605
103 suvh456_Stroud_2013 Gains:16438 Losses:16155 Delta:283
103 suvh5_Stroud_2013 Gains:10727 Losses:17839 Delta:-7112
103 ibm1_Stroud_2013 Gains:12270 Losses:13861 Delta:-1591
103 WT_ibm1_r2_CD_2012 Gains:10512 Losses:21350 Delta:-10838
103 WT_ibm1_r1_CD_2012 Gains:12755 Losses:19263 Delta:-6508
103 ibm1_6_r2_CD_2012 Gains:10852 Losses:14890 Delta:-4038
103 ibm1_6_r1_CD_2012 Gains:2289 Losses:8316 Delta:-6027
103 suvh3_Stroud_2013 Gains:10077 Losses:13092 Delta:-3015
103 JAY_1 Gains:0 Losses:0 Delta:0
103 WT_Lin_2017 Gains:5905 Losses:2830 Delta:3075
103 WTr2_Stroud_2013 Gains:5677 Losses:7563 Delta:-1886
103 WTr3_Stroud_2013 Gains:3283 Losses:2596 Delta:687
103 WTr1_Harris_2016 Gains:7904 Losses:7222 Delta:682
103 WTr2_Harris_2016 Gains:6767 Losses:8692 Delta:-1925
103 WT_Groth_2016 Gains:12564 Losses:10094 Delta:2470
103 SRR342381 Gains:16840 Losses:19846 Delta:-3006
103 SRR342391 Gains:11043 Losses:13921 Delta:-2878
103 ERR046552 Gains:11565 Losses:19279 Delta:-7714
103 ERR046553 Gains:11810 Losses:18881 Delta:-7071
103 SRR342378 Gains:15975 Losses:21759 Delta:-5784
103 ERR046560 Gains:14071 Losses:22425 Delta:-8354
103 ERR046559 Gains:13379 Losses:22498 Delta:-9119
103 ERR046575 Gains:7864 Losses:13626 Delta:-5762
103 ERR046574 Gains:8675 Losses:14955 Delta:-6280
103 SRR342385 Gains:15514 Losses:19344 Delta:-3830
103 SRR342354 Gains:10865 Losses:13049 Delta:-2184
103 SRR342353 Gains:11321 Losses:13026 Delta:-1705
103 SRR342383 Gains:8318 Losses:14553 Delta:-6235
103 SRR342348 Gains:10585 Losses:15608 Delta:-5023
103 ERR046549 Gains:5544 Losses:11181 Delta:-5637
103 SRR342382 Gains:7252 Losses:14977 Delta:-7725
103 SRR342347 Gains:9628 Losses:17461 Delta:-7833
103 ERR046548 Gains:7033 Losses:15004 Delta:-7971
103 ERR046547 Gains:8802 Losses:17745 Delta:-8943
103 ERR046546 Gains:6485 Losses:11181 Delta:-4696
103 ERR046564 Gains:12491 Losses:35456 Delta:-22965
103 ERR046563 Gains:14289 Losses:34457 Delta:-20168
103 SRR342390 Gains:13735 Losses:30789 Delta:-17054
103 SRR342380 Gains:15452 Losses:32491 Delta:-17039
103 ERR046562 Gains:14286 Losses:21007 Delta:-6721
103 ERR046561 Gains:14452 Losses:22932 Delta:-8480
103 SRR342389 Gains:14354 Losses:20065 Delta:-5711
103 SRR342379 Gains:15574 Losses:21823 Delta:-6249
103 SRR342384 Gains:15192 Losses:20279 Delta:-5087
103 SRR342349 Gains:15925 Losses:20885 Delta:-4960
103 ERR046555 Gains:14921 Losses:22582 Delta:-7661
103 ERR046554 Gains:16088 Losses:21853 Delta:-5765
103 ERR046556 Gains:10738 Losses:20134 Delta:-9396
103 cmt2_3_Zemach_2013 Gains:436 Losses:4574 Delta:-4138
103 JAY_3 Gains:2499 Losses:18251 Delta:-15752
103 cmt2_4_Zemach_2013 Gains:7544 Losses:12480 Delta:-4936
103 dcl2_1_Stroud_2013 Gains:16284 Losses:27606 Delta:-11322
103 cmt2_Stroud_2013 Gains:16206 Losses:19454 Delta:-3248
103 JAY_4 Gains:12397 Losses:17486 Delta:-5089
103 drm3_1_Stroud_2013 Gains:17378 Losses:21764 Delta:-4386
103 dcl3_1_Stroud_2013 Gains:18500 Losses:20093 Delta:-1593
103 ibm1_1GL17_Ito_2015 Gains:15573 Losses:13395 Delta:2178
103 ibm1_3GL6_2_Ito_2015 Gains:18227 Losses:18201 Delta:26
103 ibm1_3GL6_1_Ito_2015 Gains:18649 Losses:17514 Delta:1135
103 ibm1_3GL4_3_Ito_2015 Gains:16581 Losses:15921 Delta:660
103 ibm1_3GL4_1_Ito_2015 Gains:18024 Losses:17493 Delta:531
103 ibm1_1GL23_Ito_2015 Gains:16884 Losses:15467 Delta:1417
103 ibm1_1GL21_Ito_2015 Gains:16059 Losses:14895 Delta:1164
103 suvh8_Stroud_2013 Gains:15741 Losses:18764 Delta:-3023
103 suvh4_Stroud_2013 Gains:16110 Losses:14118 Delta:1992
103 DDB Gains:17402 Losses:20570 Delta:-3168
103 DDA Gains:17582 Losses:20593 Delta:-3011
103 DD14B Gains:17875 Losses:21594 Delta:-3719
103 DD11 Gains:17827 Losses:21724 Delta:-3897
103 drm12r2_Stroud_2013 Gains:15661 Losses:15413 Delta:248
103 drm12r1_Stroud_2013 Gains:9730 Losses:14078 Delta:-4348
103 cmt2cmt3_Stroud_2014 Gains:17038 Losses:15379 Delta:1659
103 drm12cmt23_Stroud_2014 Gains:19137 Losses:14470 Delta:4667
103 ddcc_Lin_2017 Gains:20970 Losses:13652 Delta:7318
103 drm12cmt2_Stroud_2014 Gains:18826 Losses:14482 Delta:4344
103 drm12cmt3_Stroud_2013 Gains:20000 Losses:17524 Delta:2476
103 cmt3_Stroud_2013 Gains:16315 Losses:15933 Delta:382
103 TD13 Gains:15137 Losses:19586 Delta:-4449
103 TD18 Gains:14795 Losses:20046 Delta:-5251
103 TDB Gains:14615 Losses:18853 Delta:-4238
103 TDA Gains:14035 Losses:20023 Delta:-5988
103 cmt3_12_Zemach_2013 Gains:9471 Losses:15880 Delta:-6409
103 cmt2_3_Stroud_2014 Gains:14990 Losses:14714 Delta:276
103 CB18 Gains:15280 Losses:18060 Delta:-2780
103 103 Gains:0 Losses:0 Delta:0
103 104 Gains:16342 Losses:18115 Delta:-1773
103 105 Gains:16230 Losses:17801 Delta:-1571
103 106 Gains:16282 Losses:18027 Delta:-1745
103 107 Gains:15203 Losses:15235 Delta:-32
103 108 Gains:13998 Losses:13422 Delta:576
103 109 Gains:13637 Losses:12961 Delta:676
103 110 Gains:13520 Losses:13315 Delta:205
103 111 Gains:12686 Losses:13546 Delta:-860
103 112 Gains:13595 Losses:12625 Delta:970
103 113 Gains:13088 Losses:12125 Delta:963
103 114 Gains:11784 Losses:11033 Delta:751
103 115 Gains:10396 Losses:10451 Delta:-55
103 116 Gains:8730 Losses:10205 Delta:-1475
103 117 Gains:7165 Losses:7771 Delta:-606
103 118 Gains:5572 Losses:5651 Delta:-79
103 119 Gains:7469 Losses:7902 Delta:-433
103 120 Gains:11096 Losses:13660 Delta:-2564
103 121 Gains:7997 Losses:8903 Delta:-906
103 122 Gains:16719 Losses:20611 Delta:-3892
103 123 Gains:17231 Losses:21709 Delta:-4478
103 124 Gains:17294 Losses:22122 Delta:-4828
103 125 Gains:8321 Losses:10041 Delta:-1720
103 126 Gains:7348 Losses:9934 Delta:-2586
103 127 Gains:4906 Losses:4769 Delta:137
103 128 Gains:3665 Losses:3724 Delta:-59
103 129 Gains:7479 Losses:8913 Delta:-1434
103 130 Gains:2897 Losses:3053 Delta:-156
103 131 Gains:1718 Losses:1946 Delta:-228
103 132 Gains:4949 Losses:6710 Delta:-1761
103 133 Gains:5151 Losses:8146 Delta:-2995
103 134 Gains:10123 Losses:15098 Delta:-4975
103 135 Gains:11092 Losses:17037 Delta:-5945
103 136 Gains:883 Losses:941 Delta:-58
103 137 Gains:851 Losses:633 Delta:218
103 138 Gains:3239 Losses:2802 Delta:437
103 139 Gains:2030 Losses:1242 Delta:788
103 140 Gains:3852 Losses:3324 Delta:528
103 141 Gains:5749 Losses:6072 Delta:-323
103 142 Gains:6613 Losses:7514 Delta:-901
103 143 Gains:7199 Losses:8262 Delta:-1063
103 144 Gains:7882 Losses:8896 Delta:-1014
103 145 Gains:14220 Losses:18655 Delta:-4435
103 146 Gains:14295 Losses:19417 Delta:-5122
103 147 Gains:14785 Losses:20370 Delta:-5585
103 148 Gains:8328 Losses:9585 Delta:-1257
103 149 Gains:14725 Losses:19780 Delta:-5055
103 150 Gains:15398 Losses:21435 Delta:-6037
103 151 Gains:15623 Losses:22259 Delta:-6636
103 152 Gains:15585 Losses:22536 Delta:-6951
103 153 Gains:15638 Losses:22682 Delta:-7044
103 154 Gains:8537 Losses:10587 Delta:-2050
103 155 Gains:10705 Losses:12790 Delta:-2085
103 156 Gains:8355 Losses:11451 Delta:-3096
103 157 Gains:9672 Losses:15059 Delta:-5387
103 158 Gains:8029 Losses:12373 Delta:-4344
103 159 Gains:8021 Losses:13228 Delta:-5207
103 160 Gains:8920 Losses:16344 Delta:-7424
103 161 Gains:8302 Losses:13647 Delta:-5345
103 162 Gains:8814 Losses:14737 Delta:-5923
103 163 Gains:8384 Losses:11531 Delta:-3147
103 164 Gains:13056 Losses:30844 Delta:-17788
103 165 Gains:13943 Losses:34017 Delta:-20074
103 166 Gains:14511 Losses:32259 Delta:-17748
103 167 Gains:14111 Losses:20677 Delta:-6566
103 168 Gains:14842 Losses:21994 Delta:-7152
103 169 Gains:15498 Losses:21897 Delta:-6399
103 170 Gains:14451 Losses:19601 Delta:-5150
103 171 Gains:15421 Losses:20822 Delta:-5401
103 172 Gains:15925 Losses:21701 Delta:-5776
103 173 Gains:16139 Losses:22368 Delta:-6229
103 174 Gains:8658 Losses:11625 Delta:-2967
103 175 Gains:8055 Losses:12569 Delta:-4514
103 176 Gains:7704 Losses:12944 Delta:-5240
103 177 Gains:11296 Losses:14263 Delta:-2967
103 178 Gains:12223 Losses:14700 Delta:-2477
103 179 Gains:12721 Losses:15307 Delta:-2586
103 180 Gains:13734 Losses:11163 Delta:2571
103 181 Gains:14705 Losses:12413 Delta:2292
103 182 Gains:16208 Losses:15202 Delta:1006
103 183 Gains:15764 Losses:13677 Delta:2087
103 184 Gains:16331 Losses:15920 Delta:411
103 185 Gains:16493 Losses:13918 Delta:2575
103 186 Gains:13660 Losses:12702 Delta:958
103 187 Gains:13963 Losses:12859 Delta:1104
103 188 Gains:14759 Losses:13369 Delta:1390
103 189 Gains:15638 Losses:15226 Delta:412
103 190 Gains:18700 Losses:20901 Delta:-2201
103 191 Gains:19156 Losses:21339 Delta:-2183
103 192 Gains:19250 Losses:21510 Delta:-2260
103 193 Gains:16277 Losses:16484 Delta:-207
103 194 Gains:16712 Losses:13797 Delta:2915
103 195 Gains:17758 Losses:13529 Delta:4229
103 196 Gains:19841 Losses:14053 Delta:5788
103 197 Gains:20660 Losses:14766 Delta:5894
103 198 Gains:14564 Losses:14128 Delta:436
103 199 Gains:14678 Losses:16630 Delta:-1952
103 200 Gains:16555 Losses:19988 Delta:-3433
103 201 Gains:16756 Losses:20538 Delta:-3782
103 202 Gains:16784 Losses:20677 Delta:-3893
103 203 Gains:0 Losses:0 Delta:0

# Add a set of new tips to the existing tips, to store the tip names and leave node names free to accumulate methylation values
#library(phytools)
#for (tip_name in annotated_tree$tip.label) {
#	tip_number = match(tip_name, annotated_tree$tip.label)
#	annotated_tree = bind.tip(tree=annotated_tree, tip.label=paste0("X_",annotated_tree$tip.label[tip_number]), edge.length=0, where=tip_number, position=0)
#}


# Traverse the tree from the root, accumulating net delta at each node

find_children <- function(this_node, this_tree) {
	node_ids = this_tree$edge[,2][this_tree$node.label[this_tree$edge[,1] - length(this_tree$tip.label)] == this_node]
	ifelse(node_ids > length(this_tree$tip.label),this_tree$node.label[node_ids - ifelse(node_ids>length(this_tree$tip.label), length(this_tree$tip.label), 0)], this_tree$tip.label[node_ids])
}

annotate_children <- function(this_node, this_tree, this_value) {
	#cat(paste0("Entry: ", this_node, " ", this_value, "\n"))
	for (this_child in find_children(this_node, this_tree)) {
		this_delta = edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_child),]$gains - edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_child),]$losses
	cat(paste0(this_node, " ", this_child, " ", this_value, " ", this_delta, "\n"))

		this_tree = annotate_children(this_child, this_tree, this_value + this_delta)
	}
	this_tree$node.label[this_tree$node.label == this_node] = this_value
	#this_tree$tip.label[this_tree$tip.label == this_node] = this_value
	this_tree
}


annotated_tree = annotate_children(this_node="ROOT", this_tree=annotated_tree, this_value=0)
annotated_tree$edge.length = branch_length_tree$edge.length
write.tree(annotated_tree, file="annotated_tree.nwk")

# This produced slightly weird results with lots of branches showing lots of moves only in one direction.  Check this!



# Correlate Schmitz results with those from previous analysis
all_samples_meth_status <- data.frame(readRDS(paste0("../Jay-SMP_variation/5-analysis/SRA035939/",project_id,"_",meth_context,"_all_samples_meth_status.rds")))




### Let's try a new way to estimate ancestral reconstruction, then repeat the above analyses

# Function to traverse tree and infer ancestral states by parsimony
# Start from tips
# Annotate tips with probabilities of each being U or M
# If reads indicate U, p=1,0
# If reads indicate M, p=0,1
# If reads indicate ?, p=0.5,0.5
# At parental node, p = mean of child nodes

library(phangorn)

# Read RAxML rooted tree
#labelled_tree=read.tree(file="5-analysis/RAxML_RootedTree.draft_tree_rooted")
#labelled_tree=read.tree(file="5-analysis/RAxML_RootedTree.draft_tree_rooted")
labelled_tree=read.tree(file="5-analysis/RAxML_bestTree.99_plus_me.nwk")
plot(labelled_tree, cex=0.5)
nodelabels(labelled_tree$node.label, cex=0.5, frame="none")

# Reroot the tree using the JAY_1 sample (Col-0)
# This page provides discussion about different approaches to rooting trees in R:
# https://www.biostars.org/p/332030/

rerooted_tree = phytools::reroot(labelled_tree, which(labelled_tree$tip.label=="JAY_1"))
# explicitly: 
#rerooted_tree = phytools::reroot(labelled_tree, 27)
# node JAY_1 shares with rest of tree: 
#rerooted_tree = phytools::reroot(labelled_tree, 136)

plot(rerooted_tree, cex=0.5)
nodelabels(rerooted_tree$node.label, cex=0.5, frame="none")

write.tree(rerooted_tree, file="rerooted_tree.nwk")

# Read in states of leaves of tree (THIS TAKES MANY MINUTES)
#observed_states = read.table(file="4-summary/tree_samples.fas.reduced", header=FALSE, row.names=1, sep=" ", stringsAsFactors = FALSE)


### Testing version
#z=as.matrix(read.table(file="5-analysis/tree_samples_data_matrix_test.txt", sep="\t"))
#rownames(z)=labelled_tree$tip.label
#tree_samples_data_matrix_backup = tree_samples_data_matrix
#tree_samples_data_matrix = z
###

# Sort tree_samples_data_matrix into the tip.label order from labelled_tree
tree_samples_data_matrix = tree_samples_data_matrix[match(labelled_tree$tip.label, rownames(tree_samples_data_matrix)),]
##### SHOULDN'T THIS BE SORT INTO ORDER OF REROOTED TREE?


# Initialise a matrix to contain the probabilities of U and M at each internal node. Initialise to -1
labelled_tree$node.Mprobs = matrix(rep(-1,len=ncol(tree_samples_data_matrix)*(labelled_tree$Nnode)), ncol=ncol(tree_samples_data_matrix))
# Don't bother estimating Uprobs as they are 1-Mprobs
#labelled_tree$node.Uprobs = matrix(rep(-1,len=ncol(tree_samples_data_matrix)*(labelled_tree$Nnode)), ncol=ncol(tree_samples_data_matrix))

# Traverse the tree from the root to tips recursively; on the way back accumulate marginal ancestral state probabilities at each node

# Accessory function to find children of a node (this version works with node numbers only, not labels)
find_children_ids <- function(this_node, this_tree) {
	#node_ids = this_tree$edge[,2][this_tree$node.label[this_tree$edge[,1] - length(this_tree$tip.label)] == this_node]
	node_ids = this_tree$edge[,2][this_tree$edge[,1] == this_node]
	#ifelse(node_ids > length(this_tree$tip.label),this_tree$node.label[node_ids - ifelse(node_ids>length(this_tree$tip.label), length(this_tree$tip.label), 0)], this_tree$tip.label[node_ids])
	node_ids
}


# Use phangorn for Fitch's parsimony reconstruction:
# convert to phyDat object
#tree_samples_data_matrix.phyDat = phyDat(tree_samples_data_matrix, type="USER", levels=c("0","1"), ambiguity = c(NA))
#tree_samples_data_matrix.phyDat = phyDat(z1, type="USER", levels=c(0,1), ambiguity = c(NA))
# convert data frame to matrix - phyDat assumes one row per taxon for matrix and one col per taxon for df
tree_samples_data_matrix.phyDat = phyDat(as.matrix(z1), type="USER", levels=c(0,1), ambiguity = c(NA))
# Reconstruct ancestral states using maximum parsimony
anc.mpr = ancestral.pars(rerooted_tree, tree_samples_data_matrix.phyDat, "MPR")

#ACCTRAN (accelerated transformation) [Farris, 1970] reduces the number of node state ambiguities by forcing the state changes to be performed as close to the root as possible, and therefore prioritises the reverse mutations.
#DELTRAN (delayed transformation) [Swofford and Maddison, 1987] reduces the number of node state ambiguities by making the changes as close to the tips as possible, hence prioritizing parallel mutations.
# The acctran method assigns edge length and internal nodes to the tree.  In the example (https://cran.r-project.org/web/packages/phangorn/vignettes/Ancestral.pdf) acctran resolved unknowns more convincingly
#anc.acctran = ancestral.pars(rerooted_tree, tree_samples_data_matrix.phyDat, "ACCTRAN")



# Make a data frame containing the node states from the ancestral reconstruction
node_states = data.frame(stringsAsFactors=FALSE)
for (node_id in names(anc.mpr)) {
	node_states=rbind.data.frame(node_states, paste(ifelse(anc.mpr[[node_id]][,2] == 0.5,"-",ifelse(anc.mpr[[node_id]][,2] == 0,"0",ifelse(anc.mpr[[node_id]][,2] == 1,"1","-"))),collapse=""), stringsAsFactors=FALSE)
}
rownames(node_states)=names(anc.mpr)

# Plot methylation status on tree for some sites to check reasonableness:
example_site = 103

# indirection for example_site:
# length(unique(attr(anc.mpr, "index"))) = 572380 (number of unique phylogenetic patterns)
# attr(anc,mpr, index") has 1085658 entries (all CG sites in GBM genes)
# attr(anc.mpr, "index") returns a vector mapping the unique patterns to the original 'all GBM gene sites' entries in z1/tree_samples_data_matrix
# colnames(tree_samples_data_matrix) returns a vector mapping the z1 entries to the original 'all CG sites' entries in tree_samples_meth_status/meth_status_sites.gr

# to plot using a 'unique phylogenetic pattern' index from anc.mpr
#plotAnc(tree=rerooted_tree, data=anc.mpr, i=example_site, cex.pie=0.2,  col=c("blue","red"))

# to plot using the ID of a GBM gene site from tree_samples_data_matrix/z1:
#plotAnc(tree=rerooted_tree, data=anc.mpr, i=attr(anc.mpr, "index")[example_site], cex.pie=0.2,  col=c("blue","red"))

# to plot using the ID of a site from 'all sites' tree_samples_meth_status/meth_status_sites.gr:
plotAnc(tree=rerooted_tree, data=anc.mpr, i=attr(anc.mpr, "index")[match(example_site,colnames(tree_samples_data_matrix))], cex.pie=0.2,  col=c("blue","red"))



# Function to list the differences between two strings, from here: https://fabiomarroni.wordpress.com/2012/01/03/extract-different-characters-between-two-strings-of-equal-length/
list.string.diff <- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE)
{
if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
if(ignore.case)
{
a <- toupper(a)
b <- toupper(b)
}
split_seqs <- strsplit(c(a, b), split = "")
only.diff <- (split_seqs[[1]] != split_seqs[[2]])
only.diff[
(split_seqs[[1]] %in% exclude) |
(split_seqs[[2]] %in% exclude)
] <- NA
diff.info<-data.frame(which(is.na(only.diff)|only.diff),
split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
if(!show.excluded) diff.info<-na.omit(diff.info)
diff.info
}

# make a data frame to store the sum of gains and losses for each edge
edge_changes = data.frame(ancestor=character(nrow(rerooted_tree$edge)), derived=character(nrow(rerooted_tree$edge)), losses=integer(nrow(rerooted_tree$edge)), gains=integer(nrow(rerooted_tree$edge)), stringsAsFactors=FALSE)

# make a data frame to store the phylogenetic branch and site of each gain and loss 
site_branch_changes = data.frame(ancestor = character(), derived=character(), position=integer(), losses=integer(), gains=integer())
#site_branch_changes = data.frame(stringsAsFactors=FALSE)

# for each edge
for (edge in 1:nrow(rerooted_tree$edge)) {
	# For ancestral node, get label from node.label in tree
	locus1 = rerooted_tree$edge[edge,1]
	# For derived node, label may come either from node.label or tip.label
	#locus2 = rerooted_tree$edge[edge,2]
	locus2 = ifelse(rerooted_tree$edge[edge,2]<=length(rerooted_tree$tip.label), rerooted_tree$tip.label[rerooted_tree$edge[edge,2]], rerooted_tree$edge[edge,2])
	z=list.string.diff(as.character(node_states[locus1,1]), as.character(node_states[locus2,1]))
#	ifelse(rerooted_tree$edge[edge,2]==rerooted_tree$Nnode + length(rerooted_tree$tip.label), "ROOT", ifelse(rerooted_tree$edge[edge,2]<=length(rerooted_tree$tip.label), rerooted_tree$tip.label[rerooted_tree$edge[edge,2]],as.character(rerooted_tree$edge[edge,2])))
#	z=list.string.diff(as.character(node_states[locus1,1]), as.character(node_states[locus2,1]))
      # store the individual changes
      site_branch_changes = rbind.data.frame(site_branch_changes, cbind("ancestor" = rep(locus1, nrow(z)), "derived" = rep(locus2, nrow(z)), "position" = z$position, "losses" = as.numeric(levels(z$poly.seq.a)[z$poly.seq.a]), "gains" = as.numeric(levels(z$poly.seq.b)[z$poly.seq.b])))

      # Count number of 1->0 transitions
			losses = sum(as.numeric(levels(z$poly.seq.a)[z$poly.seq.a]))
			# Count number of 0->1 transitions
			gains = sum(as.numeric(levels(z$poly.seq.b)[z$poly.seq.b]))

			# Store the counts on the tree
			
			# Store the bias towards gaining methylation as a proportion of the change (figure between -1 and +1)
			bias = 2*(gains/(losses+gains) - 0.5)
			cat(paste0(locus1, " ", locus2, " Gains:", gains, " Losses:", losses, " Delta:",gains-losses,"\n"))			
			edge_changes[edge,] = list(locus1, locus2, as.integer(losses), as.integer(gains))
			#annotated_tree$node.label[edge] = 
			#rm(z, losses, gains, bias)
	
}

site_branch_changes$position = as.numeric(site_branch_changes$position)
site_branch_changes$losses = as.numeric(site_branch_changes$losses)
site_branch_changes$gains = as.numeric(site_branch_changes$gains)

# summarise gains and losses by site
site_changes = aggregate(cbind(site_branch_changes$losses,site_branch_changes$gains), by=list(site_branch_changes$position), FUN=sum)
colnames(site_changes) = c("position", "losses", "gains")

ggplot(site_changes) + geom_bin2d(aes(x=losses, y=gains), binwidth=1) + scale_fill_gradient(name = "count", trans = "log")

# summarise gains and losses by site, but split by terminal and non-terminal losses
site_internal_branch_changes = site_branch_changes[!is.na(as.integer(site_branch_changes$derived)),]
site_external_branch_changes = site_branch_changes[is.na(as.integer(site_branch_changes$derived)),]

site_internal_changes = aggregate(cbind(site_internal_branch_changes$losses,site_internal_branch_changes$gains), by=list(site_internal_branch_changes$position), FUN=sum)
colnames(site_internal_changes) = c("position", "losses", "gains")
site_external_changes = aggregate(cbind(site_external_branch_changes$losses,site_external_branch_changes$gains), by=list(site_external_branch_changes$position), FUN=sum)
colnames(site_external_changes) = c("position", "losses", "gains")

ggplot(site_internal_changes) + geom_bin2d(aes(x=losses, y=gains), binwidth=1) + scale_fill_gradient(name = "count", trans = "log") +xlim(0,15) + ylim(0,15)
ggplot(site_external_changes) + geom_bin2d(aes(x=losses, y=gains), binwidth=1) + scale_fill_gradient(name = "count", trans = "log") +xlim(0,15) + ylim(0,15)


# plot histograms of numbers of internal changes per site
hist(site_internal_changes$gains, breaks=seq(-0.5,15, 1))
hist(site_internal_changes$losses, breaks=seq(-0.5,15, 1))

  
# split changes between those that start methylated and those unmethylated:
str_length(node_states["JAY_1",])
length((anc.mpr$JAY_1==1))


# Now we have number of internal changes per site, focus on these.

gene_site_count = rep(0, length(GBM_gene_loci.gr))
gene_initial_M_count = rep(0, length(GBM_gene_loci.gr))
gene_internal_losses = rep(0, length(GBM_gene_loci.gr))
gene_internal_gains = rep(0, length(GBM_gene_loci.gr))
gene_external_losses = rep(0, length(GBM_gene_loci.gr))
gene_external_gains = rep(0, length(GBM_gene_loci.gr))
gene_internal_lost_once_sites = rep(0, length(GBM_gene_loci.gr))
gene_internal_lost_more_sites = rep(0, length(GBM_gene_loci.gr))
gene_internal_gain_once_sites = rep(0, length(GBM_gene_loci.gr))
gene_internal_gain_more_sites = rep(0, length(GBM_gene_loci.gr))

for (this_gene_no in 1:length(GBM_gene_loci.gr)) {
  # identify the meth_status_sites.gr sites relevant to the gene
  this_gene_sites = queryHits(findOverlaps(meth_status_sites.gr, GBM_gene_loci.gr[this_gene_no]))

  # position in e.g. site_internal_changes is 1 of 572380 - indirect through:  attr(anc.mpr, "index")[match(example_site,colnames(tree_samples_data_matrix))]
  this_gene_positions = attr(anc.mpr, "index")[match(this_gene_sites,colnames(tree_samples_data_matrix))]
    
  # count No. methylated sites in JAY_1 accession
  #this_gene_initial_pattern = attr(tree_samples_data_matrix.phyDat, "levels")[tree_samples_data_matrix.phyDat$JAY_1[this_gene_positions]]
  # count No. methylated sites in ancestral estimate of root of MA lines
  this_gene_initial_pattern = 1-anc.mpr$`142`[,1][this_gene_positions]
  gene_initial_M_count[this_gene_no] = sum(this_gene_initial_pattern[!is.na(this_gene_initial_pattern)])

  gene_site_count[this_gene_no] = length(this_gene_initial_pattern)
    
  # count how many of its sites become methylated on an internal branch, and how many times across the tree
  gene_internal_losses[this_gene_no] = sum(site_internal_changes[site_internal_changes$position %in% this_gene_positions,]$losses)
  gene_internal_gains[this_gene_no] = sum(site_internal_changes[site_internal_changes$position %in% this_gene_positions,]$gains)
  gene_external_losses[this_gene_no] = sum(site_external_changes[site_external_changes$position %in% this_gene_positions,]$losses)
  gene_external_gains[this_gene_no] = sum(site_external_changes[site_external_changes$position %in% this_gene_positions,]$gains)

  gene_internal_lost_once_sites[this_gene_no] = sum(site_internal_changes[site_internal_changes$position %in% this_gene_positions,]$losses==1)
  gene_internal_lost_more_sites[this_gene_no] = sum(site_internal_changes[site_internal_changes$position %in% this_gene_positions,]$losses>1)
  gene_internal_gain_once_sites[this_gene_no] = sum(site_internal_changes[site_internal_changes$position %in% this_gene_positions,]$gains==1)
  gene_internal_gain_more_sites[this_gene_no] = sum(site_internal_changes[site_internal_changes$position %in% this_gene_positions,]$gains>1)
  
    # count how many of its sites become unmethylated on an internal branch, and how many times across the tree
}

# Internal and external losses, and gains, and both, are correlated among genes (related to gene length?)
cor(gene_internal_losses, gene_external_losses)
#[1] 0.8785015
cor(gene_internal_gains, gene_external_gains)
#[1] 0.8123835
cor(gene_internal_losses, gene_internal_gains)
#[1] 0.7124217
cor(gene_external_losses, gene_external_gains)
#[1] 0.7046074
cor(gene_internal_lost_once_sites, gene_internal_lost_more_sites)
#[1] 0.6844108
cor(gene_internal_lost_once_sites, gene_internal_gain_once_sites)
#[1] 0.6455213
cor(gene_internal_gain_once_sites, gene_internal_gain_more_sites)
#[1] 0.544511
cor(gene_internal_lost_more_sites, gene_internal_gain_more_sites)
#[1] 0.5682163
head(gene_internal_lost_once_sites)

sum(gene_internal_gain_once_sites>0)
#[1] 12186
sum(gene_internal_gain_more_sites>0)
#[1] 9690
sum((gene_internal_gain_more_sites>0) & (gene_internal_gain_once_sites>0))
#[1] 9266
sum((gene_internal_gain_more_sites==0) & (gene_internal_gain_once_sites>0))
#[1] 2920
sum(gene_internal_lost_once_sites>0 & gene_internal_lost_more_sites==0)
#[1] 1550

sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==0)
#[1] 289
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==1)
#[1] 281
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==2)
#[1] 292
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==3)
#[1] 218
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==4)
#[1] 197
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==5)
#[1] 164
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==6)
#[1] 147
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==7)
#[1] 121
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==8)
#[1] 110
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==9)
#[1] 114
sum(gene_internal_gain_once_sites>0 & gene_internal_gain_more_sites==0 & gene_initial_M_count==10)
#[1] 97




# Plot the phylogenetic picture of changes in the genes with at least one site showing a gain once, but no sites which gain more than once, and where initial M site count is 1
GBM_gene_loci[queryHits(findOverlaps(makeGRangesFromDataFrame(df = GBM_gene_loci, start.field = "start_site", end.field = "end_site",seqnames.field = "Chromosome"), GBM_gene_loci.gr[(gene_internal_gain_once_sites>0) & (gene_internal_gain_more_sites==0) & (gene_initial_M_count==1)])),"gene_ID"]



# plot rate of internal loss vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_internal_losses)), aes(x=gene_initial_M_count, y=gene_internal_losses/gene_initial_M_count)) + geom_boxplot(aes(group=gene_initial_M_count)))
# plot count of internal losses vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_internal_losses)), aes(x=gene_initial_M_count, y=gene_internal_losses)) + geom_boxplot(aes(group=gene_initial_M_count)) + geom_smooth(method='lm'))

# plot rate of internal gain vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_internal_gains)), aes(x=gene_initial_M_count, y=gene_internal_gains/gene_initial_M_count)) + geom_boxplot(aes(group=gene_initial_M_count)))
# plot count of internal gains vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_internal_gains)), aes(x=gene_initial_M_count, y=gene_internal_gains)) + geom_boxplot(aes(group=gene_initial_M_count)) + geom_smooth(method='lm'))

# plot rate of external loss vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_external_losses)), aes(x=gene_initial_M_count, y=gene_external_losses/gene_initial_M_count)) + geom_boxplot(aes(group=gene_initial_M_count)))
# plot count of external loss vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_external_losses)), aes(x=gene_initial_M_count, y=gene_external_losses)) + geom_boxplot(aes(group=gene_initial_M_count)) + geom_smooth(method='lm'))

# plot rate of external gain vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_external_gains)), aes(x=gene_initial_M_count, y=gene_external_gains/gene_initial_M_count)) + geom_boxplot(aes(group=gene_initial_M_count)))
# plot count of external gain vs. number of initial sites
print(ggplot(cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_external_gains)), aes(x=gene_initial_M_count, y=gene_external_gains)) + geom_boxplot(aes(group=gene_initial_M_count)) + geom_smooth(method='lm'))

# try a tern plot
ggtern() + geom_point(data=cbind(as.data.frame(gene_initial_M_count), as.data.frame(gene_internal_losses), as.data.frame(gene_internal_gains)), aes(gene_initial_M_count, gene_internal_losses, gene_internal_gains))





# Recursive function to traverse the tree from root to tip, accumulating net delta at each node

#annotate_children <- function(this_node, this_tree, this_value) {
#	cat(paste0("Entry: ", this_node, " ", this_value, "\n"))
#	# Find the children of this node, accumulate their net delta, and annotate them recursively
#	for (this_child in find_children_ids(this_node, this_tree)) {
#		# Calculate the delta from this_node to child_node by lookup from the edge_changes table
#		if (this_child <= length(this_tree$tip.label)) {
#			# child is a tip node
#			this_delta = edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_tree$tip.label[this_child]),]$gains - edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_tree$tip.label[this_child]),]$losses		
#		} else {
#			# child is an internal node
#			this_delta = edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_child),]$gains - edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_child),]$losses
#		}
#		cat(paste0(this_node, " ", this_child, " ", this_value, " ", this_delta, "\n"))
#		
#		# Pass the accumulated delta for the child node into the recurse
#		this_tree = annotate_children(this_child, this_tree, this_value + this_delta)
#	}
#	
#	# store the accumulated delta for this node in the node or tip label of the tree itself
#	if(this_node>length(this_tree$tip.label)) {
#		this_tree$node.label[this_node-length(this_tree$tip.label)] = this_value
#	} else {
#		this_tree$tip.label[this_node] = paste0(this_tree$tip.label[this_node], " ", as.character(this_value))
#	}
#	
#	# pass the updated tree back down the line
#	this_tree
#}


# Alternative version that calculates difference between each node and root node, rather than accumulating differences

# make a data frame to store the gains and losses for each node
node_delta = data.frame(ancestor=character(length(rerooted_tree$tip.label) + rerooted_tree$Nnode), derived=character(length(rerooted_tree$tip.label) + rerooted_tree$Nnode), losses=integer(length(rerooted_tree$tip.label) + rerooted_tree$Nnode), gains=integer(length(rerooted_tree$tip.label) + rerooted_tree$Nnode), stringsAsFactors=FALSE)

root_node = 103

# for each node
for (node in 1:(length(rerooted_tree$tip.label)+rerooted_tree$Nnode)) {
	# For ancestral node, get label from node.label in tree
	locus1 = root_node
	# For derived node, label may come either from node.label or tip.label
	#locus2 = rerooted_tree$edge[edge,2]
	locus2 = ifelse(node<=length(rerooted_tree$tip.label), rerooted_tree$tip.label[node], node)
	z=list.string.diff(as.character(node_states[locus1,1]), as.character(node_states[locus2,1]))
#	ifelse(rerooted_tree$edge[edge,2]==rerooted_tree$Nnode + length(rerooted_tree$tip.label), "ROOT", ifelse(rerooted_tree$edge[edge,2]<=length(rerooted_tree$tip.label), rerooted_tree$tip.label[rerooted_tree$edge[edge,2]],as.character(rerooted_tree$edge[edge,2])))
#	z=list.string.diff(as.character(node_states[locus1,1]), as.character(node_states[locus2,1]))
			# Count number of 1->0 transitions
			losses = sum(as.numeric(levels(z$poly.seq.a)[z$poly.seq.a]))
			# Count number of 0->1 transitions
			gains = sum(as.numeric(levels(z$poly.seq.b)[z$poly.seq.b]))

			# Store the counts on the tree
			
			# Store the bias towards gaining methylation as a proportion of the change (figure between -1 and +1)
			bias = 2*(gains/(losses+gains) - 0.5)
			cat(paste0(locus1, " ", locus2, " Gains:", gains, " Losses:", losses, " Delta:",gains-losses,"\n"))			
			node_delta[node,] = list(locus1, locus2, as.integer(losses), as.integer(gains))
			#annotated_tree$node.label[edge] = 
			#rm(z, losses, gains, bias)
	
}


annotate_children <- function(this_node, this_tree, this_value) { #, this_root_node) {
	cat(paste0("Entry: ", this_node, " ", this_value, "\n"))
	# Find the children of this node, accumulate their net delta, and annotate them recursively
	for (this_child in find_children_ids(this_node, this_tree)) {
		# Calculate the delta from this_node to child_node by lookup from the edge_changes table
		if (this_child <= length(this_tree$tip.label)) {
			# child is a tip node
			this_delta = edge_changes[(edge_changes$ancestor == root_node) & (edge_changes$derived == this_tree$tip.label[this_child]),]$gains - edge_changes[(edge_changes$ancestor == root_node) & (edge_changes$derived == this_tree$tip.label[this_child]),]$losses		
		} else {
			# child is an internal node
			this_delta = edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_child),]$gains - edge_changes[(edge_changes$ancestor == this_node) & (edge_changes$derived == this_child),]$losses
		}
		cat(paste0(this_node, " ", this_child, " ", this_value, " ", this_delta, "\n"))
		
		# Pass the accumulated delta for the child node into the recurse
		this_tree = annotate_children(this_child, this_tree, this_value + this_delta)
	}
	
	# store the accumulated delta for this node in the node or tip label of the tree itself
	if(this_node>length(this_tree$tip.label)) {
		this_tree$node.label[this_node-length(this_tree$tip.label)] = this_value
	} else {
		this_tree$tip.label[this_node] = paste0(this_tree$tip.label[this_node], " ", as.character(this_value))
	}
	
	# pass the updated tree back down the line
	this_tree
}

annotated_tree = rerooted_tree
annotated_tree$node.label = rep(NA,annotated_tree$Nnode)
annotated_tree = annotate_children(this_node=root_node, this_tree=annotated_tree, this_value=0)
#annotated_tree = annotate_children(this_node=root_node, this_tree=annotated_tree, this_root_node=root_node)
#annotated_tree$edge.length = rerooted_tree$edge.length

plot(annotated_tree, cex=0.5)
nodelabels(annotated_tree$node.label, cex=0.5, frame="none")
write.tree(annotated_tree, file="rerooted_annotated_tree.nwk")





# Alternative plotting using ggtree package

library(ggtree)

#MRCA(rerooted_tree, tip=c("SRR342385", "SRR342381"))
#105
# Plot the basic tree
tree_plot = ggtree(rerooted_tree,) 
# Add node points
#tree_plot = tree_plot + geom_nodepoint()
# Add tip points
#tree_plot = tree_plot + geom_tippoint()
# Add the basic scale
tree_plot = tree_plot + geom_treescale()
# Add the fancy x axis
#tree_plot = tree_plot + theme_tree2()
# Add tip labels
tree_plot = tree_plot + geom_tiplab(size=2, align=TRUE, linesize=.5) 
# Add colour and label for MA lines clade
tree_plot = tree_plot + geom_cladelabel(node=142, label="MA lines", color="green") + geom_hilight(node=142, fill="green")
# Add colour and label for SALK lines clade
tree_plot = tree_plot + geom_cladelabel(node=116, label="SALK lines", color="pink") + geom_hilight(node=116, fill="pink")
tree_plot = tree_plot + geom_cladelabel(node=132, label="SALK lines", color="pink") + geom_hilight(node=132, fill="pink")
# Add internal node labels
#tree_plot = tree_plot + geom_text(aes(label=node), hjust=-.3)
# Polar coordinates
#tree_plot = tree_plot + coord_polar(theta='y')
tree_plot


 

mmaplot = function (p, fasta, offset = 0, width = 1, color = NULL, window = NULL, 
    bg_line = TRUE, height = 0.8) 
{
    if (missingArg(fasta)) {
        x <- NULL
    }
    else if (is(fasta, "DNAbin") || is(fasta, "AAbin")) {
        x <- fasta
    }
    else if (is(fasta, "character")) {
        x <- treeio::read.fasta(fasta)
    }
    else {
        x <- NULL
    }
    if (is.null(x) && is(p, "treedata") && length(p@tip_seq)) {
        x <- p@tip_seq
        p <- ggtree(p) + geom_tiplab()
    }
    if (is.null(x)) {
        stop("multiple sequence alignment is not available...\n-> check the parameter 'fasta'...")
    }
    x <- as.matrix(x)
    if (is.null(window)) {
        window <- c(1, ncol(x))
    }
    slice <- seq(window[1], window[2], by = 1)
    x <- x[, slice]
    seqs <- lapply(1:nrow(x), function(i) {
        seq <- as.vector(as.character(x[i, ]))
        seq[seq == "?"] <- "-"
        seq[seq == "*"] <- "-"
        seq[seq == " "] <- "-"
        return(seq)
    })
    names(seqs) <- labels(x)
    if (is.null(color)) {
        alphabet <- unlist(seqs) %>% unique
        alphabet <- alphabet[alphabet != "-"]
        color <- getCols(length(alphabet))
        names(color) <- alphabet
        color <- c(color, `-` = NA)
    }
    df <- p$data
    width <- width * (df$x %>% range %>% diff)/diff(window)
    df = df[df$isTip, ]
    start <- max(df$x) * 1.02 + offset
    seqs <- seqs[df$label[order(df$y)]]
    h <- ceiling(diff(range(df$y))/length(df$y))
    xmax <- start + seq_along(slice) * width
    xmin <- xmax - width
    y <- sort(df$y)
    ymin <- y - height/2 * h
    ymax <- y + height/2 * h
    from <- to <- NULL
    lines.df <- data.frame(from = min(xmin), to = max(xmax), 
        y = y)
    if (bg_line) {
        p <- p + geom_segment(data = lines.df, aes(x = from, 
            xend = to, y = y, yend = y), size = h * 0.2)
    }
    msa <- lapply(1:length(y), function(i) {
        data.frame(name = names(seqs)[i], xmin = xmin, xmax = xmax, 
            ymin = ymin[i], ymax = ymax[i], seq = seqs[[i]])
    })
    msa.df <- do.call("rbind", msa)
    p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, 
        ymax = ymax, fill = seq), data = msa.df, inherit.aes = FALSE) + 
        scale_fill_manual(values = color) + scale_x_continuous(breaks=1:length(slice))
    breaks <- hist(seq_along(slice), breaks = 10, plot = FALSE)$breaks
    pos <- start + breaks * width
    mapping <- data.frame(from = breaks + 1, to = pos)
    attr(p, "mapping") <- mapping
    return(p)
}



for (root_sites in c(0,1,2,3,4,5,6,7,8,9,10)) {
  # Plot the methylation pattern of each lowly methylated, long, gene-body methylated gene, in phylogenetic context
  #pdf(paste0("m_pattern_mCG_lt_0.05_sites_gt_50_genes.pdf"))
  pdf(paste0("m_pattern_",root_sites,"_site_root_single_internal_change_genes_P_142.pdf"))
  plot.new()
  #for (target_gene_ID in c("AT1G01030","AT5G54170","AT5G54180","AT5G54200","AT1G12930","AT5G09210","AT5G08780","AT5G08710","AT5G08130","AT5G07860","AT5G07850","AT1G08920","AT1G07770","AT5G05850","AT5G05820","AT1G02830","AT1G02960","AT1G03080","AT1G03140","AT1G03630","AT1G03980","AT1G04280","AT1G04500","AT1G04580","AT1G04750","AT1G04820","AT1G05290","AT1G05350","AT1G06310","AT1G07210","AT1G07550","AT1G07830","AT1G08920","AT1G09330","AT1G09380","AT1G09720","AT1G09840","AT1G10340","AT1G10430","AT4G17370","AT4G17730","AT4G17960","AT4G18230","AT4G18310","AT4G18330","AT4G18540","AT4G18800")) {
  
  #for (target_gene_ID in gff.genes[(gff.genes$average_methylation<0.05) & (gff.genes$CG_sites_count>50) & (gff.genes$m_class=="Gene-body Methylated"),"gene_ID"]) {
  for (target_gene_ID in GBM_gene_loci[queryHits(findOverlaps(makeGRangesFromDataFrame(df = GBM_gene_loci, start.field = "start_site", end.field = "end_site",seqnames.field = "Chromosome"), GBM_gene_loci.gr[(gene_internal_gain_more_sites==0) & (gene_initial_M_count==root_sites)])),"gene_ID"]) {
    target_gene_loci=gff.genes[gff.genes$gene_ID==target_gene_ID,c("V1","V4","V5")]
    colnames(target_gene_loci)=c("Chromosome", "start_site", "end_site")
    target_gene_loci$Chromosome=mixedCaseChr(target_gene_loci$Chromosome)
    target_gene_loci.gr=reduce(makeGRangesFromDataFrame(df=target_gene_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome"))
    
    # Extract a data matrix corresponding to the target gene loci
    tree_samples_target_gene_meth_status = tree_samples_meth_status[subjectHits(findOverlaps(target_gene_loci.gr, meth_status_sites.gr)),]
    
    
    # Export sites associated with a specific gene as a FASTA file
    ### Probably should amend this so it writes 5' to 3' per the gene, rather than per the chromosome coordinates
    ### Alternatively, add the gene logo to the plot to show where exons are
    this_fasta = paste0("4-Summary/",target_gene_ID,"_sites.fas")
    write.table(t(tree_samples_target_gene_meth_status)[4:ncol(tree_samples_target_gene_meth_status),], file=this_fasta, quote=FALSE, sep="", col.names=FALSE, na="?", row.names=paste0(">",colnames(tree_samples_target_gene_meth_status)[4:ncol(tree_samples_target_gene_meth_status)],"\n"))
    
    #plot.new()
    print(mmaplot(tree_plot, this_fasta, width=2.7, offset = 0.1, color = setNames(c("grey","red","purple","blue"),c("i","m","p","u"))))
    title(paste0(target_gene_ID, "  (",gff.genes[gff.genes$gene_ID==target_gene_ID,"V7"]," strand)"), outer=TRUE, adj=0.75, line=-1)
    legend("topright", legend=c("Methylated","Unmethylated","Partial","Missing data"),col=c("red","blue","purple","grey"), pch=c(15,15,15,15), cex=0.3, bg="white")
    
  }
  dev.off()
}

for (root_sites in c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) {
  # Count the number of internal and external losses and gains in this class of genes
  these_genes = (gene_initial_M_count==root_sites) & (gene_internal_gain_more_sites==0)
  cat(paste(root_sites, sum(these_genes), sum(gene_initial_M_count[these_genes]), sum(gene_internal_losses[these_genes]), sum(gene_internal_gains[these_genes]), sum(gene_external_losses[these_genes]), sum(gene_external_gains[these_genes]), "\n"))
}

# where root is JAY_1:
0 112 0 178 316 339 975 
1 189 189 430 395 656 1519 
2 295 590 908 638 1626 2746 
3 261 783 1084 621 1939 2719 
4 218 872 909 561 1698 2213 
5 161 805 773 434 1550 1784 
6 163 978 851 495 1695 1989 
7 124 868 799 387 1339 1474 
8 112 896 669 352 1249 1392 
9 102 918 747 376 1278 1418 
10 82 820 647 258 1065 1028 
11 63 693 508 215 948 876 
12 69 828 562 195 926 810 
13 56 728 511 201 791 768 
14 41 574 452 173 790 652 
15 64 960 633 239 1034 905 
16 45 720 491 138 771 698 
17 40 680 360 157 682 625 
18 32 576 423 158 613 573 
19 32 608 370 120 500 502 
20 31 620 421 121 601 419 

# Where root is 142 (ancestor of MA lines):
0 64 0 55 154 104 568 
1 186 186 410 403 605 1452 
2 285 570 856 678 1600 2800 
3 223 669 837 467 1571 2119 
4 187 748 889 491 1504 1822 
5 153 765 735 430 1451 1859 
6 119 714 606 313 1113 1189 
7 88 616 538 263 947 1022 
8 97 776 572 296 1076 1326 
9 75 675 548 246 896 893 
10 75 750 503 214 932 848 
11 45 495 346 142 627 529 
12 63 756 534 188 872 815 
13 41 533 349 131 578 500 
14 44 616 434 191 743 742 
15 35 525 332 117 506 418 
16 31 496 306 94 556 493 
17 46 782 430 157 682 648 
18 27 486 294 105 438 435 
19 35 665 409 156 667 486 
20 22 440 284 83 397 395 

These figures indicate that losses are more likely to be inherited than gains (>0.5 of internal changes are losses, proportion increases with increasing levels of methylation).
However, they also indicate that gains are more abundant than losses among 'new' changes (external) in lowly methylated genes (mCG<5).



# AT THIS POINT WE RUN THE SCRIPT 5-estimate_de_novo_rate_2019-09-12.R UP TO LINE 1050 (generate a masked segmentation model)
# in order to switch focus to the Schmitz data set, and calulate sites changed between generation 0 and 30


#levels(GBM_gene_loci.gr@seqnames@values) = mixedCaseChr(as.character(levels(GBM_gene_loci.gr@seqnames@values)))
levels(GBM_gene_loci.gr@seqnames@values) = toupper(as.character(levels(GBM_gene_loci.gr@seqnames@values)))
GBM_gene_loci.gr@seqinfo@seqnames = levels(GBM_gene_loci.gr@seqnames@values)

# Split GBM segments by segment length (no CG sites)
# Calculate probability of loss / site in cluster
#no_losses_observed/(no_segments*length_of_segment) by segment length group
GBM_segments.gr = masked_segmentation_model.gr[(masked_segmentation_model.gr$segment.mStatus=="GBM") & (!is.na(masked_segmentation_model.gr$nSites))]

#levels(GBM_segments.gr@seqnames@values) = mixedCaseChr(as.character(levels(GBM_segments.gr@seqnames@values)))
levels(GBM_segments.gr@seqnames@values) = toupper(as.character(levels(GBM_segments.gr@seqnames@values)))
GBM_segments.gr@seqinfo@seqnames = levels(GBM_segments.gr@seqnames@values)

for (change_class in c("losses","gains")) {
  cat(paste0(change_class,":\n"))
  previous_segment_length = 0 
  segment_length_groups = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30,35,40,45,50,100,500)
  for (this_segment_length in segment_length_groups) {
    # Find GBM segments overlapping the genes of interest
    #this_segment_group = GBM_segments.gr[(GBM_segments.gr$nSites>previous_segment_length) & (GBM_segments.gr$nSites<=this_segment_length)]
    #cat(paste(this_segment_length, length(this_segment_group), sum(this_segment_group$nSites), length(findOverlaps(clean_m_to_u_call_ranges, this_segment_group)), length(findOverlaps(clean_m_to_u_call_ranges, this_segment_group))/sum(this_segment_group$nSites)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[,], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)))
  
    if (change_class == "losses") {
      this_segment_group = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, GBM_gene_loci.gr[(gene_internal_gain_more_sites==0) & (gene_initial_M_count>previous_segment_length) & (gene_initial_M_count<=this_segment_length)]))]
      cat(paste(this_segment_length, length(this_segment_group), sum(this_segment_group$nSites), length(findOverlaps(clean_m_to_u_call_ranges, this_segment_group)), length(findOverlaps(clean_m_to_u_call_ranges, this_segment_group))/sum(this_segment_group$nSites)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[,], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)))
    } else {
      this_segment_group = GBM_gene_loci.gr[(gene_internal_gain_more_sites==0) & (gene_initial_M_count>previous_segment_length) & (gene_initial_M_count<=this_segment_length)]
      cat(paste(this_segment_length, length(this_segment_group), sum(this_segment_group$nSites), length(findOverlaps(clean_u_to_m_call_ranges, this_segment_group)), length(findOverlaps(clean_u_to_m_call_ranges, this_segment_group))/sum(this_segment_group$nSites)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[,], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)))
    }
  
  
    # above figures are related to length of segment.  instead let's consider no. M sites in segment:
    #length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group))
  
    for (this_line_no in 1:length(gen3_30)) {
      this_line = gen3_30[this_line_no]
    
      # find the subset of m_to_u and u_to_m which are in the target zone - should speed up the overlap searches later
      #line_u_to_m_target_ranges = line_u_to_m_call_ranges[[this_line]][queryHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
      #line_m_to_u_target_ranges = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
    
      if (change_class == "losses") {
        line_losses = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], this_segment_group))]
        line_prop_losses = length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), line_losses))
        cat(paste0(" ",length(line_losses), " ", line_prop_losses))
      } else {
        line_gains = line_u_to_m_call_ranges[[this_line]][queryHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], this_segment_group))]
        line_prop_gains = length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), line_gains))
        cat(paste0(" ",length(line_gains), " ", line_prop_gains))
      # line_losses and line_prop_losses are always coming out the same, but included just in case there is a problem
      }
    }
    cat("\n")
    previous_segment_length=this_segment_length
  }
}

losses 
1 165 332 26 0.0783132530120482 121 333 13 13 11 11 18 18 10 10
2 280 862 27 0.031322505800464 352 862 17 17 20 20 17 17 8 8
3 256 941 37 0.0393198724760893 449 943 32 32 18 18 18 18 9 9
4 205 1074 36 0.0335195530726257 441 1074 23 23 20 20 17 17 8 8
5 177 1152 43 0.0373263888888889 533 1152 21 21 24 24 23 23 13 13
6 140 1184 38 0.0320945945945946 574 1186 22 22 20 20 19 19 13 13
7 105 956 40 0.0418410041841004 497 956 16 16 8 8 21 21 18 18
8 135 1408 63 0.0447443181818182 708 1409 23 23 34 34 31 31 24 24
9 97 976 52 0.0532786885245902 544 976 22 22 21 21 16 16 19 19
10 86 1087 57 0.0524379024839006 588 1087 24 24 31 31 26 26 12 12
11 55 780 24 0.0307692307692308 437 780 14 14 7 7 18 18 8 8
12 81 1279 57 0.0445660672400313 697 1279 34 34 32 32 24 24 12 12
13 48 759 32 0.0421607378129117 421 759 7 7 14 14 19 19 6 6
14 56 1086 51 0.0469613259668508 571 1086 21 21 15 15 20 20 18 18
15 51 1011 43 0.0425321463897132 552 1011 18 18 15 15 19 19 19 19
20 194 4060 209 0.0514778325123153 2407 4060 82 82 85 85 93 93 64 64
25 101 2508 111 0.0442583732057416 1511 2508 46 46 40 40 48 48 24 24
30 72 2345 107 0.0456289978678038 1457 2345 57 57 31 31 46 46 34 34
35 40 1422 70 0.0492264416315049 900 1422 32 32 31 31 26 26 22 22
40 25 1097 46 0.0419325432999088 677 1098 20 20 28 28 26 26 9 9
45 9 456 17 0.037280701754386 311 456 12 12 5 5 10 10 3 3
50 10 582 33 0.0567010309278351 412 582 14 14 13 13 17 17 10 10
100 18 1267 73 0.0576164167324388 881 1267 26 26 21 21 25 25 23 23
500 1 151 12 0.0794701986754967 109 151 4 4 1 1 5 5 3 3
gains 
1 220 0 40 Inf 9180 10396 14 14 19 19 13 13 18 18
2 321 0 75 Inf 13256 15342 29 29 37 37 32 32 29 29
3 263 0 56 Inf 11381 13532 26 26 33 33 25 25 20 20
4 214 0 47 Inf 8088 9805 17 17 26 26 23 23 21 21
5 177 0 51 Inf 7843 9469 26 26 30 30 25 25 25 25
6 143 0 48 Inf 6467 7919 17 17 19 19 23 23 15 15
7 106 0 35 Inf 4779 6129 16 16 18 18 8 8 10 10
8 124 0 53 Inf 5722 7319 21 21 25 25 21 21 25 25
9 89 0 36 Inf 4275 5438 17 17 13 13 17 17 9 9
10 84 0 28 Inf 3687 4831 8 8 10 10 15 15 13 13
11 54 0 35 Inf 2354 3144 16 16 15 15 11 11 13 13
12 83 0 32 Inf 3771 5180 17 17 16 16 14 14 12 12
13 47 0 23 Inf 2031 2820 17 17 11 11 8 8 6 6
14 54 0 23 Inf 2645 3697 12 12 7 7 9 9 8 8
15 49 0 21 Inf 2133 3077 9 9 13 13 11 11 12 12
20 193 0 103 Inf 8710 12885 54 54 50 50 45 45 36 36
25 100 0 55 Inf 4611 7291 39 39 23 23 18 18 29 29
30 83 0 39 Inf 3563 6195 26 26 26 26 14 14 14 14
35 39 0 36 Inf 1943 3432 20 20 18 18 19 19 10 10
40 24 0 11 Inf 1228 2305 5 5 10 10 5 5 4 4
45 10 0 9 Inf 360 830 2 2 3 3 4 4 4 4
50 11 0 4 Inf 526 1092 1 1 5 5 2 2 2 2
100 18 0 16 Inf 913 2157 11 11 12 12 10 10 7 7
500 2 0 3 Inf 109 348 1 1 1 1 1 1 0 0
# These data show that rates of loss per M site are greatest (3 times greater) in genes with only a single mCG, in comparison to genes with more than one mCG
# These data also show that rates of gain per U site increase steadily with more mCG sites per gene.  The first site make the most difference to the rate of mCG gain, with the second site contributing 2/3 as much and subsequent site adding substantially less.


# Using the 'buffer' methodology from 5-estimate_de_novo_rate_2019-09-12.R:
  
# iterate over buffer_length and meth_quartile
no_stat_entries = 0				
segments_finished = 0

for (buffer_length in c(10)) {
  
  for (gbm_segment_length in 0:10) {
    cat("buffer length:",buffer_length,"\n")
    # meth_quartile 0 means do the unmethylated genes; 1-k means do the relevant subset of gene-body methylated genes by length of GBM segment
    
    
    pdf(paste0("m_pattern_gbM_segments_length_P_",gbm_segment_length,".pdf"))
    plot.new()
    
    
    
    # identify segments of interest
    #target_segments = unmethylated_segments∩gene-body_methylated_genes
    if (gbm_segment_length==0) {
      # Do unmethylated genes
      # Make segments to represent relevant target genes
      #unmethylated_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=na.omit(gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))
      target_segments.gr = reduce(makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Unmethylated") ,c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))[1:1000]
      #target_segments.gr = intersect(UMR_segments.gr, unmethylated_gene_loci.gr, ignore.strand=TRUE)
      # 10,416 segments
      
      # make enough space in segment_change_stats table for one entry per target segment per line
      no_stat_entries = length(target_segments.gr) * length(gen3_30)
      
      if(segments_finished==0) {
        segment_change_stats = cbind("buffer_length" = rep(buffer_length, no_stat_entries), "gbm_segment_length" = rep(gbm_segment_length, no_stat_entries), "line" = rep(NA, no_stat_entries), "gene_ID" = rep(NA, no_stat_entries), "segment_id"=seq(from=1, to=no_stat_entries, by=1), "segment_CG_sites" = rep(NA, no_stat_entries), "L_5_prime_length" = rep(NA, no_stat_entries), "core_length" = rep(NA, no_stat_entries),  "L_3_prime_length" = rep(NA, no_stat_entries),  "L_5_prime_length_bp" = rep(NA, no_stat_entries), "core_length_bp" = rep(NA, no_stat_entries), "L_3_prime_length_bp" = rep(NA, no_stat_entries), "L_5_prime_M_U" = rep(NA, no_stat_entries), "L_5_prime_U_M" = rep(NA, no_stat_entries), "core_M_U" = rep(NA, no_stat_entries), "core_U_M" = rep(NA, no_stat_entries), "L_3_prime_M_U" = rep(NA, no_stat_entries), "L_3_prime_U_M" = rep(NA, no_stat_entries))
      } else {
        segment_change_stats = rbind(segment_change_stats, cbind("buffer_length" = rep(buffer_length, no_stat_entries), "gbm_segment_length" = rep(gbm_segment_length, no_stat_entries), "line" = rep(NA, no_stat_entries), "gene_ID" = rep(NA, no_stat_entries), "segment_id"=seq(from=1, to=no_stat_entries, by=1), "segment_CG_sites" = rep(NA, no_stat_entries), "L_5_prime_length" = rep(NA, no_stat_entries), "core_length" = rep(NA, no_stat_entries),  "L_3_prime_length" = rep(NA, no_stat_entries),  "L_5_prime_length_bp" = rep(NA, no_stat_entries), "core_length_bp" = rep(NA, no_stat_entries), "L_3_prime_length_bp" = rep(NA, no_stat_entries), "L_5_prime_M_U" = rep(NA, no_stat_entries), "L_5_prime_U_M" = rep(NA, no_stat_entries), "core_M_U" = rep(NA, no_stat_entries), "core_U_M" = rep(NA, no_stat_entries), "L_3_prime_M_U" = rep(NA, no_stat_entries), "L_3_prime_U_M" = rep(NA, no_stat_entries)))
        
      }
    } else {
      # Do gene-body methylated genes
      #target_segments.gr = reduce(makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$gene_ID %in% names(table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[!is.na(segmentation_model.gr$nSites)&(segmentation_model.gr$nSites==gbm_segment_length)&(segmentation_model.gr$segment.mStatus=="GBM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))]))) & !(gff.genes$gene_ID %in% names(table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[(segmentation_model.gr$segment.mStatus=="TEM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))]))) & (gff.genes$gene_ID %in% names(table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[(segmentation_model.gr$segment.mStatus=="GBM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))]))[table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[(segmentation_model.gr$segment.mStatus=="GBM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))])==1]),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))
      #target_segments.gr = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, GBM_gene_loci.gr[(gene_internal_gain_more_sites==0) & (gene_initial_M_count==gbm_segment_length)]))]
      target_segments.gr = GBM_gene_loci.gr[(gene_internal_gain_more_sites==0) & (gene_initial_M_count==gbm_segment_length)]
      
      # Becker: 1502, 1097, 867, 697, 591, 504, 475, 428, 409, 337
      # Schmitz: 
      
      #target_segments.gr = intersect(UMR_segments.gr, makeGRangesFromDataFrame(
      #	df = gff.genes[(gff.genes$m_class == "Gene-body Methylated") &
      #					(gff.genes$average_methylation > meth_quartiles[meth_quartile]) &
      #					(gff.genes$average_methylation <= meth_quartiles[meth_quartile + 1]),],
      #	start.field = "V4",
      #	end.field = "V5",
      #	seqnames.field = "V1"
      #), ignore.strand=TRUE)
      # 29,324 segments
      
      # Accumulate previously finished count of segments, so rownumbers work later
      segments_finished = segments_finished + no_stat_entries
      no_stat_entries = length(target_segments.gr) * length(gen3_30)
      segment_change_stats = rbind(segment_change_stats, cbind("buffer_length" = rep(buffer_length, no_stat_entries), "gbm_segment_length" = rep(gbm_segment_length, no_stat_entries), "line" = rep(NA, no_stat_entries), "gene_ID" = rep(NA, no_stat_entries), "segment_id"=seq(from=1, to=no_stat_entries, by=1), "segment_CG_sites" = rep(NA, no_stat_entries), "L_5_prime_length" = rep(NA, no_stat_entries), "core_length" = rep(NA, no_stat_entries),  "L_3_prime_length" = rep(NA, no_stat_entries),  "L_5_prime_length_bp" = rep(NA, no_stat_entries), "core_length_bp" = rep(NA, no_stat_entries), "L_3_prime_length_bp" = rep(NA, no_stat_entries), "L_5_prime_M_U" = rep(NA, no_stat_entries), "L_5_prime_U_M" = rep(NA, no_stat_entries), "core_M_U" = rep(NA, no_stat_entries), "core_U_M" = rep(NA, no_stat_entries), "L_3_prime_M_U" = rep(NA, no_stat_entries), "L_3_prime_U_M" = rep(NA, no_stat_entries)))
      
    }
    cat("No parental Ms:", gbm_segment_length, " No. segments:",length(target_segments.gr),"\n")
    
    target_segments = as.data.frame(target_segments.gr)
    target_segments$segment_id = rownames(target_segments)
    target_segments_sites_olaps = findOverlaps(CG_site_ranges, target_segments.gr)
    target_segments$CG_sites_count = table(target_segments[subjectHits(target_segments_sites_olaps),]$segment_id)[target_segments$segment_id]
    
    valid_segments = 0
    for (target_no in 1:length(target_segments.gr)) {
      
      # Initialise buffer lengths to the values we are trying out
      L_5_prime_length = buffer_length
      L_3_prime_length = buffer_length
      
      target_segment = target_segments.gr[target_no]
      
      # what is strandedness of gene?
      target_gene = gff.genes[subjectHits(findOverlaps(target_segment, gene_ranges)),]
      if(nrow(target_gene) > 1) {
        # segment overlaps multiple genes: Assume strandedness of genes are opposite, and check both ends for alignment with gene 5'
        # for now, skip these genes and move on
      } else {
        #segment only overlaps one gene
        
        # Plot the gene's methylation pattern in phylogenetic context
        # Extract a data matrix corresponding to the target gene loci
        tree_samples_target_gene_meth_status = tree_samples_meth_status[subjectHits(findOverlaps(target_segment, meth_status_sites.gr)),]
        
        # Export sites associated with a specific gene as a FASTA file
        ### Probably should amend this so it writes 5' to 3' per the gene, rather than per the chromosome coordinates
        ### Alternatively, add the gene logo to the plot to show where exons are
        this_fasta = paste0("gene_sites/",target_gene$gene_ID,"_sites.fas")
        write.table(t(tree_samples_target_gene_meth_status)[4:ncol(tree_samples_target_gene_meth_status),], file=this_fasta, quote=FALSE, sep="", col.names=FALSE, na="?", row.names=paste0(">",colnames(tree_samples_target_gene_meth_status)[4:ncol(tree_samples_target_gene_meth_status)],"\n"))
        
        #plot.new()
        print(mmaplot(tree_plot, this_fasta, width=2.7, offset = 0.1, color = setNames(c("grey","red","purple","blue"),c("i","m","p","u"))))
        title(paste0(target_gene$gene_ID, "  (",gff.genes[gff.genes$gene_ID==target_gene$gene_ID,"V7"]," strand)"), outer=TRUE, adj=0.75, line=-1)
        legend("topright", legend=c("Methylated","Unmethylated","Partial","Missing data"),col=c("red","blue","purple","grey"), pch=c(15,15,15,15), cex=0.3, bg="white")
        
        
        
        
        # based on strandedness of gene, does target segment align with 5' end of gene (or 3')?
        # If so, adjust buffer zones to approximate 2kb/1kb end-of-gene buffers
        if(target_gene$V7=="+") {
          if (as.data.frame(target_segment)$start == target_gene$V4) {
            L_5_prime_length = extended_L_5_prime_length
          }
          if (as.data.frame(target_segment)$end == target_gene$V5) {
            L_3_prime_length = extended_L_3_prime_length
          }
        } else {
          if (as.data.frame(target_segment)$start == target_gene$V4) {
            L_5_prime_length = extended_L_3_prime_length
          }
          if (as.data.frame(target_segment)$end == target_gene$V5) {
            L_3_prime_length = extended_L_5_prime_length
          }
        }
        
        #cat(target_no, " ", gff.genes[subjectHits(findOverlaps(target_segment, gene_ranges)),]$gene_ID, " ") ##, "\n")
        
        # How many CG sites overlap the target region?
        segment_CG_sites = target_segments[target_no,]$CG_sites_count
        
        # Is the segment longer than the buffer zones?
        
        if (!is.na(segment_CG_sites) & (segment_CG_sites > (L_5_prime_length + L_3_prime_length))) {
          valid_segments = valid_segments + 1
          cat(paste(target_no, target_gene$gene_ID, "\n"))
          
          # Find start and end loci of core zone and buffer zones
          
          # These are the CG sites overlapping the target_segment
          #queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no])
          
          # Start position of 5' buffer 
          L_buffer_start = as.data.frame(CG_site_ranges[min(queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no]))])$start
          # End position of 5' buffer
          L_buffer_end = as.data.frame(CG_site_ranges[min(queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no]))+L_5_prime_length-1])$start+1
          # Start position of core 
          core_start = as.data.frame(CG_site_ranges[min(queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no]))+L_5_prime_length-1])$start+2
          # End position of core
          core_end = as.data.frame(CG_site_ranges[max(queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no]))+1-L_3_prime_length])$start-2
          # Start position of 3' buffer 
          R_buffer_start = as.data.frame(CG_site_ranges[max(queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no]))+1-L_3_prime_length])$start-1
          # End position of 3' buffer
          R_buffer_end = as.data.frame(CG_site_ranges[max(queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no]))])$end
          
          
          
          ### Break out gains and losses into separate events per line
          
          this_chrom = as.data.frame(CG_site_ranges[min(queryHits(target_segments_sites_olaps[subjectHits(target_segments_sites_olaps)==target_no]))+L_5_prime_length-1])$seqnames
          
          for (this_line_no in 1:length(gen3_30)) {
            this_line = gen3_30[this_line_no]
            
            # find the subset of m_to_u and u_to_m which are in the target zone - should speed up the overlap searches later
            line_u_to_m_target_ranges = line_u_to_m_call_ranges[[this_line]][queryHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
            line_m_to_u_target_ranges = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
            
            segment_change_stats[segments_finished + target_no + length(target_segments.gr) * (this_line_no - 1),] = c(buffer_length, gbm_segment_length, this_line, target_gene$gene_ID, target_no, segment_CG_sites, L_5_prime_length, segment_CG_sites - (L_5_prime_length + L_3_prime_length), L_3_prime_length, L_buffer_end-L_buffer_start+1, core_end - core_start + 1, R_buffer_end-R_buffer_start+1, length(subjectHits(findOverlaps(line_m_to_u_target_ranges, makeGRangesFromDataFrame(df = data.frame("start"=L_buffer_start, "end"=L_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_target_ranges, makeGRangesFromDataFrame(df = data.frame("start"=L_buffer_start, "end"=L_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_m_to_u_target_ranges, makeGRangesFromDataFrame(df = data.frame("start"=core_start, "end"=core_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_target_ranges, makeGRangesFromDataFrame(df = data.frame("start"=core_start, "end"=core_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_m_to_u_target_ranges, makeGRangesFromDataFrame(df = data.frame("start"=R_buffer_start, "end"=R_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_target_ranges, makeGRangesFromDataFrame(df = data.frame("start"=R_buffer_start, "end"=R_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))))
            
            #segment_change_stats[segments_finished + target_no + length(target_segments.gr) * (this_line_no - 1),] = c(buffer_length, gbm_segment_length, this_line, target_gene$gene_ID, target_no, segment_CG_sites, L_5_prime_length, segment_CG_sites - (L_5_prime_length + L_3_prime_length), L_3_prime_length, L_buffer_end-L_buffer_start+1, core_end - core_start + 1, R_buffer_end-R_buffer_start+1, length(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=L_buffer_start, "end"=L_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=L_buffer_start, "end"=L_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=core_start, "end"=core_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=core_start, "end"=core_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=R_buffer_start, "end"=R_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=R_buffer_start, "end"=R_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))))
            
            # commenting this out as it reiterates the slow overlap steps done in the previous line
            #cat(paste(buffer_length, gbm_segment_length, this_line, target_gene$gene_ID, target_no, segment_CG_sites, L_5_prime_length, segment_CG_sites - (L_5_prime_length + L_3_prime_length), L_3_prime_length, L_buffer_end-L_buffer_start+1, core_end - core_start + 1, R_buffer_end-R_buffer_start+1, length(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=L_buffer_start, "end"=L_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=L_buffer_start, "end"=L_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=core_start, "end"=core_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=core_start, "end"=core_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=R_buffer_start, "end"=R_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))), length(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], makeGRangesFromDataFrame(df = data.frame("start"=R_buffer_start, "end"=R_buffer_end, "chrom"=this_chrom), start.field = "start", end.field = "end", seqnames.field = "chrom")))),"\n"))
            #cat(paste(buffer_length, gbm_segment_length, this_line, target_gene$gene_ID, target_no, segment_CG_sites, L_5_prime_length, segment_CG_sites - (L_5_prime_length + L_3_prime_length), L_3_prime_length, L_buffer_end-L_buffer_start+1, core_end - core_start + 1, R_buffer_end-R_buffer_start+1, "\n"))
          }
        } #else {
        # Segment is shorter than sum of two buffer zones
        #}
      }
      
      #if target_segment aligns with gene 5’ end:
      
      #            L_5_prime = extended_L_5_prime
      
      #else
      
      #            L_5_prime = L_5_prime_length
      
      #L_3_prime = L_3_prime_length
      
      #record bp lengths of L_5_prime and L_3_prime
      
      #remove L_5_prime and L_3_prime from target_segment to create core_target_segment
      
      #record core_target_segment_length (W) in sites and bp
      
      #record No. gains (g) and losses (l) in core_target_segment
      
      #record No. gains (g) and losses (l) in L_5_prime and L_3_prime
      
      #Record length of L_3_prime and L_5_prime
      
      #Record any times where there is a change in the buffer region		
    }
    cat(paste0("Segments: ",nrow(target_segments), "  Valid segments: ", valid_segments,"\n"))
    dev.off()
  }
}

segment_change_stats = as.data.frame(segment_change_stats)

segment_change_stats$core_U_M = as.numeric(levels(segment_change_stats$core_U_M)[segment_change_stats$core_U_M])
segment_change_stats$core_M_U = as.numeric(levels(segment_change_stats$core_M_U)[segment_change_stats$core_M_U])
segment_change_stats$core_length = as.numeric(levels(segment_change_stats$core_length)[segment_change_stats$core_length])
segment_change_stats$core_length_bp = as.numeric(levels(segment_change_stats$core_length_bp)[segment_change_stats$core_length_bp])
segment_change_stats$buffer_length = as.numeric(levels(segment_change_stats$buffer_length)[segment_change_stats$buffer_length])

saveRDS(segment_change_stats, file=paste0(project_id,"_",meth_context,"_segment_change_stats_by_gbM_segment_length.rds"))

pdf(paste0(project_id,"_distribution_of_M-U_core_transitions_by_gbM_segment_length.pdf"))
# Plot histograms of lengths of core segments
# Exclude segments with no stats, with U-M transitions in buffer zone
ggplot(segment_change_stats[(!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),]) + geom_histogram(aes(x=core_length_bp), binwidth=10) + facet_grid(gbm_segment_length~buffer_length) + scale_y_continuous(trans='log10') + xlab("Core segment length (bp)") + ylab("Number of segments") + theme_minimal()

ggplot(segment_change_stats[(!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),]) + geom_histogram(aes(x=core_length), binwidth=10) + facet_grid(gbm_segment_length~buffer_length) + scale_y_continuous(trans='log10') + xlab("Core segment length (CG sites)") + ylab("Number of segments") + theme_minimal()

# Plot histograms of numbers of U-M transitions in core
ggplot(segment_change_stats[(!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),]) + geom_histogram(aes(x=core_U_M), binwidth=1) + facet_grid(gbm_segment_length~buffer_length) + scale_y_continuous(trans='log10') + xlab("No. U->M changes/segment in core") + ylab("Number of segments") + theme_minimal()

# Plot distribution of proportion of U-M transitions in core
ggplot(segment_change_stats[(!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),]) + geom_histogram(aes(x=core_U_M/core_length), binwidth=.01) + facet_grid(gbm_segment_length~buffer_length) + scale_y_continuous(trans='log10') + xlim(-0.01,0.25) + xlab("Proportion of core sites changing U->M") + ylab("Number of segments") + theme_minimal()
dev.off()

# Extract table of rate estimates for each line for each set of parameter values
plot_histograms = FALSE
for (this_buffer_length in unique(segment_change_stats$buffer_length)) {
  for (this_gbm_segment_length in unique(segment_change_stats$gbm_segment_length)) {
    for (this_line in gen3_30) {
      # Calculate proportion of sites gaining methylation
      changed_site_counts = hist(segment_change_stats[(segment_change_stats$buffer_length==this_buffer_length) & (segment_change_stats$gbm_segment_length==this_gbm_segment_length) & (segment_change_stats$line == this_line) & (!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),"core_U_M"], breaks=seq(from=-0.5, to=10.5, by=1), plot=plot_histograms)$counts
      
      # Calculate proportion of sites losing methylation
      #changed_site_counts = hist(segment_change_stats[(segment_change_stats$buffer_length==this_buffer_length) & (segment_change_stats$gbm_segment_length==this_gbm_segment_length) & (segment_change_stats$line == this_line) & (!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),"core_M_U"], breaks=seq(from=-0.5, to=10.5, by=1), plot=plot_histograms)$counts
      
      tot_changes = 0
      for (no_changes in 1:10) {
        tot_changes = tot_changes + no_changes * changed_site_counts[no_changes+1]
      }
      prop_changes = tot_changes/sum(segment_change_stats[(segment_change_stats$buffer_length==this_buffer_length) & (segment_change_stats$gbm_segment_length==this_gbm_segment_length) & (segment_change_stats$line == this_line) & (!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),"core_length"])
      #cat(sum(segment_change_stats[(segment_change_stats$buffer_length==this_buffer_length) & (segment_change_stats$gbm_segment_length==this_gbm_segment_length) & (segment_change_stats$line == this_line) & (!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),"core_length"]),"\n")
      cat(c(this_buffer_length, this_gbm_segment_length, this_line, prop_changes, changed_site_counts), "\n")
    }
  }
}














# Plot tree_plot with a heatmap
#tree_plot_combo = (tree_plot + scale_y_continuous(expand=c(0, 0.3))) %>% gheatmap(genotype, offset=4, width=0.5, colnames=FALSE) %>% scale_x_ggtree() 
#tree_plot_combo + theme(legend.position="right")


# Need to check for SNPs and control 

# Simulate how low can you go on coverage and still get acceptable phylogenetic placement



# Analyse small segments:


# read in gff.genes
gff.genes = readRDS("../Jay-SMP_variation/5-analysis/SRA035939/Draft Model results/SRA035939_CG_gff.genes.rds")
gff.transposons = readRDS("../Jay-SMP_variation/5-analysis/SRA035939/Draft Model results/SRA035939_CG_gff.transposons.rds")

gff.genes$V1 = mixedCaseChr(gff.genes$V1)
gff.transposons$V1 = mixedCaseChr(gff.transposons$V1)



library(GenomicRanges)

# Make a GRanges for gene space
gene_ranges=makeGRangesFromDataFrame(df = gff.genes, start.field = "V4", end.field = "V5",seqnames.field = "V1")
# Make a GRanges for transposon space
transposon_ranges=makeGRangesFromDataFrame(df = gff.transposons, start.field = "V4", end.field = "V5",seqnames.field = "V1")

#  Make a GRanges for all CG sites
CG_site_ranges =makeGRangesFromDataFrame(df = tree_samples_meth_status, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")

# Find the overlaps
olaps_CG_sites = findOverlaps(CG_site_ranges, gene_ranges)
tolaps_CG_sites = findOverlaps(CG_site_ranges, transposon_ranges)


# Load in the draft segmentation model prepared using Schmitz gen 3 data only:
segmentation_model.gr = readRDS(paste0("../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_segmentation_model_draft3.rds"))
as(segmentation_model.gr,"GRanges")

	

# Make table of genes with M counts per sample (do 1 sample as an example)
### this example doesn't really go anywhere - delete it?

SRR342385M_ranges = makeGRangesFromDataFrame(df = tree_samples_meth_status[tree_samples_meth_status$SRR342385=="M",], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
olaps_gene_SRR342385M = findOverlaps(SRR342385M_ranges, gene_ranges)
gff.genes$SRR342385M_count = table(gff.genes[subjectHits(olaps_gene_SRR342385M),]$gene_ID)[gff.genes$gene_ID]

# Manual review:
View(gff.genes[!is.na(gff.genes$SRR342385M_count),])

ggplot(gff.genes[(gff.genes$Gen1M_count) == (gff.genes$Gen2M_count),c(1,11,12,13,19,21,22,23,24,25,26,27)]) + geom_boxplot(aes(x=Gen1M_count, y=Gen3M_count))
 
# replace NAs with 0:
gff.genes[is.na(gff.genes$Gen1M_count),]$Gen1M_count = 0
gff.genes[is.na(gff.genes$Gen2M_count),]$Gen2M_count = 0
gff.genes[is.na(gff.genes$Gen3M_count),]$Gen3M_count = 0
gff.genes[is.na(gff.genes$Gen4M_count),]$Gen4M_count = 0
gff.genes[is.na(gff.genes$Gen5M_count),]$Gen5M_count = 0
gff.genes[is.na(gff.genes$Gen6M_count),]$Gen6M_count = 0



# Find overlaps between GBM segments and genes
gene_GBM_segment_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM"),], gene_ranges)

# Pick out genes overlapping a single segment 
#gene_GBM_segment_count = table(gff.genes[subjectHits(gene_GBM_segment_olaps),]$gene_ID)[gff.genes$gene_ID]
#gene_single_GBM_segment = (table(gff.genes[subjectHits(gene_GBM_segment_olaps),]$gene_ID)[gff.genes$gene_ID] == 1)

single_GBM_segment_genes = gff.genes[(table(gff.genes[subjectHits(gene_GBM_segment_olaps),]$gene_ID)[gff.genes$gene_ID] == 1),][!is.na(gff.genes[(table(gff.genes[subjectHits(gene_GBM_segment_olaps),]$gene_ID)[gff.genes$gene_ID] == 1),"gene_ID"]),]

#gene_single_GBM_segment_ranges = makeGRangesFromDataFrame(df = gff.genes[gene_single_GBM_segment,][!is.na(gff.genes[gene_single_GBM_segment,"gene_ID"]),], start.field = "V4", end.field = "V5",seqnames.field = "V1")

single_GBM_segment_ranges = makeGRangesFromDataFrame(df = single_GBM_segment_genes, start.field = "V4", end.field = "V5",seqnames.field = "V1")

# Identify genes overlapping segments of length 1 CpG site, 2 sites, 3 sites etc
gene_GBM_segment_1_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==1),], single_GBM_segment_ranges)
gene_GBM_segment_2_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==2),], single_GBM_segment_ranges)
gene_GBM_segment_3_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==3),], single_GBM_segment_ranges)
gene_GBM_segment_4_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==4),], single_GBM_segment_ranges)
gene_GBM_segment_5_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==5),], single_GBM_segment_ranges)
gene_GBM_segment_6_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==6),], single_GBM_segment_ranges)
gene_GBM_segment_7_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==7),], single_GBM_segment_ranges)
gene_GBM_segment_8_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==8),], single_GBM_segment_ranges)
gene_GBM_segment_9_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==9),], single_GBM_segment_ranges)
gene_GBM_segment_10_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==10),], single_GBM_segment_ranges)
gene_GBM_segment_11_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==11),], single_GBM_segment_ranges)
gene_GBM_segment_12_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==12),], single_GBM_segment_ranges)
gene_GBM_segment_13_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==13),], single_GBM_segment_ranges)
gene_GBM_segment_14_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==14),], single_GBM_segment_ranges)
gene_GBM_segment_15_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==15),], single_GBM_segment_ranges)
gene_GBM_segment_16_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==16),], single_GBM_segment_ranges)
gene_GBM_segment_17_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==17),], single_GBM_segment_ranges)
gene_GBM_segment_18_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==18),], single_GBM_segment_ranges)
gene_GBM_segment_19_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==19),], single_GBM_segment_ranges)
gene_GBM_segment_20_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==20),], single_GBM_segment_ranges)
gene_GBM_segment_21_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==21),], single_GBM_segment_ranges)
gene_GBM_segment_22_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==22),], single_GBM_segment_ranges)
gene_GBM_segment_23_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==23),], single_GBM_segment_ranges)
gene_GBM_segment_24_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==24),], single_GBM_segment_ranges)
gene_GBM_segment_25_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==25),], single_GBM_segment_ranges)
gene_GBM_segment_26_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==26),], single_GBM_segment_ranges)
gene_GBM_segment_27_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==27),], single_GBM_segment_ranges)
gene_GBM_segment_28_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==28),], single_GBM_segment_ranges)
gene_GBM_segment_29_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==29),], single_GBM_segment_ranges)
gene_GBM_segment_30_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==30),], single_GBM_segment_ranges)
gene_GBM_segment_31_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==31),], single_GBM_segment_ranges)
gene_GBM_segment_32_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==32),], single_GBM_segment_ranges)
gene_GBM_segment_33_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==33),], single_GBM_segment_ranges)
gene_GBM_segment_34_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==34),], single_GBM_segment_ranges)
gene_GBM_segment_35_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==35),], single_GBM_segment_ranges)
gene_GBM_segment_36_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==36),], single_GBM_segment_ranges)
gene_GBM_segment_37_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==37),], single_GBM_segment_ranges)
gene_GBM_segment_38_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==38),], single_GBM_segment_ranges)
gene_GBM_segment_39_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==39),], single_GBM_segment_ranges)
gene_GBM_segment_40_olaps = findOverlaps(segmentation_model.gr[(as.data.frame(segmentation_model.gr)$segment.mStatus=="GBM") & (!is.na(as.data.frame(segmentation_model.gr)$nSites)) & (as.data.frame(segmentation_model.gr)$nSites==40),], single_GBM_segment_ranges)

gene_GBM_segment_1_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_1_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_2_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_2_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_3_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_3_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_4_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_4_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_5_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_5_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_6_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_6_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_7_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_7_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_8_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_8_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_9_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_9_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_10_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_10_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_11_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_11_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_12_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_12_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_13_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_13_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_14_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_14_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_15_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_15_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_16_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_16_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_17_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_17_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_18_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_18_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_19_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_19_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_20_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_20_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_21_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_21_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_22_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_22_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_23_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_23_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_24_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_24_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_25_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_25_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_26_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_26_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_27_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_27_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_28_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_28_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_29_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_29_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_30_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_30_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_31_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_31_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_32_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_32_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_33_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_33_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_34_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_34_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_35_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_35_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_36_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_36_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_37_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_37_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_38_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_38_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_39_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_39_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]
gene_GBM_segment_40_count = table(single_GBM_segment_genes[subjectHits(gene_GBM_segment_40_olaps),]$gene_ID)[single_GBM_segment_genes$gene_ID]

cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_1_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_1_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_2_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_2_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_3_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_3_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_4_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_4_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_5_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_5_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_6_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_6_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_7_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_7_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_8_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_8_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_9_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_9_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_10_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_10_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_11_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_11_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_12_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_12_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_13_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_13_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_14_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_14_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_15_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_15_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_16_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_16_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_17_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_17_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_18_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_18_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_19_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_19_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_20_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_20_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_21_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_21_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_22_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_22_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_23_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_23_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_24_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_24_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_25_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_25_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_26_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_26_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_27_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_27_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_28_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_28_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_29_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_29_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_30_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_30_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_31_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_31_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_32_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_32_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_33_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_33_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_34_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_34_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_35_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_35_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_36_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_36_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_37_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_37_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_38_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_38_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_39_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_39_olaps),]$u_to_m_count),"\n"))
cat(paste0(mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_40_olaps),]$m_to_u_count)," ", mean(single_GBM_segment_genes[subjectHits(gene_GBM_segment_40_olaps),]$u_to_m_count),"\n"))



