#!/usr/bin/env Rscript

# non-CG_methylation.R
#
# Jay Moore, John Innes Centre, Norwich
# jonathan.moore@jic.ac.uk
#
# 28/10/2019
#
# Version 1 - 1001 methylomes data analysis, adapted from Schmitz/Becker project - methylation has previously been called for all samples.  Integrate them into a single matrix
#
# Change log
#


# Data source (project_id)
#  .methratio2 files from bsmap


# Summary of functions
#
# Done:
#
# Underway:
#
# To do:
# 
#
#
# Results:
#
# 

if(!require(optparse)){
	install.packages("optparse")
	library(optparse)
}

library("optparse")

option_list = list(
	make_option(c("-p", "--project"), action="store", default="m1001", type='character',
              help="project code for location of mapping output"),
	make_option(c("-s", "--sample"), action="store", default=NA, type='character',
              help="project code for location of mapping output"),
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
    
	cat("Sample: ")
    cat(opt$sample)
    cat("\n")
    
    cat("Context: ")
    cat(opt$context)
    cat("\n")

    cat("Action: ")
    cat(opt$action)
    cat("\n")
}

# project_id, sample_id and meth_context together are used for finding raw data and naming output files
project_id=NULL
if (is.na(opt$project)) {
	# No project defined - set a default
	project_id="m1001"
	#project_id="memory_line_2"
} else {
	project_id=opt$project
}
	
sample_id=NULL
if (is.na(opt$sample)) {
	# No sample defined - quit
	stop("No sample specified")
	#project_id="memory_line_2"
} else {
	sample_id=opt$sample
}
	
meth_context=NULL
if (is.na(opt$context)) {
	# No context defined - set a default
	meth_context="CG"
	#meth_context="CHG"
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

#### This part installs packages so will be slow the first time per platform

# Sometimes rlang is needed for bioconductor to instal. This problem seems to have gone away now...
#if(!require(rlang)){
#	install.packages("rlang")
#	library(rlang)
#}

	
#### This part sets up libraries and contextual information and metadata about the project


# Platform-specific stuff:
#  Base of path to project
#  Where to find bioconductor (packages need to be pre-installed in R using singularity if running on the cluster)
on_cluster = FALSE
pathroot = ""
if (.Platform$OS.type=="windows") {
  pathroot="X:/Daniel-Zilberman"
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
  #biocLite(c("GenomicRanges", "BSgenome", "BSgenome.Athaliana.TAIR.TAIR9", "MethylSeekR", "karyoploteR"))
  BiocManager::install(c("GenomicRanges", "BSgenome", "BSgenome.Athaliana.TAIR.TAIR9", "MethylSeekR", "karyoploteR", "ggtree"))
}

#### This part installs packages so will be slow the first time per platform

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
if(!require(sqldf)){
  install.packages("sqldf")
  library(sqldf)
}

library(data.table)
library(stringr)
library(ggplot2)
library(plyr)
#library(ggpubr)
library(sqldf)


# Source directory containing alignments from bs_sequel
source_dir = paste0(pathroot,"/Projects/Jay-1001_methylomes/3-alignments/")

# Working directory where outputs will be directed to
setwd(dir = paste0(pathroot,"/Projects/Jay-1001_methylomes/5-analysis/"))

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

# H3K4me1 Plant DNase I hypersensitivity data. Downloaded from http://plantdhs.org/Download (Jiming Jiang lab) on 05/12/2017
# For H3K4me1 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.RDS")
reference_H3K4me1 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_buds_H3K4me1.RDS")

# H3K27me3 Plant DNase I hypersensitivity data. Downloaded from http://plantdhs.org/Download (Jiming Jiang lab) on 05/12/2017
# For H3K27me3 data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.RDS")
reference_H3K27me3 = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/Ath_leaf_H3K27me3.RDS")

# H2AW data from Jaemyung Choi
reference_H2AW = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Histone_modifications/h2aw.w50.gff")





# Sample info from Kawakatsu et al, 2016 Table S2
# http://neomorph.salk.edu/publications/Kawakatsu_Cell.pdf
# This has useful metadata about the samples as published in the paper
sample_info = read.table(file="../0-reference/1-s2.0-S0092867416308522-mmc3.txt", sep="\t", header=TRUE)



# Sample list is a list of the samples we have sequence data for
# fastq list is the cor4respondence between fastq files we have aligned, and samples
# Global Arabidopsis accessions from Ecker study (927 of)
sample_list1 = read.table(paste0("../0-raw_data/","GSE43857.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=43857
fastq_list1 = read.table(paste0("../0-raw_data/","ENA_SRA065807.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/SRA065807

# sample_list1$Title is a string containing 3 separate data items:
# sample_info$Name (sample_info$Ecotype_id)<-Replicate_id>
# Let's split them up
sample_list1$Name = str_split_fixed(sample_list1$Title, "\\(",2)[,1]
sample_list1$Ecotype_id = as.numeric(str_split_fixed(str_split_fixed(sample_list1$Title, "\\(",2)[,2], "\\)", 2)[,1])
sample_list1$Replicate_id = str_split_fixed(str_split_fixed(sample_list1$Title, "\\(",2)[,2], "-", 2)[,2]
# sample_list1$Name is sample_info$Name
# Merge sample_list1 with sample_info by Name

# Swedish Arabidopsis accessions from separate study (284 of)
sample_list2 = read.table(paste0("../0-raw_data/","GSE54292.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=54292
fastq_list2 = read.table(paste0("../0-raw_data/","ENA_SRP035593.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/PRJNA236110
# sample_list2$Title is sample_info$Sample, and Name is not present
#sample_list2$Sample = sample_list$Title

# Merge sample_list2 with sample_info by Sample


# Put the two sets of accessions together
# This step never quite worked, due to the differences in annotation between the different labs, so did the merge manually in Excel and import this in the block of code below
#sample_list = rbind.data.frame(sample_list1, sample_list2)
fastq_list = rbind.data.frame(fastq_list1, fastq_list2)

rm(sample_list1, fastq_list1, sample_list2, fastq_list2)


# The above sample_list routine was somewhat unsatisfactory.  Instead, we will read in the synthesis created manually in Excel, using the above resources as a starting point:
sample_list = read.table(paste0("../0-raw_data/","sample_info_synthesis.txt"), sep="\t", header=TRUE)
sample_list$Name = str_split_fixed(sample_list$Title, "\\(",2)[,1]
#sample_list1$Ecotype_id = as.numeric(str_split_fixed(str_split_fixed(sample_list1$Title, "\\(",2)[,2], "\\)", 2)[,1])
sample_list$Ecotype_id = as.numeric(str_split_fixed(sample_list$Sample, "_", 2)[,1])


### this merge fails - sample_info has no Title
sample_info_merged = merge(sample_info, sample_list, by=c("Sample"), all=TRUE)
### Maybe it should merge on c(Name, Ecotype_id), but not certain whether this is unique - probably not

sample_files = data.frame()
for (this_accession in sample_list$SRA.Accession) {
	these_files = fastq_list[fastq_list$experiment_accession == this_accession,"fastq_ftp"]
	for (this_file in these_files) {
		some_files=strsplit(this_file,";")
		sample_files = rbind.data.frame(sample_files, c("Accession" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Accession"]),  "Title" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Title"]), "SRA_Accession" = this_accession, "FASTQ_file" = some_files[[1]][[1]]), stringsAsFactors=FALSE)
		if (length(some_files[[1]]) > 1) {
			sample_files = rbind.data.frame(sample_files, c("Accession" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Accession"]), "Title" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Title"]), "SRA_Accession" = this_accession, "FASTQ_file" = some_files[[1]][[2]]), stringsAsFactors=FALSE)
			# it looks like a lot of the paired-end files are represented as 'SINGLE' but should be 'PAIRED' - replace where necessary
			if (fastq_list[fastq_list$fastq_ftp == this_file,"library_layout"] == "SINGLE") {
				fastq_list[fastq_list$fastq_ftp == this_file,"library_layout"] = "PAIRED"
			}
		} else {
			#cat(fastq_list[fastq_list$experiment_accession == this_accession,"library_layout"])
		}
	}
}
colnames(sample_files) = c("Accession", "Title", "SRA_Accession", "FASTQ_file")


# Load the reference genome annotation, and anotate each gene and transposon based on the CG methylome
    library(GenomicRanges)	
	
	#gff.genes=readRDS(paste0(reference_dir,"SRA035939_CG_gff.genes.rds"))
	#gff.transposons=readRDS(paste0(reference_dir,"SRA035939_CG_gff.transposons.rds"))

	# Load the annotation
	annot_gff = read.delim(reference_gff, header=F, comment.char="#")	
	#gff.exons = read.delim(reference_exons, header=F, comment.char="#")
	#gff.introns = read.delim(reference_introns, header=F, comment.char="#")
	#gff.5UTR = read.delim(reference_5UTR, header=F, comment.char="#")
	#gff.3UTR = read.delim(reference_3UTR, header=F, comment.char="#")

	# Grab the portion relating to genes
	gff.genes = annot_gff[annot_gff[,3]=="gene",]
	# Grab the portion relating to transposons
	gff.transposons = annot_gff[annot_gff[,3]=="transposable_element",]

	rm(annot_gff)
	
	# Convert chromosome names to uppercase to match previous objects
	#gff.genes$V1=toupper(gff.genes$V1)
	#gff.transposons$V1=toupper(gff.transposons$V1)
	#gff.exons$V1=toupper(gff.exons$V1)
	#gff.introns$V1=toupper(gff.introns$V1)
	#gff.5UTR$V1=toupper(gff.5UTR$V1)
	#gff.3UTR$V1=toupper(gff.3UTR$V1)
	
	# The following line sometimes fails if stringi has not been installed correctly.  Not sure why.  Fix is to install.packages("stringi") and to allow RStudio to restart R.
	# Grab gene IDs from descriptive field
	gff.genes$gene_ID=str_split_fixed(str_split_fixed(gff.genes$V9, ';',3)[,1],"=",2)[,2]
	gff.transposons$gene_ID=str_split_fixed(str_split_fixed(gff.transposons$V9, ';',3)[,1],"=",2)[,2]
	#gff.exons$gene_ID=str_split_fixed(str_split_fixed(gff.exons$V9, ';',3)[,1],"=",2)[,2]
	#gff.introns$gene_ID=str_split_fixed(str_split_fixed(gff.introns$V9, ';',3)[,1],"=",2)[,2]
	#gff.5UTR$gene_ID=str_split_fixed(str_split_fixed(gff.5UTR$V9, ';',3)[,1],"=",2)[,2]
	#gff.3UTR$gene_ID=str_split_fixed(str_split_fixed(gff.3UTR$V9, ';',3)[,1],"=",2)[,2]

	#Gene names are tricky - need to find the element containing "symbol="
	gff.genes=within(gff.genes, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,2],"=",2)[,2]})
	gff.transposons=within(gff.transposons, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,3],"=",2)[,2]})

	# Fix the gene names in exon, intron and UTRs
	#gff.exons$exon_ID = gff.exons$gene_ID
	#gff.exons$gene_ID = substr(gff.exons$gene_ID,1,9)
	#gff.introns$intron_ID = gff.introns$gene_ID
	#gff.introns$gene_ID = substr(gff.introns$gene_ID,1,9)
	#gff.5UTR$UTR_ID = gff.5UTR$gene_ID
	#gff.5UTR$gene_ID = substr(gff.5UTR$gene_ID,1,9)
	#gff.3UTR$UTR_ID = gff.3UTR$gene_ID
	#gff.3UTR$gene_ID = substr(gff.3UTR$gene_ID,1,9)
	
	# Add the exon and intron numbers
	#gff.exons$exon_no = as.numeric(str_split_fixed(gff.exons$exon_ID, ':',3)[,3])
	#gff.introns$exon_no = as.numeric(str_split_fixed(gff.introns$intron_ID, ':',3)[,3])
	
	# Find out how many exons each gene has
	#install.packages("sqldf")
	#library(sqldf)
	#gff.genes=merge(gff.genes, sqldf('SELECT genes.gene_id, MAX(exon_no) AS no_exons FROM [gff.exons] exons, [gff.genes] genes WHERE genes.gene_ID=exons.gene_ID GROUP BY genes.gene_id'), by="gene_ID", all=TRUE)
	 
	# Find how many variable sites each gene has
	#gene_info = gff.genes
	#varloc_genes = NULL
	
	# Make a GRanges for gene space
	gene_ranges=makeGRangesFromDataFrame(df = gff.genes, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	# Make a GRanges for transposon space
	transposon_ranges=makeGRangesFromDataFrame(df = gff.transposons, start.field = "V4", end.field = "V5",seqnames.field = "V1")

	
	
	
	# Load the BSGenome version of the reference 
	library(BSgenome)

	# List the available genome assemblies
	#available.genomes()

	# TAIR9 and TAIR10 correspond to the same genome assembly so there is no need for a BSgenome pkg for TAIR10 :-)
	# https://stat.ethz.ch/pipermail/bioconductor/2010-December/036938.html
	library(BSgenome.Athaliana.TAIR.TAIR9)

	# Find the sequence lengths for each chromosome
	sLengths=seqlengths(Athaliana)


	# Use the MethylSeekR package to segment the methylome
	library(MethylSeekR)

	# MethylSeekR expects an input file with 4 columns:  Chrom, Locus, T and M.  T is total reads, M is methylated reads.
	
	# Create an input file for MethylSeekR based on the merged methylomes of all valid samples
	
	# Sum the M (Cov_C) and T (Cov_C+Cov_T) counts for all valid samples 

	# Create empty tables to accumulate each of the 'across samples' columns data
	#across_samples_TM=matrix(0, nrow=nrow(all_samples_meth_status), ncol=2)

	# We are only dealing with one sample at a time so:
	#valid_samples = sample_id
	
	#for(sample_name in valid_samples) {
	#	if(opt$verbose) {cat(paste0("Adding methylation read counts to across-samples per-site table: ",sample_name,"\n"))}
		# Merge each of the sample counts into the relevant accumulation table, accounting for possible NA values

	#	across_samples_TM[,1]=ifelse(is.na(coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T),across_samples_TM[,1],across_samples_TM[,1]+coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T)

	#	across_samples_TM[,2]=ifelse(is.na(coverage_data[coverage_data$Sample==sample_name,]$Cov_C),across_samples_TM[,2],across_samples_TM[,2]+coverage_data[coverage_data$Sample==sample_name,]$Cov_C)

	#} # end for each sample

	
	### Read in the SNPs data
	#snps.gr <- readSNPTable(FileName=snpFname, seqLengths=sLengths)
    # We will set sample_C and sample_T to 0 for all C loci where there is an annotated mutation
	


	meth_by_segment <- function (m, segment_model, meth_context, num.cores = 1, myGenomeSeq, seqLengths, nSite.smoothing = 3, mMeth.classification=c(0.7), mMeth.classes=c("LMR","FMR")) 
	{
		# Replaced nCG with nSites for generality
		
		# m is a genomic ranges object representing a methylome
		# segment_model is an arbitrary genomic ranges object representing a set of query segments of interest

		# Count CG sites in each segment
		# Added fixed=FALSE so that ambiguity codes will work properly
		nSites = vcountPattern(meth_context, getSeq(Athaliana, resize(segment_model, width(segment_model), fix = "start"), as.character = FALSE), fixed=FALSE)
	
		# Overlap the methylome with the segments
		ov <- findOverlaps(m, segment_model)

		# Count Ts and Ms per segment
		T = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), sum)
		M = tapply(values(m[queryHits(ov)])[, 2], subjectHits(ov), sum)
		nSites.segmentation = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), length)
		# The next step carries out a robust smoothing of the methylation values across the segment.  It takes k adjacent sites (default=3) and calculates M/T for a sliding window of that size, then takes the median value of this proportion across the segment. k=3 give 'minimal' robust smoothing eliminating isolated outliers.
		# This might be problematic for us, as in cases where each site is either fully methylated or unmethylated, it will tend to generate sequences like this for a window size of 3:
		# 0,0,0,1/3,1/3,1/3,1/3,2/3,2/3,2/3,1,1,1,1,1  - this will lead to clusters of median methylation at 1/3 or 2/3 except for mostly unmethylated or methylated segments.
		median.meth = tapply(as.vector(runmean(Rle(values(m[queryHits(ov)])[,2]/values(m[queryHits(ov)])[, 1]), nSite.smoothing, endrule = "constant")), subjectHits(ov), median)
		median.meth = pmax(0, median.meth)
		#if (!all.equal(as.numeric(names(T)), 1:length(segment_model))) {
		#	message("error in calculating methylation levels for PMDs")
		#}	
		# Original code used median.meth to make a 'type' determination for the segment.  We will use pmeth instead
		#type=ifelse(median.meth<mMeth.classification[1],mMeth.classes[1], mMeth.classes[2])
		type=ifelse(M/T<mMeth.classification[1],mMeth.classes[1], mMeth.classes[2])
		#type=ifelse(median.meth<mMeth.classification[1],0,ifelse(median.meth<mMeth.classification[2],1,ifelse(median.meth<mMeth.classification[3],2,3)))
		#type = c("UMR", "LMR")[as.integer(nCG < nCG.classification) + 1]

		DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
	}




doing_non_cg = TRUE
m.sel=0.09
not_processed_samples = c()
processed_samples = c()
all_samples_meth_status = NA
no_samples = 0

# make empty structures for the various things we are accumulating
gene_gbm_overlap = NULL
gene_tem_overlap = NULL
gene_umr_overlap = NULL
gene_gbm_mCG_sites = NULL
gene_tem_mCG_sites = NULL

for (this_sample in sample_list$SRA.Accession) {
  cat(paste0(this_sample,"\n"))
  
  if (doing_non_cg) {
    # read in CHG methylome
    meth_CHG.gr <- readMethylome(FileName=paste0(this_sample,"/",project_id,"_CHG_",this_sample,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)

    # read in CHH methylome
    meth_CHH.gr <- readMethylome(FileName=paste0(this_sample,"/",project_id,"_CHH_",this_sample,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)

    # Set aside genes with no CHG sites
    set_aside_no_CHG.gr = setdiff(gene_ranges, gene_ranges[unique(subjectHits(findOverlaps(meth_CHG.gr, gene_ranges)))])
    #cat(paste0("set_aside_no_CHG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
    # 335 in first sample

    # Set aside genes with no CHH sites
    set_aside_no_CHH.gr = setdiff(gene_ranges, gene_ranges[unique(subjectHits(findOverlaps(meth_CHH.gr, gene_ranges)))])
    #cat(paste0("set_aside_no_CHH.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
    # 267 in first sample

    CHG_genes = unique(subjectHits(findOverlaps(meth_CHG.gr, gene_ranges)))
    gene_ranges_CHG = gene_ranges[CHG_genes]
    #cat(paste0("gene_ranges_CHG contains ",length(gene_ranges_CHG)," segments after trimming segments with no CHG sites.\n"))

    CHH_genes = unique(subjectHits(findOverlaps(meth_CHH.gr, gene_ranges)))
    gene_ranges_CHH = gene_ranges[CHH_genes]
    #cat(paste0("gene_ranges_CHH contains ",length(gene_ranges_CHH)," segments after trimming segments with no CHH sites.\n"))
  }
  
  # read in segmentation
  segmentation_model = read.table(file = paste0(this_sample,"/",project_id,"_","CG","_",this_sample,"_segmentation_model_draft3.tsv"), sep="\t", header=TRUE)
  
  segmentation_model$Chromosome = paste0("Chr", segmentation_model$Chromosome)
  
  # Capitalise chromosome names in segmentation_model.gr
  #levels(segmentation_model.gr@seqnames)=toupper(levels(segmentation_model.gr@seqnames))
  #segmentation_model.gr@seqinfo@seqnames=levels(segmentation_model.gr@seqnames)

  #levels(variant_call_ranges@seqnames)=toupper(levels(variant_call_ranges@seqnames))
  #variant_call_ranges@seqinfo@seqnames=levels(variant_call_ranges@seqnames)

  GBM_segments.gr = makeGRangesFromDataFrame(df = segmentation_model[substr(segmentation_model$Type,1,3)=="GBM",], start.field = "Start", end.field = "End",seqnames.field = "Chromosome")
  TEM_segments.gr = makeGRangesFromDataFrame(df = segmentation_model[substr(segmentation_model$Type,1,3)=="TEM",], start.field = "Start", end.field = "End",seqnames.field = "Chromosome")
  UMR_segments.gr = makeGRangesFromDataFrame(df = segmentation_model[substr(segmentation_model$Type,1,3)=="UMR",], start.field = "Start", end.field = "End",seqnames.field = "Chromosome")
  length(GBM_segments.gr)
  length(TEM_segments.gr)
  length(UMR_segments.gr)

  # Identify GenomicRegions covering genes 

  if (doing_non_cg) {  
    # overlap CHG methylome with gene loci and assign methylation % to each gene
    meth_CHG = meth_by_segment(meth_CHG.gr, segment_model=gene_ranges_CHG, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))$pmeth
 
    meth_CHG = cbind(gff.genes[CHG_genes,"gene_ID"], meth_CHG)

    colnames(meth_CHG) = c("gene_ID",this_sample)

    # overlap CHH methylome with gene loci and assign methylation % to each gene
    meth_CHH = meth_by_segment(meth_CHH.gr, segment_model=gene_ranges_CHH, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))$pmeth
  
    meth_CHH = cbind(gff.genes[CHH_genes,"gene_ID"], meth_CHH)

    colnames(meth_CHH) = c("gene_ID",this_sample)
  }

  # overlap gbM segments with genes and assign overlap statistic to each gene
  gbm_genes = gff.genes[unique(subjectHits(findOverlaps(GBM_segments.gr, gene_ranges))),"gene_ID"]
  gbm_genes = cbind(gbm_genes, rep(1, length(gbm_genes)))
  colnames(gbm_genes) = c("gene_ID", this_sample)
  
  # overlap UMR segments with genes and assign overlap statistic to each gene
  umr_genes = gff.genes[unique(subjectHits(findOverlaps(UMR_segments.gr, gene_ranges))),"gene_ID"]
  umr_genes = cbind(umr_genes, rep(1, length(umr_genes)))
  colnames(umr_genes) = c("gene_ID", this_sample)
  
  # overlap TEM segments with genes and assign overlap statistic to each gene
  tem_genes = gff.genes[unique(subjectHits(findOverlaps(TEM_segments.gr, gene_ranges))),"gene_ID"]
  tem_genes = cbind(tem_genes, rep(1, length(tem_genes)))
  colnames(tem_genes) = c("gene_ID", this_sample)
  
  # overlap TEM segments with genes and assign overlap statistic to each gene - length of overlapping segment
  
  # Overlap genes with all_samples_meth_status and count No. of mCGs per gene, No. of sites per gene, prop of mCGs per gene
  
  # overlap TEM segments and gbM segments with CG sites and count No. mCGs in each segment in each gene
  # identify mCG sites overlapping gbM segments and TEM segments
  this_meth = all_samples_meth_status[all_samples_meth_status[,colnames(all_samples_meth_status)==this_sample]=="M",c("Chromosome", "Locus")]
  this_meth=cbind(this_meth, this_meth$Locus+1)
  colnames(this_meth)=c("chromosome","start","end")
  this_meth.gr = makeGRangesFromDataFrame(df = this_meth[!is.na(this_meth$start),])
  
  GBM_mCG_sites = this_meth.gr[queryHits(findOverlaps(this_meth.gr, GBM_segments.gr))]
  TEM_mCG_sites = this_meth.gr[queryHits(findOverlaps(this_meth.gr, TEM_segments.gr))]
  this_gene_GBM_mCG_sites = table(gff.genes[subjectHits(findOverlaps(GBM_mCG_sites, gene_ranges)),"gene_ID"])
  this_gene_TEM_mCG_sites = table(gff.genes[subjectHits(findOverlaps(TEM_mCG_sites, gene_ranges)),"gene_ID"])
  this_gene_GBM_mCG_sites = as.data.frame(this_gene_GBM_mCG_sites)
  colnames(this_gene_GBM_mCG_sites) = c("gene_ID", this_sample)
  this_gene_TEM_mCG_sites = as.data.frame(this_gene_TEM_mCG_sites)
  colnames(this_gene_TEM_mCG_sites) = c("gene_ID", this_sample)
  
  # add each of the statistics vectors to the relevant all samples matrix
  if (no_samples == 0) {
    if (doing_non_cg) {
      gene_meth_CHG = as.data.frame(meth_CHG)
      gene_meth_CHH = as.data.frame(meth_CHH)
  	}
    gene_gbm_overlap = as.data.frame(gbm_genes)
    gene_tem_overlap = as.data.frame(tem_genes)
    gene_umr_overlap = as.data.frame(umr_genes)
    gene_gbm_mCG_sites = this_gene_GBM_mCG_sites
    gene_tem_mCG_sites = this_gene_TEM_mCG_sites
    
  } else {
    if (doing_non_cg) {
      gene_meth_CHG = merge(gene_meth_CHG, meth_CHG, by="gene_ID", all=TRUE)
      gene_meth_CHH = merge(gene_meth_CHH, meth_CHH, by="gene_ID", all=TRUE)
	  }
    gene_gbm_overlap = merge(gene_gbm_overlap, gbm_genes, by="gene_ID", all=TRUE)
    gene_tem_overlap = merge(gene_tem_overlap, tem_genes, by="gene_ID", all=TRUE)
    gene_umr_overlap = merge(gene_umr_overlap, umr_genes, by="gene_ID", all=TRUE)
    gene_gbm_mCG_sites = merge(gene_gbm_mCG_sites, this_gene_GBM_mCG_sites, by="gene_ID", all=TRUE)
    gene_tem_mCG_sites = merge(gene_tem_mCG_sites, this_gene_TEM_mCG_sites, by="gene_ID", all=TRUE)
  }   
  
  # These lines not needed for R 4.0 and higher
  #if (doing_non_cg) {
    #gene_meth_CHG[,no_samples+2] = as.numeric(levels(gene_meth_CHG[,no_samples+2])[gene_meth_CHG[,no_samples+2]])
    #gene_meth_CHH[,no_samples+2] = as.numeric(levels(gene_meth_CHH[,no_samples+2])[gene_meth_CHH[,no_samples+2]])
  #}
  #gene_gbm_overlap[,no_samples+2] = as.numeric(levels(gene_gbm_overlap[,no_samples+2])[gene_gbm_overlap[,no_samples+2]])
  #gene_tem_overlap[,no_samples+2] = as.numeric(levels(gene_tem_overlap[,no_samples+2])[gene_tem_overlap[,no_samples+2]])
  #gene_umr_overlap[,no_samples+2] = as.numeric(levels(gene_umr_overlap[,no_samples+2])[gene_umr_overlap[,no_samples+2]])
  #gene_gbm_mCG_sites[,no_samples+2] = as.numeric(levels(gene_gbm_mCG_sites[,no_samples+2])[gene_gbm_mCG_sites[,no_samples+2]])
  #gene_tem_mCG_sites[,no_samples+2] = as.numeric(levels(gene_tem_mCG_sites[,no_samples+2])[gene_tem_mCG_sites[,no_samples+2]])

  # These lines are the version for R 4.0 and higher
  if (doing_non_cg) {
    gene_meth_CHG[,no_samples+2] = as.numeric(gene_meth_CHG[,no_samples+2])
    gene_meth_CHH[,no_samples+2] = as.numeric(gene_meth_CHH[,no_samples+2])
  }
  gene_gbm_overlap[,no_samples+2] = as.numeric(gene_gbm_overlap[,no_samples+2])
  gene_tem_overlap[,no_samples+2] = as.numeric(gene_tem_overlap[,no_samples+2])
  gene_umr_overlap[,no_samples+2] = as.numeric(gene_umr_overlap[,no_samples+2])
  #gene_gbm_mCG_sites[,no_samples+2] = as.numeric(levels(gene_gbm_mCG_sites[,no_samples+2])[gene_gbm_mCG_sites[,no_samples+2]])
  #gene_tem_mCG_sites[,no_samples+2] = as.numeric(levels(gene_tem_mCG_sites[,no_samples+2])[gene_tem_mCG_sites[,no_samples+2]])
  
  no_samples = no_samples + 1
  processed_samples=c(processed_samples, this_sample)
#  this_sample_meth_status = data.frame(readRDS(this_file))
  if (no_samples == 1) {
#    all_samples_meth_status = this_sample_meth_status
  } else {
#	all_samples_meth_status = merge(all_samples_meth_status, this_sample_meth_status, by=c("Chromosome","Locus","Strand"), all=TRUE)
  }
#} else {
#  not_processed_samples=c(not_processed_samples, this_sample)
#}
}

# Replace NAs with 0
gene_meth_CHG[is.na(gene_meth_CHG)] = 0
gene_meth_CHH[is.na(gene_meth_CHH)] = 0
gene_gbm_overlap[is.na(gene_gbm_overlap)] = 0
gene_tem_overlap[is.na(gene_tem_overlap)] = 0
gene_umr_overlap[is.na(gene_umr_overlap)] = 0
gene_gbm_mCG_sites[is.na(gene_gbm_mCG_sites)] = 0
gene_tem_mCG_sites[is.na(gene_tem_mCG_sites)] = 0

write.table(gene_meth_CHG, file=paste0(project_id, "_gene_meth_CHG.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_meth_CHH, file=paste0(project_id, "_gene_meth_CHH.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_gbm_overlap, file=paste0(project_id, "_gene_gbm_overlap.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_overlap, file=paste0(project_id, "_gene_tem_overlap.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_umr_overlap, file=paste0(project_id, "_gene_umr_overlap.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_gbm_mCG_sites, file=paste0(project_id, "_gene_gbm_mCG_sites.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_mCG_sites, file=paste0(project_id, "_gene_tem_mCG_sites.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

gene_gbm_mCG_sites$gene_ID = as.character(gene_gbm_mCG_sites$gene_ID)
gene_tem_mCG_sites$gene_ID = as.character(gene_tem_mCG_sites$gene_ID)

# Dump the data table to file for a convenience cache
saveRDS(all_samples_meth_status, file=paste0(project_id,"_",meth_context,"_all_samples_meth_status_0.005_0.05_10_25.rds"))
# saved for batches 1-29 so far
#all_samples_meth_status = readRDS(file=paste0(project_id,"_",meth_context,"_all_samples_meth_status_0.005_0.05_10_25.rds"))


# Restore the data from the convenience dumps
#gene_meth_CHG = read.table(file=paste0(project_id, "_gene_meth_CHG.txt"), sep="\t", header=TRUE)
#gene_meth_CHH = read.table(file=paste0(project_id, "_gene_meth_CHH.txt"), sep="\t", header=TRUE)
#gene_gbm_overlap = read.table(file=paste0(project_id, "_gene_gbm_overlap.txt"), sep="\t", header=TRUE)
#gene_tem_overlap = read.table(file=paste0(project_id, "_gene_tem_overlap.txt"), sep="\t", header=TRUE)
#gene_umr_overlap = read.table(file=paste0(project_id, "_gene_umr_overlap.txt"), sep="\t", header=TRUE)
#gene_gbm_mCG_sites = read.table(file=paste0(project_id, "_gene_gbm_mCG_sites.txt"), sep="\t", header=TRUE)
#gene_tem_mCG_sites = read.table(file=paste0(project_id, "_gene_tem_mCG_sites.txt"), sep="\t", header=TRUE)
#all_samples_meth_status = readRDS(file=paste0("../5-analysis_2020-01-17/",project_id,"_",meth_context,"_all_samples_meth_status_0.005_0.05_10_25.rds"))


# Read the SNP array originating from 1001 genomes website https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf
# parse_1001genome_VCF.pl/.sh together should have generated a separate file of mutation sites in C and CG contexts for each sample ID, so read them in one sample at a time.
# 1001 genomes SNP files are named for the accession No. in 1001 genome project. We are assuming they correspond directly with the Ecotype_ID we pull from the substring of the sample Title string
samples_without_genomes = c()
samples_with_genomes = c()
no_ecotypes = 0
for (this_sample in processed_samples) {
  # Find the ecotype corresponding to this sample
  this_ecotype = sample_list[sample_list$SRA.Accession==this_sample, "Ecotype_id"]
  this_file = paste0("../0-reference/",this_ecotype,"_CG_vars.txt")
  if (file.exists(this_file)) {
	if(opt$verbose) {cat(paste0("Loading ",this_ecotype," SNPs\n"))}
    no_ecotypes = no_ecotypes + 1
    samples_with_genomes=c(samples_with_genomes, this_sample)
	# read in the known CG sites with SNPs for this sample.  Use unique to deduplicate cases where the same site is listed more than once
	this_sample_snps = unique(read.table(this_file, header=FALSE, sep="\t"))
	colnames(this_sample_snps) = c("Chromosome", "Locus")
	#cat(nrow(this_sample_snps[(this_sample_snps$Chromosome=="Chr5") & (this_sample_snps$Locus>=20048325) & (this_sample_snps$Locus<=20049775),]))
	#merge snps and sample. merge the merged set with sample to get rid of any rows exclusively from 
	mutant_sites = merge(all_samples_meth_status[,c("Chromosome", "Locus", this_sample)], cbind(this_sample_snps, "mutation"=rep(1,nrow(this_sample_snps))), by = c("Chromosome", "Locus"), all.x=TRUE, all.y=FALSE) 
    all_samples_meth_status[, this_sample] = ifelse(is.na(mutant_sites$mutation), as.character(all_samples_meth_status[, this_sample]), "X")
  } else {
    samples_without_genomes=c(samples_without_genomes, this_sample)
  }
}

# Dump the data table to file for a convenience cache
saveRDS(all_samples_meth_status, file=paste0(project_id,"_",meth_context,"_all_samples_meth_status_0.005_0.05_10_25_SNPs.rds"))
# saved for batches 1-29 so far
#all_samples_meth_status = readRDS(file=paste0(project_id,"_",meth_context,"_all_samples_meth_status_0.005_0.05_10_25_SNPs.rds"))

# work out how many genes are both gbM and TEM
gene_both_overlap = NULL
gene_both_overlap_filtered = NULL
gene_both_overlap_resolved = NULL
gene_both_overlap_resolved_gbm = NULL
gene_both_overlap_resolved_tem = NULL
gene_gbm_calls = NULL
gene_tem_calls = NULL
no_samples = 0
for (this_sample in sample_list$SRA.Accession) {
  cat(paste0(this_sample,"\n"))
  gene_stats = merge(cbind("gene_ID"=gene_gbm_overlap[,"gene_ID"], "gene_gbm_overlap"=gene_gbm_overlap[,colnames(gene_gbm_overlap)==this_sample]), cbind("gene_ID"=gene_tem_overlap[,"gene_ID"], "gene_tem_overlap"=gene_tem_overlap[,colnames(gene_tem_overlap)==this_sample]), by="gene_ID", all=TRUE)
  gene_stats = merge(gene_stats, cbind("gene_ID"=gene_gbm_mCG_sites[,"gene_ID"], "gene_gbm_mCG_sites"=gene_gbm_mCG_sites[,colnames(gene_gbm_mCG_sites)==this_sample]), by="gene_ID", all=TRUE)
  gene_stats = merge(gene_stats, cbind("gene_ID"=gene_tem_mCG_sites[,"gene_ID"], "gene_tem_mCG_sites"=gene_tem_mCG_sites[,colnames(gene_tem_mCG_sites)==this_sample]), by="gene_ID", all=TRUE)
  gene_stats[,2] = as.numeric(gene_stats[,2])
  gene_stats[,3] = as.numeric(gene_stats[,3])
  gene_stats[,4] = as.numeric(gene_stats[,4])
  gene_stats[,5] = as.numeric(gene_stats[,5])
  
  gene_stats$both_overlap = gene_stats$gene_gbm_overlap * gene_stats$gene_tem_overlap

  gene_stats$filtered = gene_stats$both_overlap*ifelse(gene_stats$gene_gbm_mCG_sites>2,1,0)*ifelse(gene_stats$gene_tem_mCG_sites>2,1,0)
  gbm_tem_ratio = 4
  gene_stats$resolved = gene_stats$both_overlap*(ifelse(gene_stats$gene_gbm_mCG_sites>gene_stats$gene_tem_mCG_sites*gbm_tem_ratio,1,0) + ifelse(gene_stats$gene_tem_mCG_sites>gene_stats$gene_gbm_mCG_sites*gbm_tem_ratio,1,0))*(ifelse((gene_stats$gene_gbm_mCG_sites>2) | (gene_stats$gene_tem_mCG_sites>2),1,0))*(ifelse(!((gene_stats$gene_gbm_mCG_sites>2) & (gene_stats$gene_tem_mCG_sites>2)),1,0))
  gene_stats$resolved_gbm = gene_stats$both_overlap*(ifelse(gene_stats$gene_gbm_mCG_sites>gene_stats$gene_tem_mCG_sites*gbm_tem_ratio,1,0))*(ifelse((gene_stats$gene_gbm_mCG_sites>2) | (gene_stats$gene_tem_mCG_sites>2),1,0))*(ifelse(!((gene_stats$gene_gbm_mCG_sites>2) & (gene_stats$gene_tem_mCG_sites>2)),1,0))
  gene_stats$resolved_tem = gene_stats$both_overlap*(ifelse(gene_stats$gene_tem_mCG_sites>gene_stats$gene_gbm_mCG_sites*gbm_tem_ratio,1,0))*(ifelse((gene_stats$gene_gbm_mCG_sites>2) | (gene_stats$gene_tem_mCG_sites>2),1,0))*(ifelse(!((gene_stats$gene_gbm_mCG_sites>2) & (gene_stats$gene_tem_mCG_sites>2)),1,0))
  # replace NAs with 0s, so that they don't proliferate when we start comparing columns
  gene_stats[is.na(gene_stats)] = 0
  # collect the calls that are unambiguous, or meet the filtering criteria
  gene_stats$gbm_calls = gene_stats$gene_gbm_overlap * ((1-gene_stats$both_overlap) + (gene_stats$resolved_gbm))
  gene_stats$tem_calls = gene_stats$gene_tem_overlap * ((1-gene_stats$both_overlap) + (gene_stats$resolved_tem))
  
  if (no_samples == 0) {
    gene_both_overlap = as.data.frame(gene_stats[,c(1,6)])
    gene_both_overlap_filtered = as.data.frame(gene_stats[,c(1,7)])
    gene_both_overlap_resolved = as.data.frame(gene_stats[,c(1,8)])
    gene_both_overlap_resolved_gbm = as.data.frame(gene_stats[,c(1,9)])
    gene_both_overlap_resolved_tem = as.data.frame(gene_stats[,c(1,10)])
    gene_gbm_calls = as.data.frame(gene_stats[,c(1,11)])
    gene_tem_calls = as.data.frame(gene_stats[,c(1,12)])
  } else {
    gene_both_overlap = merge(gene_both_overlap, gene_stats[,c(1,6)], by="gene_ID", all=TRUE)
    gene_both_overlap_filtered = merge(gene_both_overlap_filtered, gene_stats[,c(1,7)], by="gene_ID", all=TRUE)
    gene_both_overlap_resolved = merge(gene_both_overlap_resolved, gene_stats[,c(1,8)], by="gene_ID", all=TRUE)
    gene_both_overlap_resolved_gbm = merge(gene_both_overlap_resolved_gbm, gene_stats[,c(1,9)], by="gene_ID", all=TRUE)
    gene_both_overlap_resolved_tem = merge(gene_both_overlap_resolved_tem, gene_stats[,c(1,10)], by="gene_ID", all=TRUE)
    gene_gbm_calls = merge(gene_gbm_calls, gene_stats[,c(1,11)], by="gene_ID", all=TRUE)
    gene_tem_calls = merge(gene_tem_calls, gene_stats[,c(1,12)], by="gene_ID", all=TRUE)
  }   
  colnames(gene_both_overlap)[no_samples+2] = this_sample
  colnames(gene_both_overlap_filtered)[no_samples+2] = this_sample
  colnames(gene_both_overlap_resolved)[no_samples+2] = this_sample
  colnames(gene_both_overlap_resolved_gbm)[no_samples+2] = this_sample
  colnames(gene_both_overlap_resolved_tem)[no_samples+2] = this_sample
  colnames(gene_gbm_calls)[no_samples+2] = this_sample
  colnames(gene_tem_calls)[no_samples+2] = this_sample
  
  
  # These lines not needed for R 4.0 and higher
  #if (doing_non_cg) {
  #gene_both_overlap[,no_samples+2] = as.numeric(levels(gene_both_overlap[,no_samples+2])[gene_both_overlap[,no_samples+2]])

  # These lines are the version for R 4.0 and higher
  gene_both_overlap[,no_samples+2] = as.numeric(gene_both_overlap[,no_samples+2])
  gene_both_overlap_filtered[,no_samples+2] = as.numeric(gene_both_overlap_filtered[,no_samples+2])
  gene_both_overlap_resolved[,no_samples+2] = as.numeric(gene_both_overlap_resolved[,no_samples+2])
  gene_both_overlap_resolved_gbm[,no_samples+2] = as.numeric(gene_both_overlap_resolved_gbm[,no_samples+2])
  gene_both_overlap_resolved_tem[,no_samples+2] = as.numeric(gene_both_overlap_resolved_tem[,no_samples+2])
  gene_gbm_calls[,no_samples+2] = as.numeric(gene_gbm_calls[,no_samples+2])
  gene_tem_calls[,no_samples+2] = as.numeric(gene_tem_calls[,no_samples+2])
  
    no_samples = no_samples + 1
}


# plot distributions of numbers of gbM genes, TEM genes and both per accession
gbm_gene_counts = colSums(gene_gbm_overlap[,2:1212], na.rm=TRUE)
tem_gene_counts = colSums(gene_tem_overlap[,2:1212], na.rm=TRUE)
both_counts = colSums(gene_both_overlap[2:1212], na.rm=TRUE)
both_filtered_counts = colSums(gene_both_overlap_filtered[2:1212], na.rm=TRUE)
both_resolved_counts = colSums(gene_both_overlap_resolved[2:1212], na.rm=TRUE)
both_resolved_gbm_counts = colSums(gene_both_overlap_resolved_gbm[2:1212], na.rm=TRUE)
both_resolved_tem_counts = colSums(gene_both_overlap_resolved_tem[2:1212], na.rm=TRUE)

pdf(file=paste0(project_id, "_gene_gbm_tem_stats.pdf"))
print(ggplot(as.data.frame(gbm_gene_counts)) + geom_histogram(aes(x=gbm_gene_counts), binwidth = 100))
print(ggplot(as.data.frame(tem_gene_counts)) + geom_histogram(aes(x=tem_gene_counts), binwidth = 100))
print(ggplot(as.data.frame(cbind(gbm_gene_counts, tem_gene_counts))) + geom_point(aes(x=gbm_gene_counts, y=tem_gene_counts)))
print(ggplot(as.data.frame(both_counts)) + geom_histogram(aes(x=both_counts), binwidth = 10))
print(ggplot(as.data.frame(cbind(gbm_gene_counts, both_counts))) + geom_point(aes(x=gbm_gene_counts, y=both_counts)))
print(ggplot(as.data.frame(cbind(tem_gene_counts, both_counts))) + geom_point(aes(x=tem_gene_counts, y=both_counts)))
print(ggplot(as.data.frame(both_filtered_counts)) + geom_histogram(aes(x=both_filtered_counts), binwidth = 10))
print(ggplot(as.data.frame(both_resolved_counts)) + geom_histogram(aes(x=both_resolved_counts), binwidth = 10))
print(ggplot(as.data.frame(cbind(both_counts, both_filtered_counts))) + geom_point(aes(x=both_counts, y=both_filtered_counts)))
print(ggplot(as.data.frame(cbind(both_counts, both_resolved_counts))) + geom_point(aes(x=both_counts, y=both_resolved_counts)))
print(ggplot(as.data.frame(cbind(both_resolved_gbm_counts, both_resolved_tem_counts))) + geom_point(aes(x=both_resolved_gbm_counts, y=both_resolved_tem_counts)))


dev.off()

cor(cbind(gbm_gene_counts, both_counts), use="complete.obs")
#gbm_gene_counts both_counts
#gbm_gene_counts       1.0000000   0.2897491
#both_counts           0.2897491   1.0000000
cor(cbind(tem_gene_counts, both_counts), use="complete.obs")
#tem_gene_counts both_counts
#tem_gene_counts       1.0000000   0.7091566
#both_counts           0.7091566   1.0000000
cor(cbind(tem_gene_counts, gbm_gene_counts), use="complete.obs")
#tem_gene_counts gbm_gene_counts
#tem_gene_counts       1.0000000      -0.1927565
#gbm_gene_counts      -0.1927565       1.0000000
cor(cbind(both_counts, both_filtered_counts), use="complete.obs")
#both_counts both_filtered_counts
#both_counts             1.000000             0.451117
#both_filtered_counts    0.451117             1.000000
cor(cbind(both_counts, both_resolved_counts), use="complete.obs")
#both_counts both_resolved_counts
#both_counts            1.0000000            0.9018878
#both_resolved_counts   0.9018878            1.0000000



# before we export this for use in GWAS analysis, let's check out whether some of these calls are shite due to low coverage, and if so, replace them with NA

# Also, we don't currently have an active annotation of unmethylated genes - we just assume it's the ones which are not annotated as methylated.  Once we introduce the concept of n data, we will need a positive state of unmethylated

# First check that the dimensions of the things we've already made are equal
sum(colnames(gene_gbm_calls)==colnames(gene_tem_calls))
#[1] 1212
ncol(gene_gbm_calls)
#[1] 1212
sum(colnames(gene_gbm_calls)!=colnames(gene_tem_calls))
#[1] 0
sum(gene_gbm_calls$gene_ID==gene_tem_calls$gene_ID)
#[1] 32925
nrow(gene_gbm_calls)
#[1] 32925
sum(gene_gbm_calls$gene_ID==gene_both_overlap_filtered$gene_ID)
#[1] 32925
sum(colnames(gene_gbm_calls)==colnames(gene_both_overlap_filtered))
#[1] 1212


# We make a single matrix of all genes in all samples. We assign all to U, then we add in the G, T and Bs from above.
accession_gene_calls = as.data.frame(matrix("U", nrow=nrow(gene_gbm_calls), ncol=ncol(gene_gbm_calls)))
colnames(accession_gene_calls) = colnames(gene_gbm_calls)
rownames(accession_gene_calls) = rownames(gene_gbm_calls)
accession_gene_calls$gene_ID = gene_gbm_calls$gene_ID

accession_gene_calls[gene_gbm_calls==1] = "G"
accession_gene_calls[gene_tem_calls==1] = "T"
accession_gene_calls[gene_both_overlap_filtered==1] = "B"

#accession_gene_calls still has ~400 missing genes cf. gff.genes. These are (mostly) the C and M genes. Add them in, to get gene loci, then get rid of them again:
accession_gene_calls = merge(accession_gene_calls, gff.genes, by="gene_ID", all=TRUE)
accession_gene_calls = accession_gene_calls[!(substr(accession_gene_calls$gene_ID,3,3) %in% c("C","M")),]  

relevant_gene_ranges = makeGRangesFromDataFrame(df = accession_gene_calls, start.field = "V4", end.field = "V5",seqnames.field = "V1")
levels(relevant_gene_ranges@seqnames@values) = substr(as.character(levels(relevant_gene_ranges@seqnames@values)),4,4)
relevant_gene_ranges@seqinfo@seqnames = levels(relevant_gene_ranges@seqnames@values)

# load in all the alternative, trimmed segmentation models, discard any segments with less than 3 sites, and for each gene in each sample, if it is gbm or both, check that it has a gbM segment, else set it to NA
model_dir = "segmentation_models_v2"

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=3

ignore_list=c()
accession_gene_calls_backup = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    # These are the gbM segments, trimmed to exclude empty segments, with more than 3 CG sites
    # Now find any genes marked as gbm or both which do not overlap with one of these segments
    gene_model_overlaps = findOverlaps(relevant_gene_ranges, this_model_trimmed)
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were G as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="G"] = "X"
    # for the ones which overlap a gene and were originally G, set them back to G
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[unique(queryHits(gene_model_overlaps)),2] = "G"
    sample_calls[-unique(queryHits(gene_model_overlaps)),2] = "X"
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]=="G", "G", "I"), ifelse(sample_calls[,2]=="G", ifelse(sample_calls[,1]=="B","B","I"), ifelse(sample_calls[,1]=="B","I",sample_calls[,1])))
    
    #model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    #CG_site_gbM_counts[subjectHits(model_site_olaps)] = CG_site_gbM_counts[subjectHits(model_site_olaps)] + 1
        
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# What do the Col-0 calls look like now?
View(accession_gene_calls[,c("gene_ID","SRX446038", "SRX446039", "SRX446040", "SRX248644")])

# check the U genes now. For each gene in each accession, where it is marked as U, check that it has at least 3 sites in all_samples_meth_status set to U, else set it to NA

accession_gene_calls_backup2 = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup2

min_U_sites = 3

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))

    relevant_gene_site_olaps = findOverlaps(relevant_gene_ranges, CG_site_ranges[!is.na(all_samples_meth_status[,this_sample]) & (all_samples_meth_status[,this_sample]=="U")])
    sample_gene_U_counts = table(queryHits(relevant_gene_site_olaps))
      
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were U as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="U"] = "X"
    # for the ones which overlap sufficient U sites, and were originally U, set them back to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[as.numeric(names(sample_gene_U_counts)),2] = sample_gene_U_counts
    sample_calls[-as.numeric(names(sample_gene_U_counts)),2] = 0
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]>=min_U_sites, "U", "I"), sample_calls[,1])

    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}


# Export the matrices to use in GWAS analysis
# These are the original matrices, where wee assume complete data, no genes overflowed with empty gbM segments etc.
write.table(gene_gbm_calls, file=paste0(project_id,"_gene_gbm_calls.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_calls, file=paste0(project_id,"_gene_tem_calls.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_both_overlap_filtered, file=paste0(project_id,"_gene_both_calls.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#gene_gbm_calls = read.table(file=paste0(project_id, "_gene_gbm_calls.txt"), sep="\t", header=TRUE)
#gene_tem_calls = read.table(file=paste0(project_id, "_gene_tem_calls.txt"), sep="\t", header=TRUE)
#gene_both_overlap_filtered = read.table(file=paste0(project_id, "_gene_both_calls.txt"), sep="\t", header=TRUE)

# This is the new approach where we filter out any gbM genes without 3 M sites, or U genes without 3 U sites
write.table(accession_gene_calls[,1:ncol(gene_gbm_calls)], file=paste0(project_id,"_gene_gbm_tem_umr_calls_0.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#  Now we try some different variations where we vary the gbM and U site count thresholds, and also, whether we treat very short gbM segments as missing data or U
# We consider the above method variant 0
# Variant 1 is as per variant 0 but with threshold of 2 rather than 3 for no sites per gbM segment

# load in all the alternative, trimmed segmentation models, discard any segments with less than 3 sites, and for each gene in each sample, if it is gbm or both, check that it has a gbM segment, else set it to NA
model_dir = "segmentation_models_v2"

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=2

ignore_list=c()
#accession_gene_calls_backup = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    # These are the gbM segments, trimmed to exclude empty segments, with more than 2 CG sites
    # Now find any genes marked as gbm or both which do not overlap with one of these segments
    gene_model_overlaps = findOverlaps(relevant_gene_ranges, this_model_trimmed)
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were G as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="G"] = "X"
    # for the ones which overlap a gene and were originally G, set them back to G
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[unique(queryHits(gene_model_overlaps)),2] = "G"
    sample_calls[-unique(queryHits(gene_model_overlaps)),2] = "X"
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]=="G", "G", "I"), ifelse(sample_calls[,2]=="G", ifelse(sample_calls[,1]=="B","B","I"), ifelse(sample_calls[,1]=="B","I",sample_calls[,1])))
    
    #model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    #CG_site_gbM_counts[subjectHits(model_site_olaps)] = CG_site_gbM_counts[subjectHits(model_site_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# What do the Col-0 calls look like now?
View(accession_gene_calls[,c("gene_ID","SRX446038", "SRX446039", "SRX446040", "SRX248644")])

# check the U genes now. For each gene in each accession, where it is marked as U, check that it has at least 3 sites in all_samples_meth_status set to U, else set it to NA

accession_gene_calls_backup3 = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup3

min_U_sites = 3

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    
    relevant_gene_site_olaps = findOverlaps(relevant_gene_ranges, CG_site_ranges[!is.na(all_samples_meth_status[,this_sample]) & (all_samples_meth_status[,this_sample]=="U")])
    sample_gene_U_counts = table(queryHits(relevant_gene_site_olaps))
    
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were U as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="U"] = "X"
    # for the ones which overlap sufficient U sites, and were originally U, set them back to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[as.numeric(names(sample_gene_U_counts)),2] = sample_gene_U_counts
    sample_calls[-as.numeric(names(sample_gene_U_counts)),2] = 0
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]>=min_U_sites, "U", "I"), sample_calls[,1])
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}


# Export the matrices to use in GWAS analysis
# These are the original matrices, where wee assume complete data, no genes overflowed with empty gbM segments etc.
write.table(gene_gbm_calls, file=paste0(project_id,"_gene_gbm_calls_1.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_calls, file=paste0(project_id,"_gene_tem_calls_1.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_both_overlap_filtered, file=paste0(project_id,"_gene_both_calls_1.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#gene_gbm_calls = read.table(file=paste0(project_id, "_gene_gbm_calls.txt"), sep="\t", header=TRUE)
#gene_tem_calls = read.table(file=paste0(project_id, "_gene_tem_calls.txt"), sep="\t", header=TRUE)
#gene_both_overlap_filtered = read.table(file=paste0(project_id, "_gene_both_calls.txt"), sep="\t", header=TRUE)

# This is the new approach where we filter out any gbM genes without 2 M sites, or U genes without 3 U sites
write.table(accession_gene_calls[,1:ncol(gene_gbm_calls)], file=paste0(project_id,"_gene_gbm_tem_umr_calls_1.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Variant 2 is as per variant 0 but with threshold of 1 rather than 3 for no sites per gbM segment

# load in all the alternative, trimmed segmentation models, discard any segments with less than 3 sites, and for each gene in each sample, if it is gbm or both, check that it has a gbM segment, else set it to NA
model_dir = "segmentation_models_v2"

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=1

ignore_list=c()
#accession_gene_calls_backup = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    # These are the gbM segments, trimmed to exclude empty segments, with more than 2 CG sites
    # Now find any genes marked as gbm or both which do not overlap with one of these segments
    gene_model_overlaps = findOverlaps(relevant_gene_ranges, this_model_trimmed)
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were G as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="G"] = "X"
    # for the ones which overlap a gene and were originally G, set them back to G
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[unique(queryHits(gene_model_overlaps)),2] = "G"
    sample_calls[-unique(queryHits(gene_model_overlaps)),2] = "X"
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]=="G", "G", "I"), ifelse(sample_calls[,2]=="G", ifelse(sample_calls[,1]=="B","B","I"), ifelse(sample_calls[,1]=="B","I",sample_calls[,1])))
    
    #model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    #CG_site_gbM_counts[subjectHits(model_site_olaps)] = CG_site_gbM_counts[subjectHits(model_site_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# What do the Col-0 calls look like now?
View(accession_gene_calls[,c("gene_ID","SRX446038", "SRX446039", "SRX446040", "SRX248644")])

# check the U genes now. For each gene in each accession, where it is marked as U, check that it has at least 3 sites in all_samples_meth_status set to U, else set it to NA

accession_gene_calls_backup5 = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup5

min_U_sites = 3

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    
    relevant_gene_site_olaps = findOverlaps(relevant_gene_ranges, CG_site_ranges[!is.na(all_samples_meth_status[,this_sample]) & (all_samples_meth_status[,this_sample]=="U")])
    sample_gene_U_counts = table(queryHits(relevant_gene_site_olaps))
    
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were U as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="U"] = "X"
    # for the ones which overlap sufficient U sites, and were originally U, set them back to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[as.numeric(names(sample_gene_U_counts)),2] = sample_gene_U_counts
    sample_calls[-as.numeric(names(sample_gene_U_counts)),2] = 0
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]>=min_U_sites, "U", "I"), sample_calls[,1])
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# variant 5 uses variant 2 calls, so save them for later
variant_2_calls = accession_gene_calls

# Export the matrices to use in GWAS analysis
# These are the original matrices, where wee assume complete data, no genes overflowed with empty gbM segments etc.
write.table(gene_gbm_calls, file=paste0(project_id,"_gene_gbm_calls_2.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_calls, file=paste0(project_id,"_gene_tem_calls_2.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_both_overlap_filtered, file=paste0(project_id,"_gene_both_calls_2.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#gene_gbm_calls = read.table(file=paste0(project_id, "_gene_gbm_calls.txt"), sep="\t", header=TRUE)
#gene_tem_calls = read.table(file=paste0(project_id, "_gene_tem_calls.txt"), sep="\t", header=TRUE)
#gene_both_overlap_filtered = read.table(file=paste0(project_id, "_gene_both_calls.txt"), sep="\t", header=TRUE)

# This is the new approach where we filter out any gbM genes without 1 M sites, or U genes without 3 U sites
write.table(accession_gene_calls[,1:ncol(gene_gbm_calls)], file=paste0(project_id,"_gene_gbm_tem_umr_calls_2.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# Variant 3 is as per variant 0 but instead of classifying short gbM segments as missing (I) we class them as U

# load in all the alternative, trimmed segmentation models, discard any segments with less than 3 sites, and for each gene in each sample, if it is gbm or both, check that it has a gbM segment, else set it to NA
model_dir = "segmentation_models_v2"

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=3

ignore_list=c()
#accession_gene_calls_backup = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    # These are the gbM segments, trimmed to exclude empty segments, with more than 2 CG sites
    # Now find any genes marked as gbm or both which do not overlap with one of these segments
    gene_model_overlaps = findOverlaps(relevant_gene_ranges, this_model_trimmed)
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were G as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="G"] = "X"
    # for the ones which overlap a gene and were originally G, set them back to G
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[unique(queryHits(gene_model_overlaps)),2] = "G"
    sample_calls[-unique(queryHits(gene_model_overlaps)),2] = "X"
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]=="G", "G", "I"), ifelse(sample_calls[,2]=="G", ifelse(sample_calls[,1]=="B","B","I"), ifelse(sample_calls[,1]=="B","I",sample_calls[,1])))
    
    this_model_short = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))<min_segment_sites)]
    # These are the short gbM segments
    # Now find any genes marked as G in the original annotation which overlap with one of these segments
    gene_short_overlaps = findOverlaps(relevant_gene_ranges, this_model_short)
    sample_calls = accession_gene_calls[,this_sample]
    sample_calls = cbind(sample_calls, accession_gene_calls_backup[,this_sample])
    # for the ones which overlap a short segment, were G but got changed to I, set them to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[unique(queryHits(gene_short_overlaps)),3] = "S"
    sample_calls[-unique(queryHits(gene_short_overlaps)),3] = "Z"
    accession_gene_calls[,this_sample] = ifelse((sample_calls[,1]=="I") & (sample_calls[,2]=="G") & (sample_calls[,3]=="S"), "U", sample_calls[,1])
    
    #model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    #CG_site_gbM_counts[subjectHits(model_site_olaps)] = CG_site_gbM_counts[subjectHits(model_site_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# What do the Col-0 calls look like now?
View(accession_gene_calls[,c("gene_ID","SRX446038", "SRX446039", "SRX446040", "SRX248644")])

# check the U genes now. For each gene in each accession, where it is marked as U, check that it has at least 3 sites in all_samples_meth_status set to U, else set it to NA

accession_gene_calls_backup6 = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup6

min_U_sites = 3

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    
    relevant_gene_site_olaps = findOverlaps(relevant_gene_ranges, CG_site_ranges[!is.na(all_samples_meth_status[,this_sample]) & (all_samples_meth_status[,this_sample]=="U")])
    sample_gene_U_counts = table(queryHits(relevant_gene_site_olaps))
    
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were U as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="U"] = "X"
    # for the ones which overlap sufficient U sites, and were originally U, set them back to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[as.numeric(names(sample_gene_U_counts)),2] = sample_gene_U_counts
    sample_calls[-as.numeric(names(sample_gene_U_counts)),2] = 0
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]>=min_U_sites, "U", "I"), sample_calls[,1])
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}


# Export the matrices to use in GWAS analysis
# These are the original matrices, where wee assume complete data, no genes overflowed with empty gbM segments etc.
write.table(gene_gbm_calls, file=paste0(project_id,"_gene_gbm_calls_3.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_calls, file=paste0(project_id,"_gene_tem_calls_3.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_both_overlap_filtered, file=paste0(project_id,"_gene_both_calls_3.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#gene_gbm_calls = read.table(file=paste0(project_id, "_gene_gbm_calls.txt"), sep="\t", header=TRUE)
#gene_tem_calls = read.table(file=paste0(project_id, "_gene_tem_calls.txt"), sep="\t", header=TRUE)
#gene_both_overlap_filtered = read.table(file=paste0(project_id, "_gene_both_calls.txt"), sep="\t", header=TRUE)

# This is the new approach where we report as U any gbM genes without 3 M sites, and report as missing U genes without 3 U sites
write.table(accession_gene_calls[,1:ncol(gene_gbm_calls)], file=paste0(project_id,"_gene_gbm_tem_umr_calls_3.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# Variant 4 is as per variant 0 but instead of classifying gbM segments with a single site as missing (I) we class them as U

# load in all the alternative, trimmed segmentation models, discard any segments with less than 3 sites, and for each gene in each sample, if it is gbm or both, check that it has a gbM segment, else set it to NA
model_dir = "segmentation_models_v2"

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=3

ignore_list=c()
#accession_gene_calls_backup = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    # These are the gbM segments, trimmed to exclude empty segments, with more than 2 CG sites
    # Now find any genes marked as gbm or both which do not overlap with one of these segments
    gene_model_overlaps = findOverlaps(relevant_gene_ranges, this_model_trimmed)
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were G as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="G"] = "X"
    # for the ones which overlap a gene and were originally G, set them back to G
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[unique(queryHits(gene_model_overlaps)),2] = "G"
    sample_calls[-unique(queryHits(gene_model_overlaps)),2] = "X"
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]=="G", "G", "I"), ifelse(sample_calls[,2]=="G", ifelse(sample_calls[,1]=="B","B","I"), ifelse(sample_calls[,1]=="B","I",sample_calls[,1])))
    
    this_model_short = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))<(min_segment_sites-1))]
    # These are the short gbM segments
    # Now find any genes marked as G in the original annotation which overlap with one of these segments
    gene_short_overlaps = findOverlaps(relevant_gene_ranges, this_model_short)
    sample_calls = accession_gene_calls[,this_sample]
    sample_calls = cbind(sample_calls, accession_gene_calls_backup[,this_sample])
    # for the ones which overlap a short segment, were G but got changed to I, set them to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[unique(queryHits(gene_short_overlaps)),3] = "S"
    sample_calls[-unique(queryHits(gene_short_overlaps)),3] = "Z"
    accession_gene_calls[,this_sample] = ifelse((sample_calls[,1]=="I") & (sample_calls[,2]=="G") & (sample_calls[,3]=="S"), "U", sample_calls[,1])
    
    #model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    #CG_site_gbM_counts[subjectHits(model_site_olaps)] = CG_site_gbM_counts[subjectHits(model_site_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# What do the Col-0 calls look like now?
View(accession_gene_calls[,c("gene_ID","SRX446038", "SRX446039", "SRX446040", "SRX248644")])

# check the U genes now. For each gene in each accession, where it is marked as U, check that it has at least 3 sites in all_samples_meth_status set to U, else set it to NA

accession_gene_calls_backup7 = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup7

min_U_sites = 3

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    
    relevant_gene_site_olaps = findOverlaps(relevant_gene_ranges, CG_site_ranges[!is.na(all_samples_meth_status[,this_sample]) & (all_samples_meth_status[,this_sample]=="U")])
    sample_gene_U_counts = table(queryHits(relevant_gene_site_olaps))
    
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were U as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="U"] = "X"
    # for the ones which overlap sufficient U sites, and were originally U, set them back to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[as.numeric(names(sample_gene_U_counts)),2] = sample_gene_U_counts
    sample_calls[-as.numeric(names(sample_gene_U_counts)),2] = 0
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]>=min_U_sites, "U", "I"), sample_calls[,1])
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}


# Export the matrices to use in GWAS analysis
# These are the original matrices, where wee assume complete data, no genes overflowed with empty gbM segments etc.
write.table(gene_gbm_calls, file=paste0(project_id,"_gene_gbm_calls_4.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_calls, file=paste0(project_id,"_gene_tem_calls_4.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_both_overlap_filtered, file=paste0(project_id,"_gene_both_calls_4.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#gene_gbm_calls = read.table(file=paste0(project_id, "_gene_gbm_calls.txt"), sep="\t", header=TRUE)
#gene_tem_calls = read.table(file=paste0(project_id, "_gene_tem_calls.txt"), sep="\t", header=TRUE)
#gene_both_overlap_filtered = read.table(file=paste0(project_id, "_gene_both_calls.txt"), sep="\t", header=TRUE)

# This is the new approach where we report as U any gbM genes with a single M, missing if more than 1 but less than 3 M sites, and report as missing U genes without 3 U sites
write.table(accession_gene_calls[,1:ncol(gene_gbm_calls)], file=paste0(project_id,"_gene_gbm_tem_umr_calls_4.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)



# Variant 5 As the previous permissive annotation, but applying the constraint that a gbM segment has to be entirely contained within the body of a gene for it to be used. If the gbM segment extends to one end (or both ends) of the gene we class it as missing data, unless the gene is classed as gbM in annotation 2, in which case it is classed as gbM.

# load in all the original, untrimmed segmentation models, and for each gene in each sample, if it is gbm or both, check that its gbM segment doesn't extend to the end of the gene, else check its variant 2 setting is gbM, else set it to NA
model_dir = "segmentation_models"

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=1

ignore_list=c()
#accession_gene_calls_backup = accession_gene_calls
accession_gene_calls = accession_gene_calls_backup

# variant_2_calls were saved above, but just in case not, can load them in:
#variant_2_calls = read.table(file=paste0(project_id,"_gene_gbm_tem_umr_calls_3.txt"), sep="\t", quote=FALSE, header=TRUE)

#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 10, str_length(this_model_file)-30)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model[substr(this_model$Type,1,4)=="GBM_",], seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    # These are the gbM segments, trimmed to exclude empty segments, with more than 2 CG sites
    # Now find any genes marked as gbm or both which do not overlap with one of these segments
    gene_model_overlaps = findOverlaps(relevant_gene_ranges, this_model_trimmed)
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were G as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="G"] = "X"
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    # for each of these check whether their gbM segment(s) touch the gene boundary
    gene_starts = accession_gene_calls$V4
    gene_ends = accession_gene_calls$V5
    segment_starts = findOverlaps(relevant_gene_ranges, this_model_trimmed, select="first")
    segment_ends = findOverlaps(relevant_gene_ranges, this_model_trimmed, select="last")

    #sample_calls[i,2] = ifelse(sample_calls[,1]=="X",ifelse(is.na(segment_starts), I, ifelse(this_model_trimmed[segment_starts]==gene_starts, ifelse(variant_2_calls=="G","G","I"),"G")))
    for (i in 1:nrow(sample_calls)) {
      if (is.na(sample_calls[i,1])) {
        sample_calls[i,2] = "I"
      } else {
        if (sample_calls[i,1]=="X") {
          if (is.na(segment_starts[i])) {
            sample_calls[i,2] = "I"
          } else {
            if(this_model_trimmed[segment_starts[i]]@ranges@start == gene_starts[i]) {
              if(variant_2_calls[i,this_sample]=="G") {
                sample_calls[i,2] = "G"
              } else {
                sample_calls[i,2] = "I"
              }
            } else {
              if (is.na(segment_ends[i])) {
                sample_calls[i,2] = "I"
              } else {
                if((this_model_trimmed[segment_starts[i]]@ranges@start+this_model_trimmed[segment_starts[i]]@ranges@width-1) == gene_ends[i]) {
                  if(variant_2_calls[i,this_sample]=="G") {
                    sample_calls[i,2] = "G"
                  } else {
                    sample_calls[i,2] = "I"
                  }
                } else {
                  sample_calls[i,2] = "G"
                }
              }
            }
          }
        }
      }  
    }

    # for the ones which overlap a gene and were originally G, set them back to G
    #sample_calls[-unique(queryHits(gene_model_overlaps)),2] = "X"
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", sample_calls[,2], sample_calls[,1])
    
    #model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    #CG_site_gbM_counts[subjectHits(model_site_olaps)] = CG_site_gbM_counts[subjectHits(model_site_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# What do the Col-0 calls look like now?
View(accession_gene_calls[,c("gene_ID","SRX446038", "SRX446039", "SRX446040", "SRX248644")])

# check the U genes now. For each gene in each accession, where it is marked as U, check that it has at least 3 sites in all_samples_meth_status set to U, else set it to NA
#I DID NOT BOTHER DOING THIS PART FOR VARIANT 5 - I WENT STRAIGHT TO THE EXPORT BIT BELOW
accession_gene_calls_backup8 = accession_gene_calls


min_U_sites = 3
accession_gene_calls = accession_gene_calls_backup8
#CG_site_gbM_counts = rep(0,length(CG_site_ranges))
#no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    
    relevant_gene_site_olaps = findOverlaps(relevant_gene_ranges, CG_site_ranges[!is.na(all_samples_meth_status[,this_sample]) & (all_samples_meth_status[,this_sample]=="U")])
    sample_gene_U_counts = table(queryHits(relevant_gene_site_olaps))
    
    sample_calls = accession_gene_calls[,this_sample]
    # mark all the ones which were U as X (being updated)
    sample_calls[accession_gene_calls[,this_sample]=="U"] = "X"
    # for the ones which overlap sufficient U sites, and were originally U, set them back to U
    sample_calls=cbind(sample_calls, rep(NA, length(sample_calls)))
    sample_calls[as.numeric(names(sample_gene_U_counts)),2] = sample_gene_U_counts
    sample_calls[-as.numeric(names(sample_gene_U_counts)),2] = 0
    accession_gene_calls[,this_sample] = ifelse(sample_calls[,1]=="X", ifelse(sample_calls[,2]>=min_U_sites, "U", "I"), sample_calls[,1])
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}


# Export the matrices to use in GWAS analysis
# These are the original matrices, where wee assume complete data, no genes overflowed with empty gbM segments etc.
write.table(gene_gbm_calls, file=paste0(project_id,"_gene_gbm_calls_5.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_tem_calls, file=paste0(project_id,"_gene_tem_calls_5.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(gene_both_overlap_filtered, file=paste0(project_id,"_gene_both_calls_5.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#gene_gbm_calls = read.table(file=paste0(project_id, "_gene_gbm_calls.txt"), sep="\t", header=TRUE)
#gene_tem_calls = read.table(file=paste0(project_id, "_gene_tem_calls.txt"), sep="\t", header=TRUE)
#gene_both_overlap_filtered = read.table(file=paste0(project_id, "_gene_both_calls.txt"), sep="\t", header=TRUE)

# This is the new approach where we report as U any gbM genes with a single M, missing if more than 1 but less than 3 M sites, and report as missing U genes without 3 U sites
write.table(accession_gene_calls[,1:ncol(gene_gbm_calls)], file=paste0(project_id,"_gene_gbm_tem_umr_calls_5.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)






# plot distributions of numbers of gbM genes, TEM genes and both per accession
gbm_gene_counts2 = colSums(accession_gene_calls[,2:1212]=="G", na.rm=TRUE)
tem_gene_counts2 = colSums(accession_gene_calls[,2:1212]=="T", na.rm=TRUE)
both_counts2 = colSums(accession_gene_calls[,2:1212]=="B", na.rm=TRUE)
#both_filtered_counts = colSums(gene_both_overlap_filtered[2:1212], na.rm=TRUE)
#both_resolved_counts = colSums(gene_both_overlap_resolved[2:1212], na.rm=TRUE)
#both_resolved_gbm_counts = colSums(gene_both_overlap_resolved_gbm[2:1212], na.rm=TRUE)
#both_resolved_tem_counts = colSums(gene_both_overlap_resolved_tem[2:1212], na.rm=TRUE)
umr_gene_counts2 = colSums(accession_gene_calls[,2:1212]=="U", na.rm=TRUE)

pdf(file=paste0(project_id, "_gene_gbm_tem_stats2.pdf"))
print(ggplot(as.data.frame(gbm_gene_counts2)) + geom_histogram(aes(x=gbm_gene_counts2), binwidth = 100))
print(ggplot(as.data.frame(tem_gene_counts2)) + geom_histogram(aes(x=tem_gene_counts2), binwidth = 100))
print(ggplot(as.data.frame(cbind(gbm_gene_counts2, tem_gene_counts2))) + geom_point(aes(x=gbm_gene_counts2, y=tem_gene_counts2)))
print(ggplot(as.data.frame(both_counts2)) + geom_histogram(aes(x=both_counts2), binwidth = 10))
print(ggplot(as.data.frame(cbind(gbm_gene_counts2, both_counts2))) + geom_point(aes(x=gbm_gene_counts2, y=both_counts2)))
print(ggplot(as.data.frame(cbind(tem_gene_counts2, both_counts2))) + geom_point(aes(x=tem_gene_counts2, y=both_counts2)))
#print(ggplot(as.data.frame(both_filtered_counts)) + geom_histogram(aes(x=both_filtered_counts), binwidth = 10))
#print(ggplot(as.data.frame(both_resolved_counts)) + geom_histogram(aes(x=both_resolved_counts), binwidth = 10))
#print(ggplot(as.data.frame(cbind(both_counts, both_filtered_counts))) + geom_point(aes(x=both_counts, y=both_filtered_counts)))
#print(ggplot(as.data.frame(cbind(both_counts, both_resolved_counts))) + geom_point(aes(x=both_counts, y=both_resolved_counts)))
#print(ggplot(as.data.frame(cbind(both_resolved_gbm_counts, both_resolved_tem_counts))) + geom_point(aes(x=both_resolved_gbm_counts, y=both_resolved_tem_counts)))


dev.off()

cor(cbind(gbm_gene_counts, both_counts), use="complete.obs")
#gbm_gene_counts both_counts
#gbm_gene_counts       1.0000000   0.2897491
#both_counts           0.2897491   1.0000000
cor(cbind(tem_gene_counts, both_counts), use="complete.obs")
#tem_gene_counts both_counts
#tem_gene_counts       1.0000000   0.7091566
#both_counts           0.7091566   1.0000000
cor(cbind(tem_gene_counts, gbm_gene_counts), use="complete.obs")
#tem_gene_counts gbm_gene_counts
#tem_gene_counts       1.0000000      -0.1927565
#gbm_gene_counts      -0.1927565       1.0000000
cor(cbind(both_counts, both_filtered_counts), use="complete.obs")
#both_counts both_filtered_counts
#both_counts             1.000000             0.451117
#both_filtered_counts    0.451117             1.000000
cor(cbind(both_counts, both_resolved_counts), use="complete.obs")
#both_counts both_resolved_counts
#both_counts            1.0000000            0.9018878
#both_resolved_counts   0.9018878            1.0000000








# alternative approach to SNPs using hdf5 file
#BiocManager::install("rhdf5")
# No need for this, it's probably already installed by another package
#library(rhdf5)

#SNP_file = "../0-reference/1001_SNP_MATRIX/imputed_snps_binary.hdf5"

#h5ls(SNP_file)
#SNP_accessions = h5read(SNP_file, name="accessions") 
#SNP_positions = h5read(SNP_file, name="positions") 
#SNP_snps = rhdf5::h5read(SNP_file, name="snps") 


# alternative approach uses VCF file of SNP calls 
#install.packages("vcfR")
#library(vcfR)


# Load coverage data for each site to get averages per gene
# coverage data for each CG site with coverage > 0 is available in meth_CG.gr
# need to get total No. CG sites per gene first
CG_site_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
gene_CG_site_counts = as.data.frame(table(queryHits(findOverlaps(gene_ranges, CG_site_ranges))))
colnames(gene_CG_site_counts) = c("gene_no", "site_count")
gene_CG_site_counts$gene_no = as.numeric(levels(gene_CG_site_counts$gene_no)[gene_CG_site_counts$gene_no])
gene_CG_site_counts = cbind(gene_CG_site_counts, "gene_ID"=gff.genes[gene_CG_site_counts$gene_no,]$gene_ID)

#gff.genes = merge(gff.genes, gene_CG_site_counts, by="gene_ID", all=TRUE)

no_samples = 0

# make empty structures for the various things we are accumulating
gene_coverage = NULL

for (this_sample in sample_list$SRA.Accession) {
  cat(paste0(this_sample,"\n"))

  # read in CG methylome
  meth_CG.gr <- readMethylome(FileName=paste0(this_sample,"/",project_id,"_CG_",this_sample,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)
  
  # Set aside genes with no CG sites
  set_aside_no_CG.gr = setdiff(gene_ranges, gene_ranges[unique(subjectHits(findOverlaps(meth_CG.gr, gene_ranges)))])
  cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
  # 335 in first sample
    
  CG_gene_overlaps = findOverlaps(meth_CG.gr, gene_ranges)
  CG_genes = unique(subjectHits(CG_gene_overlaps))

  gene_site_coverage = as.data.frame(cbind(gff.genes[subjectHits(CG_gene_overlaps),"gene_ID"], as.data.frame(meth_CG.gr)$T[queryHits(CG_gene_overlaps)]))
  colnames(gene_site_coverage)=c("gene_ID", "coverage")
  gene_site_coverage$coverage = as.numeric(gene_site_coverage$coverage)
  sample_gene_coverage = aggregate(coverage ~ gene_ID, data=gene_site_coverage, sum, na.rm=TRUE)
  sample_gene_coverage = merge(sample_gene_coverage, gene_CG_site_counts, by="gene_ID", all=TRUE)
  sample_gene_coverage$coverage = sample_gene_coverage$coverage/sample_gene_coverage$site_count
  sample_gene_coverage = sample_gene_coverage[,1:2]

  gene_ranges_CG = gene_ranges[CG_genes]
  cat(paste0("gene_ranges_CG contains ",length(gene_ranges_CG)," segments after trimming segments with no CG sites.\n"))
  
  
  

  # add each of the statistics vectors to the relevant all samples matrix
  if (no_samples == 0) {
    gene_coverage = sample_gene_coverage
  } else {
    gene_coverage = merge(gene_coverage, sample_gene_coverage, by="gene_ID", all=TRUE)
  }   
  
  # These lines not needed for R 4.0 and higher
  #if (doing_non_cg) {
  #gene_meth_CHG[,no_samples+2] = as.numeric(levels(gene_meth_CHG[,no_samples+2])[gene_meth_CHG[,no_samples+2]])
  #gene_meth_CHH[,no_samples+2] = as.numeric(levels(gene_meth_CHH[,no_samples+2])[gene_meth_CHH[,no_samples+2]])
  #}
  #gene_gbm_overlap[,no_samples+2] = as.numeric(levels(gene_gbm_overlap[,no_samples+2])[gene_gbm_overlap[,no_samples+2]])
  #gene_tem_overlap[,no_samples+2] = as.numeric(levels(gene_tem_overlap[,no_samples+2])[gene_tem_overlap[,no_samples+2]])
  #gene_umr_overlap[,no_samples+2] = as.numeric(levels(gene_umr_overlap[,no_samples+2])[gene_umr_overlap[,no_samples+2]])
  #gene_gbm_mCG_sites[,no_samples+2] = as.numeric(levels(gene_gbm_mCG_sites[,no_samples+2])[gene_gbm_mCG_sites[,no_samples+2]])
  #gene_tem_mCG_sites[,no_samples+2] = as.numeric(levels(gene_tem_mCG_sites[,no_samples+2])[gene_tem_mCG_sites[,no_samples+2]])
  
  # These lines are the version for R 4.0 and higher
  gene_coverage[,no_samples+2] = as.numeric(gene_coverage[,no_samples+2])
  colnames(gene_coverage)[no_samples+2] = this_sample
  
  no_samples = no_samples + 1
  processed_samples=c(processed_samples, this_sample)
  #  this_sample_meth_status = data.frame(readRDS(this_file))
  if (no_samples == 1) {
    #    all_samples_meth_status = this_sample_meth_status
  } else {
    #	all_samples_meth_status = merge(all_samples_meth_status, this_sample_meth_status, by=c("Chromosome","Locus","Strand"), all=TRUE)
  }
  #} else {
  #  not_processed_samples=c(not_processed_samples, this_sample)
  #}
}


SORT OUT THE PLOTTING AND CORRELATIONS
# plot distributions of numbers of gbM genes, TEM genes and both per accession
gene_coverage_means = colSums(gene_coverage[,2:1212], na.rm=TRUE)/nrow(gene_coverage)

pdf(file=paste0(project_id, "_gene_gbm_tem_stats.pdf"))
print(ggplot(as.data.frame(gbm_gene_counts)) + geom_histogram(aes(x=gbm_gene_counts), binwidth = 100))
print(ggplot(as.data.frame(tem_gene_counts)) + geom_histogram(aes(x=tem_gene_counts), binwidth = 100))
print(ggplot(as.data.frame(cbind(gbm_gene_counts, tem_gene_counts))) + geom_point(aes(x=gbm_gene_counts, y=tem_gene_counts)))
print(ggplot(as.data.frame(both_counts)) + geom_histogram(aes(x=both_counts), binwidth = 10))
print(ggplot(as.data.frame(cbind(gbm_gene_counts, both_counts))) + geom_point(aes(x=gbm_gene_counts, y=both_counts)))
print(ggplot(as.data.frame(cbind(tem_gene_counts, both_counts))) + geom_point(aes(x=tem_gene_counts, y=both_counts)))
print(ggplot(as.data.frame(both_filtered_counts)) + geom_histogram(aes(x=both_filtered_counts), binwidth = 10))
print(ggplot(as.data.frame(both_resolved_counts)) + geom_histogram(aes(x=both_resolved_counts), binwidth = 10))
print(ggplot(as.data.frame(cbind(both_counts, both_filtered_counts))) + geom_point(aes(x=both_counts, y=both_filtered_counts)))
print(ggplot(as.data.frame(cbind(both_counts, both_resolved_counts))) + geom_point(aes(x=both_counts, y=both_resolved_counts)))
print(ggplot(as.data.frame(cbind(both_resolved_gbm_counts, both_resolved_tem_counts))) + geom_point(aes(x=both_resolved_gbm_counts, y=both_resolved_tem_counts)))


dev.off()

cor(cbind(gbm_gene_counts, both_counts), use="complete.obs")
#gbm_gene_counts both_counts
#gbm_gene_counts       1.0000000   0.2897491
#both_counts           0.2897491   1.0000000
cor(cbind(tem_gene_counts, both_counts), use="complete.obs")
#tem_gene_counts both_counts
#tem_gene_counts       1.0000000   0.7091566
#both_counts           0.7091566   1.0000000
cor(cbind(tem_gene_counts, gbm_gene_counts), use="complete.obs")
#tem_gene_counts gbm_gene_counts
#tem_gene_counts       1.0000000      -0.1927565
#gbm_gene_counts      -0.1927565       1.0000000
cor(cbind(both_counts, both_filtered_counts), use="complete.obs")
#both_counts both_filtered_counts
#both_counts             1.000000             0.451117
#both_filtered_counts    0.451117             1.000000
cor(cbind(both_counts, both_resolved_counts), use="complete.obs")
#both_counts both_resolved_counts
#both_counts            1.0000000            0.9018878
#both_resolved_counts   0.9018878            1.0000000


# Export the matrices to use in GWAS analysis
write.table(gene_coverage, file=paste0(project_id,"_gene_coverage.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#gene_coverage = read.table(file=paste0(project_id, "_gene_coverage.txt"), sep="\t", header=TRUE)




# Correlate variance in methylation with coverage, among median globally methylated accessions
mean_mCG_per_gene_per_sample = read.table(file="mean_mCG_per_gene_per_sample.txt", sep="\t", header=TRUE)

sample_mean_gene_mCG = colMeans(mean_mCG_per_gene_per_sample[,2:ncol(mean_mCG_per_gene_per_sample)], na.rm=TRUE)
mean_sample_gene_mCG = mean(sample_mean_gene_mCG)
sd_sample_gene_mCG = sd(sample_mean_gene_mCG)
# pick out the accessions with mean gene mCG within 1SD of the mean (802 accessions out of 1221)
median_samples = names(sample_mean_gene_mCG[(sample_mean_gene_mCG>mean_sample_gene_mCG - sd_sample_gene_mCG) & (sample_mean_gene_mCG<mean_sample_gene_mCG + sd_sample_gene_mCG)] )

gene_coverage[,c("gene_ID",median_samples)]
mean_mCG_per_gene_per_sample[,c("gene_ID",median_samples)]

coverage_long = melt(as.data.table(gene_coverage), id.vars="gene_ID", measure.vars=median_samples, na.rm=FALSE)
colnames(coverage_long) = c("gene_ID", "sample", "coverage")
mCG_long = melt(as.data.table(mean_mCG_per_gene_per_sample), id.vars="gene_ID", measure.vars=median_samples, na.rm=FALSE)
colnames(mCG_long) = c("gene_ID", "sample", "mean_mCG")
coverage_mCG = merge(coverage_long, mCG_long, by=c("gene_ID", "sample"), all=TRUE)

max_plot_coverage = 20
plot_binsize = 2
page_size = 36
plot_genes = unique(coverage_mCG[coverage_mCG$coverage<=max_plot_coverage,"gene_ID"])
pdf("m1001_gene_mCG_variance_with_coverage.pdf")
for (i in 1:round(nrow(plot_genes)+round(page_size/2)/page_size)) {
cat(paste(i,"\n"))
    gene_list = coverage_mCG[(coverage_mCG$coverage<=max_plot_coverage) & (coverage_mCG$gene_ID %in% plot_genes[((i-1)*page_size+1):(i*page_size),]$gene_ID),]

p = gene_list %>%
  mutate( bin=cut_width(coverage, width=plot_binsize, boundary=0) ) %>%
  ggplot( aes(x=bin, y=mean_mCG) ) +geom_boxplot(coef=1e30) +theme_minimal() +facet_wrap(~gene_ID) +theme_minimal()
print(p)
}
dev.off()

# How many data points are excluded at each coverage threshold?
coverage_hist=hist(coverage_mCG$coverage, xlim=c(1,100), breaks=seq(0,1000000,by=1))


# Disregarding bins with coverage below 6, which bin has the lowest number of accessions (for each gene)?
# Use this level of coverage to sample in each bin, and calculate variance

gene_gbm_calls_long = melt(as.data.table(gene_gbm_calls), id.vars="gene_ID", measure.vars=median_samples, na.rm=FALSE)
colnames(gene_gbm_calls_long) = c("gene_ID", "sample", "gbm_call")
coverage_mCG = merge(coverage_mCG, gene_gbm_calls_long , by=c("gene_ID", "sample"), all=TRUE)

adequate_coverage_threshold = 6
max_adequate_coverage=20
minimum_sufficient_accessions = 10
adequate_coverage_mCG = coverage_mCG[coverage_mCG$coverage>=adequate_coverage_threshold,]
coverage_mCG_stats = NULL
stat_count = 0
for (this_gene in unique(adequate_coverage_mCG[adequate_coverage_mCG$gbm_call==1,]$gene_ID)) {
  bin_sizes = hist(adequate_coverage_mCG[adequate_coverage_mCG$gene_ID==this_gene,]$coverage, xlim=c(0,max_plot_coverage), breaks=seq(0,1000000,by=plot_binsize))
  # adequate_coverage_threshold might not have very many accessions, so step the threshold up until more than minimum_sufficient_accessions are found
  min_bin_size = min(bin_sizes$counts[(bin_sizes$breaks>=adequate_coverage_threshold) & (bin_sizes$breaks<=max_adequate_coverage)], na.rm=TRUE)
  adequate_bins = (bin_sizes$counts >= minimum_sufficient_accessions) & (bin_sizes$counts >= min_bin_size) & (bin_sizes$counts <= max_adequate_coverage)
  cat(paste0(this_gene, " ", min_bin_size, "\n"))
  for (this_bin in bin_sizes$breaks[adequate_bins]) {
    sample_items=sample(seq(1,nrow(adequate_coverage_mCG[(adequate_coverage_mCG$gene_ID==this_gene) & (adequate_coverage_mCG$coverage>(this_bin)) & (adequate_coverage_mCG$coverage<=(this_bin+plot_binsize))])), min_bin_size, replace=FALSE)
    mCG_variance = var(adequate_coverage_mCG[(adequate_coverage_mCG$gene_ID==this_gene) & (adequate_coverage_mCG$coverage>(this_bin)) & (adequate_coverage_mCG$coverage<=(this_bin+plot_binsize)),][sample_items,"mean_mCG"])
    if (stat_count == 0) {
      coverage_mCG_stats = data.frame("gene_ID"=this_gene, "coverage"=this_bin, "mCG_variance"=mCG_variance)
    } else {
      coverage_mCG_stats = rbind.data.frame(coverage_mCG_stats, data.frame("gene_ID"=this_gene, "coverage"=this_bin, "mCG_variance"=mCG_variance))
    }
    stat_count = stat_count + 1
  }
}
                                                                                                                                                         
ggplot(coverage_mCG_stats) + geom_point(aes(x=coverage, y=mean_mCG)) +facet_wrap(~gene_ID, scales="free_y") +xlim(0,20)

page_size = 36
plot_genes = unique(coverage_mCG[coverage_mCG$coverage<=max_plot_coverage,"gene_ID"])
pdf("m1001_gene_mCG_variance_by_coverage.pdf")
for (i in 1:round(nrow(plot_genes)+round(page_size/2)/page_size)) {
  cat(paste(i,"\n"))
  gene_list = coverage_mCG_stats[(coverage_mCG_stats$coverage<=max_plot_coverage) & (coverage_mCG_stats$gene_ID %in% plot_genes[((i-1)*page_size+1):(i*page_size),]$gene_ID),]
  
  #  p = gene_list %>%
  #  mutate( bin=cut_width(coverage, width=plot_binsize, boundary=0) ) %>%
  #  ggplot( aes(x=bin, y=mean_mCG) ) +geom_point() +theme_minimal() +facet_wrap(~gene_ID) +theme_minimal()
  p = ggplot(gene_list, aes(x=coverage, y=mean_mCG) ) +geom_point() +theme_minimal() +facet_wrap(~gene_ID, scales="free_y") +theme_minimal()+xlim(6,20) +ylab("mCG variance")
  print(p)
}
dev.off()


# We now want to generate phenotype matrices at a range of minimum coverage thresholds
# For genes in accessions where coverage is below threshold, we should replace their methylation values with NAs.
# For genes in accessions where the gene contains a TEM segment, we should replace their methylation values with NAs.

for (min_coverage in c(8,10,12,14,16)) {
  pheno_values = NULL
  pheno_values_gbm = NULL
  pheno_values_tem = NULL
  no_samples = 0
  for (this_sample in sample_list$SRA.Accession) {
    cat(paste0(min_coverage," ",this_sample,"\n"))
    gene_stats = merge(cbind("gene_ID"=mean_mCG_per_gene_per_sample[,"gene_ID"], "mean_mCG"=mean_mCG_per_gene_per_sample[,colnames(mean_mCG_per_gene_per_sample)==this_sample]), cbind("gene_ID"=gene_coverage[,"gene_ID"], "coverage"=gene_coverage[,colnames(gene_coverage)==this_sample]), by="gene_ID", all=TRUE)
    gene_stats = merge(gene_stats, cbind("gene_ID"=gene_gbm_calls[,"gene_ID"], "gbm_call"=gene_gbm_calls[,colnames(gene_gbm_calls)==this_sample]), by="gene_ID", all=TRUE)
    gene_stats = merge(gene_stats, cbind("gene_ID"=gene_tem_calls[,"gene_ID"], "tem_call"=gene_tem_calls[,colnames(gene_tem_calls)==this_sample]), by="gene_ID", all=TRUE)
    
    if (ncol(gene_stats)==5) {
      gene_stats[,2] = as.numeric(gene_stats[,2])
      gene_stats[,3] = as.numeric(gene_stats[,3])
      gene_stats[,4] = as.numeric(gene_stats[,4])
      gene_stats[,5] = as.numeric(gene_stats[,5])
      
      gene_stats[is.na(gene_stats$coverage),]$mean_mCG = NA
      gene_stats[(!is.na(gene_stats$coverage)) & (gene_stats$coverage<min_coverage),]$mean_mCG = NA
    
      gene_stats$gbm_mCG = gene_stats$mean_mCG*gene_stats$gbm_call
      gene_stats$tem_mCG = gene_stats$mean_mCG*gene_stats$tem_call
      
      if (no_samples == 0) {
        pheno_values = as.data.frame(gene_stats[,c(1,2)])
        pheno_values_gbm = as.data.frame(gene_stats[,c(1,6)])
        pheno_values_tem = as.data.frame(gene_stats[,c(1,7)])
      } else {
        pheno_values = merge(pheno_values, gene_stats[,c(1,2)], by="gene_ID", all=TRUE)
        pheno_values_gbm = merge(pheno_values_gbm, gene_stats[,c(1,6)], by="gene_ID", all=TRUE)
        pheno_values_tem = merge(pheno_values_tem, gene_stats[,c(1,7)], by="gene_ID", all=TRUE)
      }   
      colnames(pheno_values)[no_samples+2] = this_sample
      colnames(pheno_values_gbm)[no_samples+2] = this_sample
      colnames(pheno_values_tem)[no_samples+2] = this_sample
      
      # These lines not needed for R 4.0 and higher
      #if (doing_non_cg) {
      #gene_both_overlap[,no_samples+2] = as.numeric(levels(gene_both_overlap[,no_samples+2])[gene_both_overlap[,no_samples+2]])
    
      # These lines are the version for R 4.0 and higher
      pheno_values[,no_samples+2] = as.numeric(pheno_values[,no_samples+2])
      pheno_values_gbm[,no_samples+2] = as.numeric(pheno_values_gbm[,no_samples+2])
      pheno_values_tem[,no_samples+2] = as.numeric(pheno_values_tem[,no_samples+2])
      
      no_samples = no_samples + 1
    }
  }
  write.table(pheno_values, file=paste0("mean_mCG_per_gene_per_sample_cov_",min_coverage,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(pheno_values_gbm, file=paste0("mean_mCG_per_gbm_gene_per_sample_cov_",min_coverage,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(pheno_values_gbm, file=paste0("mean_mCG_per_tem_gene_per_sample_cov_",min_coverage,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}


     








# Some samples may need to be excluded from this analysis
excluded_samples = c()
valid_samples=c()
if (project_id=="SRA035939") {
  excluded_samples = c("SRR342380","SRR342390")
} else if (project_id=="PRJEB2678") {
  excluded_samples = c("ERR046563","ERR046564")
}
#valid_samples=rownames(sample_info)[!(rownames(sample_info) %in% excluded_samples)]
# Start from samples_with_genomes - processed_samples for which we have SNP/indel information
valid_samples=samples_with_genomes[!(samples_with_genomes %in% excluded_samples)]
### Need to check that Col-0 is in valid list even though it has 'no genome' as no need to apply SNPs and indels

# This never gets used again so commented out:
#valid_samples_meth_status=all_samples_meth_status[,!names(all_samples_meth_status) %in% excluded_samples]


# Find out which sites are generally always methylated, and which unmethylated, across the genome

# This part sums up methylation calls across lines to give an 'average' methylation level per site irrespective of line
# Partial calls are assigned an arbitrary methylation level of 0.25
average_methylation = matrix(0, nrow=nrow(all_samples_meth_status), ncol=2)
for(sample_name in valid_samples) {
  if(opt$verbose) {cat(paste0("Summing methylation levels in sample ",sample_name,"\n"))}
  average_methylation[,1]=average_methylation[,1]+ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse(all_samples_meth_status[,sample_name]=="U",0,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",0.25,0))))
  average_methylation[,2]=average_methylation[,2]+ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse(all_samples_meth_status[,sample_name]=="U",1,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",1,0))))
}

# Plot distribution of average methylation level of cytosines across samples
hist(average_methylation[,1]/average_methylation[,2], breaks=100)

