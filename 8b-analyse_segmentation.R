#!/usr/bin/env Rscript

# 8-non_CG_methylation.R
#
# Jay Moore, John Innes Centre, Norwich
# jonathan.moore@jic.ac.uk
#
# 05/11/2019
#
# Version 1 - 1001 methylomes data analysis, adapted from Schmitz/Becker project - load non-CG methylation and call status for each sample
#
# Change log
#


# Data source (project_id)
#  .methratio2 files from bsmap


# Summary of functions
#
# Done:
# Read GFF files output from bs_sequel for each sample in a project - each file contains C and T calls at all sites in a meth_context
# Align and merge the C and T calls at sites across samples
# Plot summary data about samples, sites
# Assign methylation status at each site in each sample
#
# This script calls methylation status of a C site based on:
# Fisher's exact test -> "I" (identify sites where coverage is too low to make a clear distinction between methylated and unmethylated status given the conversion efficiency and sequencing and alignment errors for the sample)
# Binomial test with B-H FDH adjusted p-values >= 0.01 -> "U" (identify sites where the number of "C" calls is too few to make a call of methylated for the sample)
# #C >= #T -> differentiate between "M" and "P" (impose an expectation that methylated sites will be represented by at least half of aligned reads being converted, else call Partial methylation)
# Coverage < 5 sigma -> "I" (assume coverage is approxinmately Gaussian distribution, and ignore sites where coverage is too high to be meaningful for the sample)
#
# Once all sites have been called, the script combines adjacent C and G sites identified from the reference and addresses the concordance of calls in adjacent CG dinucleotides.
#
# Bring CG site combined calls into same structure as CHG/CHH calls
#
# Underway:
#
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
	project_id="Mp1"
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
	meth_context="NONCG"
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
  pathroot="X:/Daniel-Zilberman/"
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
#if(!require(mixtools)){
#  install.packages("mixtools")
#  library(mixtools)
#}
if(!require(rootSolve)){
  install.packages("rootSolve")
  library(rootSolve)
}
#if(!require(ggpubr)){
#  install.packages("ggpubr")
#  library(ggpubr)
#}

library(data.table)
library(stringr)
library(ggplot2)
library(plyr)
#library(ggpubr)


# Source directory containing alignments from bs_sequel
source_dir = paste0(pathroot,"/Projects/Jay-1001_methylomes/3-alignments/")

# Working directory where outputs will be directed to
setwd(dir = paste0(pathroot,"/Projects/Jay-1001_methylomes/5-analysis/"))

raw_data_dir = "../0-raw_data/"
reference_dir = "../0-reference/"

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



# Started with global Arabidopsis accessions from Ecker study (927 of)
sample_list1 = read.table(paste0(raw_data_dir,"GSE43857.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=43857
fastq_list1 = read.table(paste0(raw_data_dir,"ENA_SRA065807.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/SRA065807

# Swedish Arabidopsis accessions from separate study (284 of)
sample_list2 = read.table(paste0(raw_data_dir,"GSE54292.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=54292
fastq_list2 = read.table(paste0(raw_data_dir,"ENA_SRP035593.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/PRJNA236110

sample_list = rbind.data.frame(sample_list1, sample_list2)
fastq_list = rbind.data.frame(fastq_list1, fastq_list2)

# The above sample_list routine was somewhat unsatisfactory.  Instead, we will read in the synthesis created manually in Excel, using the above resources as a starting point:
sample_list = read.table(paste0(raw_data_dir,"sample_info_synthesis.txt"), sep="\t", header=TRUE)
sample_list$Name = str_split_fixed(sample_list$Title, "\\(",2)[,1]
#sample_list1$Ecotype_id = as.numeric(str_split_fixed(str_split_fixed(sample_list1$Title, "\\(",2)[,2], "\\)", 2)[,1])
sample_list$Ecotype_id = as.numeric(str_split_fixed(sample_list$Sample, "_", 2)[,1])


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


# get the list of runs for the sample
# concatenate the methylation data for the runs
all_samples_ct=NA
run_no = 0
for (this_run in fastq_list[fastq_list$experiment_accession==sample_id,"run_accession"]) {
	run_no = run_no + 1
	if(opt$verbose) {cat(paste0("Loading ",this_run,"\n"))}
	x=read.table(file=paste0(source_dir, this_run, "_TAIR10.non-cg.methratio2"), header=TRUE, sep="\t")
    # needs sep="\t" because some entries (e.g. first two in mitochondria) have no entry in context string for some reason
	
	# merge the data for this sample with the accumulated data	
	colnames(x)[7]=paste0(this_run,"_C")
	x$T_count = x[,6]-x[,7]
	colnames(x)[10]=paste0(this_run,"_T")
	
	# Added in the context column for meth_context=NONCG because we need to split out CHG and CHH later
	if (run_no==1) {
		all_samples_ct=x[c("chr","pos","strand","context",paste0(this_run,"_C"),paste0(this_run,"_T"))]
	} else {
		all_samples_ct=merge(x[c("chr","pos","strand","context",paste0(this_run,"_C"),paste0(this_run,"_T"))], all_samples_ct, by=c("chr","pos","strand","context"), all=TRUE)
	}
	
	rm(x)
}

# make new columns for the totals
all_samples_ct[,paste0(sample_id,"_C")] = rep(0, nrow(all_samples_ct))
all_samples_ct[,paste0(sample_id,"_T")] = rep(0, nrow(all_samples_ct))
#colnames(all_samples_ct)[4 + run_no*2]=paste0(sample_id,"_C")
#colnames(all_samples_ct)[4 + run_no*2]=paste0(sample_id,"_T")

# accumulate the totals for all runs
run_no = 0
for (this_run in fastq_list[fastq_list$experiment_accession==sample_id,"run_accession"]) {
	run_no = run_no + 1
	if(opt$verbose) {cat(paste0("Accumulating ",this_run,"\n"))}

	all_samples_ct[,paste0(sample_id,"_C")] = all_samples_ct[,paste0(sample_id,"_C")] + ifelse(is.na(all_samples_ct[,paste0(this_run,"_C")]),0,all_samples_ct[,paste0(this_run,"_C")])
	all_samples_ct[,paste0(sample_id,"_T")] = all_samples_ct[,paste0(sample_id,"_T")] + ifelse(is.na(all_samples_ct[,paste0(this_run,"_T")]),0,all_samples_ct[,paste0(this_run,"_T")])
}

colnames(all_samples_ct)[1]="Chromosome"
colnames(all_samples_ct)[2]="Locus"
colnames(all_samples_ct)[3]="Strand"

all_samples_ct = all_samples_ct[,c(1, 2, 3, 4, 5 + run_no*2, 6 + run_no*2)]
#all_samples_ct = data.table(all_samples_ct,key=c("Chromosome", "Locus","Strand"))

	
# Dump the data table to file for a convenience cache
saveRDS(all_samples_ct, file=paste0(project_id,"_",meth_context,"_",sample_id,"_all_samples_ct.rds"))


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
	
    # Find the ecotype corresponding to this sample
    this_ecotype = sample_list[sample_list$SRA.Accession==sample_id, "Ecotype_id"]
    this_snp_file = paste0(reference_dir,this_ecotype,"_C_vars.txt")
    if (file.exists(this_snp_file)) {
	  if(opt$verbose) {cat(paste0("Loading ",this_ecotype," SNPs\n"))}
      #no_ecotypes = no_ecotypes + 1
      #samples_with_genomes=c(samples_with_genomes, this_sample)
	  # read in the known C sites with SNPs for this sample.  Use unique to deduplicate cases where the same site is listed more than once
	  this_sample_snps = unique(read.table(this_snp_file, header=FALSE, sep="\t"))
	  colnames(this_sample_snps) = c("Chromosome", "Locus")
	  #merge snps and sample to get rid of any rows exclusively from snp file, then replace sample values with 0 where a SNP was found
	  mutant_sites = merge(all_samples_ct[,c("Chromosome", "Locus")], cbind(this_sample_snps, "mutation"=rep(1,nrow(this_sample_snps))), by = c("Chromosome", "Locus"), all.x=TRUE, all.y=FALSE) 
      all_samples_ct[, paste0(sample_id,"_C")] = ifelse(is.na(mutant_sites$mutation), as.character(all_samples_ct[, paste0(sample_id,"_C")]), 0)
      all_samples_ct[, paste0(sample_id,"_T")] = ifelse(is.na(mutant_sites$mutation), as.character(all_samples_ct[, paste0(sample_id,"_T")]), 0)
    }


	# Write out a copy of the SNP-masked methylome in the format that MethylSeekR likes to read in (Cov_C+Cov_T, Cov_C)
	#write.table(cbind(paste0(substr(as.character(all_samples_meth_status$Chromosome),1,1),tolower(substr(as.character(all_samples_meth_status$Chromosome),2,3)),substr(as.character(all_samples_meth_status$Chromosome),4,4)),all_reps_meth_status$Locus,across_samples_TM),file=paste0(project_id,"_",meth_context,"_TM_read_counts_across_samples.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	#all_samples_ct = data.frame(all_samples_ct)
    write.table(cbind(as.character(all_samples_ct$Chromosome), format(all_samples_ct$Locus, scientific=FALSE), as.numeric(all_samples_ct[,paste0(sample_id,"_C")])+as.numeric(all_samples_ct[,paste0(sample_id,"_T")]), as.numeric(all_samples_ct[,paste0(sample_id,"_C")])),file=paste0(project_id,"_",meth_context,"_",sample_id,"_TM_read_counts_across_samples.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	
	# Also write out separate files containing the SNP-masked CHG and CHH methylomes
	# First make a vector defining the context of each site
	noncg_context = ifelse((substr(all_samples_ct$context,3,3)=="C" & substr(all_samples_ct$context,5,5)=="G") | (substr(all_samples_ct$context,3,3)=="G" & substr(all_samples_ct$context,1,1)=="C"), "CHG", "CHH")
    write.table(cbind(as.character(all_samples_ct[noncg_context=="CHG","Chromosome"]), format(all_samples_ct[noncg_context=="CHG","Locus"], scientific=FALSE), as.numeric(all_samples_ct[noncg_context=="CHG",paste0(sample_id,"_C")])+as.numeric(all_samples_ct[noncg_context=="CHG",paste0(sample_id,"_T")]), as.numeric(all_samples_ct[noncg_context=="CHG",paste0(sample_id,"_C")])),file=paste0(project_id,"_CHG_",sample_id,"_TM_read_counts_across_samples.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(cbind(as.character(all_samples_ct[noncg_context=="CHH","Chromosome"]), format(all_samples_ct[noncg_context=="CHH","Locus"], scientific=FALSE), as.numeric(all_samples_ct[noncg_context=="CHH",paste0(sample_id,"_C")])+as.numeric(all_samples_ct[noncg_context=="CHH",paste0(sample_id,"_T")]), as.numeric(all_samples_ct[noncg_context=="CHH",paste0(sample_id,"_C")])),file=paste0(project_id,"_CHH_",sample_id,"_TM_read_counts_across_samples.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

	
	# Write out a copy of the SNP-masked methylome in Bismark .cov format for visualising in SeqMonk
	#write.table(cbind(as.character(all_reps_meth_status$Chromosome),all_reps_meth_status$Locus,ifelse(meth_context=="CG",1,0)+all_reps_meth_status$Locus,across_samples_TM[,2]/across_samples_TM[,1],across_samples_TM[,2],across_samples_TM[,1]-across_samples_TM[,2]),file=paste0(project_id,"_",meth_context,"_CT_read_counts_across_samples.cov"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    write.table(cbind(toupper(as.character(all_samples_ct$Chromosome)), str_trim(format(all_samples_ct$Locus, scientific=FALSE)), str_trim(format(all_samples_ct$Locus + ifelse(meth_context=="CG",1,0), scientific=FALSE)), as.numeric(all_samples_ct[,paste0(sample_id,"_C")])/(as.numeric(all_samples_ct[,paste0(sample_id,"_C")])+as.numeric(all_samples_ct[,paste0(sample_id,"_T")])), as.numeric(all_samples_ct[,paste0(sample_id,"_C")]), as.numeric(all_samples_ct[,paste0(sample_id,"_T")])),file=paste0(project_id,"_",meth_context,"_",sample_id,"_CT_read_counts_across_samples.cov"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


	# Read the methylome data in as a Granges object
	#meth.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_",sample_id,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)

	
	### Remove the CpG sites from the SNP loci
	#meth.gr <- removeSNPs(meth.gr, snps.gr)
	
	#plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr1", num.cores=1)
	#plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr2", num.cores=1)
	#plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr3", num.cores=1)
	#plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr4", num.cores=1)
	#plotAlphaDistributionOneChr(m=meth.gr, chr.sel="Chr5", num.cores=1)
	#plotAlphaDistributionOneChr(m=meth.gr, chr.sel="ChrM", num.cores=1)
	#plotAlphaDistributionOneChr(m=meth.gr, chr.sel="ChrC", num.cores=1)

	# Posterior means of alpha<1 indicate a bias towards very low or very high methylation.  Only ChrM has a bimodal distribution, with significant density above 1.0, indicative of likely presence of partially methylated domains (PMDs) (although in our case it is more likely due to multiple copies of differently methylated mitochondrial genomes being analysed together).  Accoringly we skip the next part.

	# Define some loci to mask out from the genome (mitochondria, chloroplast, genes found previously to have coverage anomalies)
	mask_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	mask_loci=rbind(mask_loci,list(Chromosome="ChrM",start_site=1,end_site=sLengths[6]),list(Chromosome="ChrC",start_site=1,end_site=sLengths[7]))

	mask_loci.gr=makeGRangesFromDataFrame(df=mask_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")

	# Define genomicranges representing the whole genome
	genome_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	levels(genome_loci$Chromosome)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
	for (chrom_no in 1:length(sLengths)) {
		genome_loci[chrom_no,]=list(Chromosome=toupper(names(sLengths[chrom_no])),start_site=1,end_site=sLengths[chrom_no])
	}
	genome_loci.gr=makeGRangesFromDataFrame(df=genome_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")

	
	# The following function is from the package MethylSeekR
	# The following changes have been made:
	# nCG.classification threshold on no CG sites to differentiate between LMR and UMR segments replaced by an array of 'median methylation' cutoff levels (may want to change this to use pmeth instead)
	# This was done because in Arabidopsis, there is not a clear separation in segment size between unmethylated and low-methylated segments, but there is a clear separation in methylation level, in contrast to mammals.
	# In Arabidopsis, CG sites tend to be fully methylated, or unmethylated, so using the median values introduces sampling artefacts
	# median.meth replaced by pmeth in scatterplot
	# nCG replaced by nSites for generality across contexts
	segmentUMRsLMRs <- function (m, meth.cutoff = 0.5, nCpG.cutoff = 3, PMDs = NA, pdfFilename = NULL, 
    num.cores = 1, myGenomeSeq, seqLengths, nCpG.smoothing = 3, 
    minCover = 5) 
	{
		#nCG.classification <- 30
		#mMeth.classification=c(0.2,0.4,0.75)
		mMeth.classification=c(0.2)
		message("identifying UMRs and LMRs")
		m = m[values(m)[, 1] >= minCover]
		nCGsPerChr = table(as.character(seqnames(m)))
		chrs = names(nCGsPerChr)[nCGsPerChr >= nCpG.smoothing]
		res <- mclapply(chrs, function(chr) {
			sel <- which(as.character(seqnames(m)) == chr)
			mean.meth <- runmean(Rle(values(m)[sel, 2]/values(m)[sel,1]), k = nCpG.smoothing, endrule = "constant")
			indx <- mean.meth < meth.cutoff
			runValue(indx)[runLength(indx) < nCpG.cutoff & runValue(indx) == TRUE] = FALSE
			runValue(indx)[runLength(indx) < nCpG.cutoff & runValue(indx) == FALSE] = TRUE
			tmp.ids <- rep(1:length(runLength(indx)), runLength(indx))
			tmp <- split(1:length(sel), tmp.ids)
			tmp <- tmp[runValue(indx) == TRUE]
			if (length(tmp) > 0) {
				coords <- cbind(sapply(tmp, min), sapply(tmp, max))
				starts <- round((start(m)[sel[pmax(1, coords[, 1] - 1)]] + start(m)[sel[coords[, 1]]])/2)
				ends <- round((start(m)[sel[coords[, 2]]] + start(m)[sel[pmin(length(sel), coords[, 2] + 1)]])/2)
				hmr.gr = GRanges(seqnames = unique(seqnames(m[sel])), 
					strand = "*", ranges = IRanges(starts, ends), 
					seqlengths = seqLengths)
			}
			else {
				hmr.gr = GRanges(, seqlengths = seqLengths)
			}
			hmr.gr
		}, mc.cores = num.cores)
		segments.gr = do.call(c, unname(res))
		if (class(PMDs) == "GRanges") {
			segments.gr = subsetByOverlaps(segments.gr, PMDs[values(PMDs)$type == "notPMD"])
		}
		nSites = vcountPattern("CG", getSeq(myGenomeSeq, resize(segments.gr, width(segments.gr), fix = "start"), as.character = FALSE))
		ov <- findOverlaps(m, segments.gr)
		T = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), sum)
		M = tapply(values(m[queryHits(ov)])[, 2], subjectHits(ov), sum)
		nSites.segmentation = tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), length)
		median.meth = tapply(as.vector(runmean(Rle(values(m[queryHits(ov)])[,2]/values(m[queryHits(ov)])[, 1]), nCpG.smoothing, endrule = "constant")), subjectHits(ov), median)
		median.meth = pmax(0, median.meth)
		if (!all.equal(as.numeric(names(T)), 1:length(segments.gr))) {
			message("error in calculating methylation levels for PMDs")
		}
		type=ifelse(median.meth<mMeth.classification[1],"UMR", "LMR")
		#type=ifelse(median.meth<mMeth.classification[1],0,ifelse(median.meth<mMeth.classification[2],1,ifelse(median.meth<mMeth.classification[3],2,3)))
		#type = c("UMR", "LMR")[as.integer(nCG < nCG.classification) + 1]
		values(segments.gr) = DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
		jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		if (!is.null(pdfFilename)) {
			pdf(pdfFilename, height = 5, width = 5)
		}
		smoothScatter(log2(values(segments.gr)$nSites), 100 * values(segments.gr)$pmeth, 
			colramp = jet.colors, xlab = "log2 number of sites in segment", 
			ylab = "mean methylation (%)")
		#abline(v = log2(nCG.classification), lty = 5)
		abline(h = 100 * mMeth.classification[1], lty = 5)
		abline(h = 100 * mMeth.classification[2], lty = 5)
		abline(h = 100 * mMeth.classification[3], lty = 5)
		if (!is.null(pdfFilename)) 
			dev.off()
		segments.gr
	}
	
	
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
		if (!all.equal(as.numeric(names(T)), 1:length(segment_model))) {
			message("error in calculating methylation levels for PMDs")
		}	
		# Original code used median.meth to make a 'type' determination for the segment.  We will use pmeth instead
		#type=ifelse(median.meth<mMeth.classification[1],mMeth.classes[1], mMeth.classes[2])
		type=ifelse(M/T<mMeth.classification[1],mMeth.classes[1], mMeth.classes[2])
		#type=ifelse(median.meth<mMeth.classification[1],0,ifelse(median.meth<mMeth.classification[2],1,ifelse(median.meth<mMeth.classification[3],2,3)))
		#type = c("UMR", "LMR")[as.integer(nCG < nCG.classification) + 1]

		DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
	}

	# Function to return NAs for cases where segmentation model doen's overlap any sites in the given context
	NA_by_segment <- function (segment_model)
	{
		segment_count = length(segment_model)
		nSites.segmentation = rep(as.integer(NA),segment_count)
		nSites = rep(as.integer(NA),segment_count)
		T = array(NA,segment_count)
		M = array(NA,segment_count)
		pmeth = array(NA,segment_count)
		median.meth = rep(as.numeric(NA),segment_count)
		type = array(NA,segment_count)

		DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
	}
	
	# Define a corresponding convenience function to save an annotated genomic ranges object
	# Function adapted from MethylSeekR package
	saveUMRLMRSegments <- function (segs, GRangesFilename = NULL, TableFilename = NULL) 
	{
		# Replaced nCG with nSites for generality
		if (!is.null(GRangesFilename)) 
			saveRDS(segs, GRangesFilename)
		if (!is.null(TableFilename)) {
			D = data.frame(chr = as.character(seqnames(segs)), start = start(segs), 
				end = end(segs), type = values(segs)$type, nSites.segmentation = values(segs)$nSites.segmentation, 
				nSites.seq = values(segs)$nSites, mean.meth = 100 * values(segs)$pmeth, 
				median.meth = 100 * values(segs)$median.meth)
			write.table(D, file = TableFilename, quote = FALSE, sep = "\t", 
			row.names = FALSE)
		}
	}
	
	# Define a convenience function to convert chromosome names to mixed-case (e.g. Chr9, ChrC)
	mixedCaseChr <- function(s, strict = FALSE) {
		paste0(toupper(substr(s,1,1)),tolower(substr(s,2,3)),toupper(substr(s,4,4)))
	}

	fixGRChomCase <- function(gr, strict = FALSE) {
    	#mixedCaseChr <- function(s, strict = FALSE) {
	    #	paste0(toupper(substr(s,1,1)),tolower(substr(s,2,3)),toupper(substr(s,4,4)))
	    #}
		levels(gr@seqnames@values) = mixedCaseChr(as.character(levels(gr@seqnames@values)))
		gr@seqinfo@seqnames = levels(gr@seqnames@values)
		gr
	}

	
	# Set chromosome names to mixed case in genome segments object
	levels(genome_loci.gr@seqnames@values) = mixedCaseChr(as.character(levels(genome_loci.gr@seqnames@values)))
	genome_loci.gr@seqinfo@seqnames = levels(genome_loci.gr@seqnames@values)
	levels(mask_loci.gr@seqnames@values) = mixedCaseChr(as.character(levels(mask_loci.gr@seqnames@values)))
	mask_loci.gr@seqinfo@seqnames = levels(mask_loci.gr@seqnames@values)

	
	#methylated_loci.gr=reduce(makeGRangesFromDataFrame(df=rbind(gff.transposons[,c("V1","V4","V5")],gff.genes[gff.genes$m_class=="Heterochromatic",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1"))
	#gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))
	#unmethylated_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))

	# Concatenate the CHG and CHH methylomes (already done), and sort them, by position and chromosome
	# CHG and CHH are already concatenated in our NONCG context file, so we just load it.  Left in sort, sortseqlevels and c in case they do anything useful 
	meth_CHG_CHH.gr <- sort(sortSeqlevels(c(readMethylome(FileName=paste0(project_id,"_",meth_context,"_",sample_id,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths))))


	#methylated_loci.gr=reduce(makeGRangesFromDataFrame(df=rbind(gff.transposons[,c("V1","V4","V5")],gff.genes[gff.genes$m_class=="Heterochromatic",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1"))
	#gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))
	#unmethylated_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))


#	m_n_non_CG_optimisation=data.frame(m=numeric(), n=numeric(), UMR=integer(), LMR=integer(), fp=numeric(), tp=numeric())

	# Loop through values for m and n parameters 
	#m_CHG.seq=seq(0.002,0.1, by=0.002)
#	m_non_CG.seq=seq(0.05,0.5, by=0.05)
	#n_CHG.seq=seq(1,20, by=1)
#	n_non_CG.seq=seq(50,500, by=50)
#	for (m_non_CG.sel in m_non_CG.seq) {
#		for (n_non_CG.sel in n_non_CG.seq) {
#			#n.sel=4
#			#UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_CHG.gr, meth.cutoff=m_CHG.sel, nCpG.cutoff=n_CHG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
#			#UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_CHH.gr, meth.cutoff=m_CHG.sel, nCpG.cutoff=n_CHG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
#			UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_CHG_CHH.gr, meth.cutoff=m_non_CG.sel, nCpG.cutoff=n_non_CG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
#
#			# Capitalise chromosome names in segmentation object
#			levels(UMRLMRsegments_non_CG.gr@seqnames)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
#			UMRLMRsegments_non_CG.gr@seqinfo@seqnames=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
#
#			# Find overlaps between segments and positive and negative 'control' loci (annotated genes and transposons)
#			fp_hits=findOverlaps(methylated_loci.gr, UMRLMRsegments_non_CG.gr)
#			fp_overlaps <- pintersect(methylated_loci.gr[queryHits(fp_hits)], UMRLMRsegments_non_CG.gr[subjectHits(fp_hits)])
#			fp_overlap_prop <- sum(width(fp_overlaps)) / sum(width(methylated_loci.gr))
#			tp_hits=findOverlaps(unmethylated_loci.gr, UMRLMRsegments_non_CG.gr)
#			tp_overlaps <- pintersect(unmethylated_loci.gr[queryHits(tp_hits)], UMRLMRsegments_non_CG.gr[subjectHits(tp_hits)])
#			tp_overlap_prop <- sum(width(tp_overlaps)) / sum(width(unmethylated_loci.gr))
#			
#			# Find overlaps between segments and variable loci
#			#variable_sites_UMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments_CHG.gr[UMRLMRsegments_CHG.gr@elementMetadata$type=="UMR"])
#			#variable_sites_LMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments_CHG.gr[UMRLMRsegments_CHG.gr@elementMetadata$type=="LMR"])
#			m_n_non_CG_optimisation=rbind(m_n_non_CG_optimisation, list(m=m_non_CG.sel, n=n_non_CG.sel, UMR=0, LMR=0, fp=fp_overlap_prop, tp=tp_overlap_prop))
#			# How many variable sites were captured?
#			#cat(paste0("m.sel=",m_CHG.sel," n.sel=",n_CHG.sel," Variable sites in UMRs: ",length(variable_sites_UMR_olaps@from),"  LMRs: ",length(variable_sites_LMR_olaps@from),"\n"))
#		}
#	}
#	m_n_non_CG_optimisation$variable_sites_captured=m_n_non_CG_optimisation$UMR+m_n_non_CG_optimisation$LMR

#	# Plot ROC curves for each value of n, and estimate AUC
#	pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_m_n_non_CG_optimisation_ROC.pdf"))
#	print(ggplot(m_n_non_CG_optimisation, aes(x=fp, y=tp)) + geom_point() +geom_text(aes(label=paste0(m)),hjust=0, vjust=0) +geom_line() +facet_wrap(~ n))
#	dev.off()
	
#	# Estimate AUC from trapezium approximation
#	cat(paste0("Estimated AUC values from ROC curves generated by varying m parameter for each value of n in MethylSeekR:\n"))
#	best_n_non_CG = 0
#	prev_best_auc_non_CG = 0
#	for (n_non_CG.sel in n_non_CG.seq) {
#		m_n_non_CG_optimisation_auc=0
#		prev_tp=0
#		prev_fp=0
#		for (m_non_CG.sel in m_non_CG.seq) {
#			m_n_non_CG_optimisation_auc = m_n_non_CG_optimisation_auc + (1-m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$fp)*(m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$tp-prev_tp) + ((m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$tp-prev_tp)*(m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$fp-prev_fp))/2
#			prev_tp=m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$tp
#			prev_fp=m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$fp
#		}
#		if (m_n_non_CG_optimisation_auc > prev_best_auc_non_CG) {
#			best_n_non_CG = n_non_CG.sel
#			prev_best_auc_non_CG = m_n_non_CG_optimisation_auc
#		}
#		cat(paste0(" n=",n_non_CG.sel," AUC=",m_n_non_CG_optimisation_auc,"\n"))
#	}
#	cat(paste0("n=",best_n_non_CG," maximises AUC (",prev_best_auc_non_CG,").\n"))


# n=50 AUC=0.855976191899861
# n=100 AUC=0.869442588656929
# n=150 AUC=0.874541502003507
# n=200 AUC=0.876019948216069
# n=250 AUC=0.874809539568385
# n=300 AUC=0.873804294288875
# n=350 AUC=0.87223443141152
# n=400 AUC=0.868125233054618
# n=450 AUC=0.866712745995736
# n=500 AUC=0.864528135434016

#n=200 maximises AUC (0.876019948216069).

# Inspection of curve for n=200 (and all the rest) show that around m=0.15 maxes out tp whilst minimising fp


# From Schmitz:
	
 #n=50 AUC=0.853405558876841
 #n=100 AUC=0.880384843230491
 #n=150 AUC=0.890148958302221
 #n=200 AUC=0.896297635230469
 #n=250 AUC=0.90066653773841
 #n=300 AUC=0.904636365469396
 #n=350 AUC=0.907922263804411
 #n=400 AUC=0.910975041031063
 #n=450 AUC=0.913507235932189
 #n=500 AUC=0.916029402508772

 #n=500 maximises AUC (0.916029402508772).
 # best_n_non_CG = 500
 # Inspection of ROC curves indicates that tp rate maxes out (~1.0) around m=0.15, with fp rate minimised at .25 for n=500.

 
 
	# Segment using best_n and m=0.7 to find appropriate cutoff for m.sel (0.9 gave weird results - almost no segments)
#	UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_CHG_CHH.gr, meth.cutoff=0.7, nCpG.cutoff=best_n_non_CG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_", sample_id, "_MethylSeekR_non_CG_segmentation_landscape_",best_n_non_CG,"_0_7.pdf"))

#	saveUMRLMRSegments(segs=UMRLMRsegments_non_CG.gr, GRangesFilename=paste0("UMRsLMRs_non_CG_0.7_",best_n_non_CG,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_", sample_id, "_MethylSeekR_segments_non_CG_0.7_",best_n_non_CG,".tsv"))
	# UMRLMRsegments_non_CG.gr <- readRDS(paste0("UMRsLMRs_non_CG_0.7_",best_n_non_CG,".gr.rds"))
	
	
	# More formally:
	# Plot density distribution of segment mCH*
#	pdf(paste0(project_id,"_",meth_context,"_", sample_id, "_MethylSeekR_non_CG_segmentation_landscape_",best_n_non_CG,"_0_7_density.pdf"))
#	print(ggplot(as.data.frame(UMRLMRsegments_non_CG.gr@elementMetadata@listData), aes(x=pmeth)) + geom_histogram(binwidth=0.01))
#	dev.off()

	# Visual inspection of m1001_NONCG_SRX445897_MethylSeekR_non_CG_segmentation_landscape_200_0_7.pdf shows at least two seperate populations - unmethylated segments with n>2^10 and m<0.05; and methylation segments with 1<n<2^12 and m>=0.05
	#From Schmitz:  Visual inspection of SRA035939_CHH_MethylSeekR_non_CG_segmentation_landscape_500_0_9.pdf shows three clear populations of segments.  Short methylated segments (n<2^5, m>0.15), and longer unmethylated and low-methylated segments (n>2^10, m<0.1; m>0.1).

	# Fit mixture of Gaussians to nonCG density to identify cutoff
	
	# Initially tried fitting 3 Gaussians to whole data set, but third Gaussian was used to cover the scattering of high-pMeth sgments (m>0.3).  Accordingly, fit to segments with m<0.3 only
	
	#install.packages("MASS")  # for fitting negative binomial models to data
	#library(MASS)
#	install.packages("mixtools")
#	library(mixtools)
#	m_non_CG_mixmdl = normalmixEM(as.data.frame(UMRLMRsegments_non_CG.gr)[(!is.na(as.data.frame(UMRLMRsegments_non_CG.gr)$pmeth)) & (as.data.frame(UMRLMRsegments_non_CG.gr)$pmeth<0.3),]$pmeth, k=3)
	# Had to exclude sites with NA CHG methylation or the fitting algorithm exploded	
	
	# Print the characteristics of the fit curves
#	cat(paste0("Fit lambda values (mixture proportions):\n"))
#	m_non_CG_mixmdl$lambda
	# Schmitz data:   
	# Becker data: 
	
#	cat(paste0("Fit mu values (means):\n"))
#	m_non_CG_mixmdl$mu
	# Schmitz data:
	# Becker data:  
	
#	cat(paste0("Fit sigma values (st.devs):\n"))
#	m_non_CG_mixmdl$sigma
	# Schmitz data: 
	# Becker data:  
	
	# Find the intersection of the two larger components in the mixture - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
#	library(rootSolve)
#	mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
#	dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

#	m_non_CG_low_model=2
#	m_non_CG_high_model=m_non_CG_low_model+1
#	segment_m_non_CG_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=m_non_CG_mixmdl$mu[m_non_CG_low_model], sd1=m_non_CG_mixmdl$sigma[m_non_CG_low_model], m2=m_non_CG_mixmdl$mu[m_non_CG_high_model], sd2=m_non_CG_mixmdl$sigma[m_non_CG_high_model], p1=m_non_CG_mixmdl$lambda[m_non_CG_low_model], p2=m_non_CG_mixmdl$lambda[m_non_CG_high_model])
#	segment_m_non_CG_cutoff
    # Most heavily methylated accession SRX445897 0.045130730 but 3 Gaussians doesn't look like great fit
	# Schmitz data: 0.05920002 but 3 Gaussians doesn't look like a great fit
	# Becker data: 

#	pdf(paste0(project_id,"_",meth_context,"_", sample_id, "_m_non_CG_segmentation_m_non_CG_cutoff_density_fitted.pdf"))
#	print(plot(m_non_CG_mixmdl,which=2, breaks=seq(0,1,0.001)))
#	dev.off()

# 	m_non_CG.sel = segment_m_non_CG_cutoff[length(segment_m_non_CG_cutoff)]
	# Run optimised segmentation
#	UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_CHG_CHH.gr, meth.cutoff=m_non_CG.sel, nCpG.cutoff=best_n_non_CG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

	# Save segmentation
#	saveUMRLMRSegments(segs=UMRLMRsegments_non_CG.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".tsv"))

 
 
	# Fitting Gaussians to separate the components of the pmeth distribution was generally unsuccessful.  It was able to identify the cutoff around pmeth=0.067, but this led to highly fragmented unmethylated regions, with very many very short methylated segments that were unconvincing.  Manually setting the cutoff to pmeth=0.15 in the segmentation, on the basis of visual inspection of the methylation landscape, produced a much better segmentation result.  We will manually set best_n_non_CG=200

	m_non_CG.sel=0.15
    best_n_non_CG=200
	
	# Clean up memory space
#	rm(UMRLMRsegments_non_CG.gr)
#	gc()
	
	UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_CHG_CHH.gr, meth.cutoff=m_non_CG.sel, nCpG.cutoff=best_n_non_CG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_non_CG_segmentation_landscape_",best_n_non_CG,"_",m_non_CG.sel,".pdf"))
	pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_CG_segmentation_landscape_",best_n_non_CG,"_",m_non_CG.sel,"_density.pdf"))
	print(ggplot(as.data.frame(UMRLMRsegments_non_CG.gr@elementMetadata@listData), aes(x=pmeth)) + geom_histogram(binwidth=0.01))
	dev.off()

	# Save segmentation
	saveUMRLMRSegments(segs=UMRLMRsegments_non_CG.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".tsv"))
	
	
	#MG_segments.gr = setdiff(genome_loci.gr, mask_loci.gr)
	m_non_CG_segments.gr = setdiff(setdiff(genome_loci.gr, UMRLMRsegments_non_CG.gr), mask_loci.gr)
	UMR_non_CG_segments.gr = setdiff(UMRLMRsegments_non_CG.gr, mask_loci.gr)

 
 	# Switch to CG context to complete the segmentation
	meth_context="CG"
	
	# Free up some space
	rm(meth_CHG_CHH.gr)
	gc()
	
	
	# get the list of runs for the sample
    # concatenate the methylation data for the runs
    all_samples_ct=NA
    run_no = 0
    for (this_run in fastq_list[fastq_list$experiment_accession==sample_id,"run_accession"]) {
	  run_no = run_no + 1
	  if(opt$verbose) {cat(paste0("Loading ",this_run,"\n"))}
	  x=read.table(file=paste0(source_dir, this_run, "_TAIR10.cg.methratio2"), header=TRUE, sep="\t")
      # needs sep="\t" because some entries (e.g. first two in mitochondria) have no entry in context string for some reason
	
	  # merge the data for this sample with the accumulated data	
	  colnames(x)[7]=paste0(this_run,"_C")
	  x$T_count = x[,6]-x[,7]
	  colnames(x)[10]=paste0(this_run,"_T")
	
      if (run_no==1) {
	  	  all_samples_ct=x[c("chr","pos","strand",paste0(this_run,"_C"),paste0(this_run,"_T"))]
	  } else {
		all_samples_ct=merge(x[c("chr","pos","strand",paste0(this_run,"_C"),paste0(this_run,"_T"))], all_samples_ct, by=c("chr","pos","strand"), all=TRUE)
	  }	  
	  rm(x)
    }

    # make new columns for the totals
    all_samples_ct[,paste0(sample_id,"_C")] = rep(0, nrow(all_samples_ct))
    all_samples_ct[,paste0(sample_id,"_T")] = rep(0, nrow(all_samples_ct))
    #colnames(all_samples_ct)[4 + run_no*2]=paste0(sample_id,"_C")
    #colnames(all_samples_ct)[4 + run_no*2]=paste0(sample_id,"_T")

    # accumulate the totals for all runs
    run_no = 0
    for (this_run in fastq_list[fastq_list$experiment_accession==sample_id,"run_accession"]) {
	  run_no = run_no + 1
	  if(opt$verbose) {cat(paste0("Accumulating ",this_run,"\n"))}

  	  all_samples_ct[,paste0(sample_id,"_C")] = all_samples_ct[,paste0(sample_id,"_C")] + ifelse(is.na(all_samples_ct[,paste0(this_run,"_C")]),0,all_samples_ct[,paste0(this_run,"_C")])
	  all_samples_ct[,paste0(sample_id,"_T")] = all_samples_ct[,paste0(sample_id,"_T")] + ifelse(is.na(all_samples_ct[,paste0(this_run,"_T")]),0,all_samples_ct[,paste0(this_run,"_T")])
    }

    colnames(all_samples_ct)[1]="Chromosome"
    colnames(all_samples_ct)[2]="Locus"
    colnames(all_samples_ct)[3]="Strand"

    all_samples_ct = all_samples_ct[,c(1, 2, 3, 4 + run_no*2, 5 + run_no*2)]
    #all_samples_ct = data.table(all_samples_ct,key=c("Chromosome", "Locus","Strand"))

	
    # Dump the data table to file for a convenience cache
    saveRDS(all_samples_ct, file=paste0(project_id,"_",meth_context,"_",sample_id,"_all_samples_ct.rds"))

	### Read in the SNPs data
	#snps.gr <- readSNPTable(FileName=snpFname, seqLengths=sLengths)
    # We will set sample_C and sample_T to 0 for all CG loci where there is an annotated mutation
	
    # Find the ecotype corresponding to this sample
    this_ecotype = sample_list[sample_list$SRA.Accession==sample_id, "Ecotype_id"]
    this_snp_file = paste0(reference_dir,this_ecotype,"_CG_vars.txt")
    if (file.exists(this_snp_file)) {
	  if(opt$verbose) {cat(paste0("Loading ",this_ecotype," SNPs\n"))}
      #no_ecotypes = no_ecotypes + 1
      #samples_with_genomes=c(samples_with_genomes, this_sample)
	  # read in the known C sites with SNPs for this sample.  Use unique to deduplicate cases where the same site is listed more than once
	  this_sample_snps = unique(read.table(this_snp_file, header=FALSE, sep="\t"))
	  colnames(this_sample_snps) = c("Chromosome", "Locus")
	  #merge snps and sample to get rid of any rows exclusively from snp file, then replace sample values with 0 where a SNP was found
	  mutant_sites = merge(all_samples_ct[,c("Chromosome", "Locus")], cbind(this_sample_snps, "mutation"=rep(1,nrow(this_sample_snps))), by = c("Chromosome", "Locus"), all.x=TRUE, all.y=FALSE) 
      all_samples_ct[, paste0(sample_id,"_C")] = ifelse(is.na(mutant_sites$mutation), all_samples_ct[, paste0(sample_id,"_C")], 0)
      all_samples_ct[, paste0(sample_id,"_T")] = ifelse(is.na(mutant_sites$mutation), all_samples_ct[, paste0(sample_id,"_T")], 0)
    }


	# Write out a copy of the SNP-masked methylome in the format that MethylSeekR likes to read in (Cov_C+Cov_T, Cov_C)
	#all_samples_ct = data.frame(all_samples_ct)
    write.table(cbind(as.character(all_samples_ct$Chromosome), format(all_samples_ct$Locus, scientific=FALSE), as.numeric(all_samples_ct[,paste0(sample_id,"_C")])+as.numeric(all_samples_ct[,paste0(sample_id,"_T")]), as.numeric(all_samples_ct[,paste0(sample_id,"_C")])),file=paste0(project_id,"_",meth_context,"_",sample_id,"_TM_read_counts_across_samples.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	
	# Write out a copy of the SNP-masked methylome in Bismark .cov format for visualising in SeqMonk
	write.table(cbind(toupper(as.character(all_samples_ct$Chromosome)), str_trim(format(all_samples_ct$Locus, scientific=FALSE)), str_trim(format(all_samples_ct$Locus + ifelse(meth_context=="CG",1,0), scientific=FALSE)), as.numeric(all_samples_ct[,paste0(sample_id,"_C")])/(as.numeric(all_samples_ct[,paste0(sample_id,"_C")])+as.numeric(all_samples_ct[,paste0(sample_id,"_T")])), as.numeric(all_samples_ct[,paste0(sample_id,"_C")]), as.numeric(all_samples_ct[,paste0(sample_id,"_T")])),file=paste0(project_id,"_",meth_context,"_",sample_id,"_CT_read_counts_across_samples.cov"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

	
	# Read the CG methylome data in as a Granges object
	meth_CG.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_",sample_id,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)

	# Identify the parts of the mCGome which overlap with unmethylated mnonCGome
	potential_CG_segments.gr = meth_CG.gr[unique(subjectHits(findOverlaps(UMR_non_CG_segments.gr, meth_CG.gr)))]

 	# Try m.sel=0.95 to check out whole methylome landscape
	m.sel=0.95
	n.sel=1  # We segment using minimum segment size of 1 so that we can extend our mCHG segments by a single CG site only, if necessary (minimise likelihood of extending a mCHG segment into an adjacent mCG but non-mCHG segment).
	UMRLMRsegments_CG.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_CHG_segmentation_landscape_",n.sel,"_",m.sel,".pdf"))
	
	# Save a copy for visualisation
	saveUMRLMRSegments(segs=UMRLMRsegments_CG.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_",m.sel,"_",n.sel,".tsv"))

	# Visual inspection of landscape shows clear cutoff at 0.4 between methylated and unmethylated segments among non-mCHG loci
	m.sel=0.4
	n.sel=1  # We segment using minimum segment size of 1 so that we can extend our m_non_CG segments by a single CG site only, if necessary (minimise likelihood of extending a m_non_CG segment into an adjacent mCG but non-m_non_CG segment).
	UMRLMRsegments_CG.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=UMRLMRsegments_CG.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_",m.sel,"_",n.sel,".tsv"))

	# Find the mCG segments, rather than the unmethylated, and filter by masked genome:
	# First numbers = with 'Unknown' genes masked out, second numbers with only CHRM and CHRC masked
	mCG_segments_1.gr = setdiff(setdiff(genome_loci.gr, mask_loci.gr),UMRLMRsegments_CG.gr)
	length(mCG_segments_1.gr)
    # , 51769
	# Schmitz: 41863, 41808 mCG segments

	# Adjust segments to remove any segments which don't overlap with any CG sites
	mCG_segments_1.gr = mCG_segments_1.gr[unique(subjectHits(findOverlaps(meth_CG.gr, mCG_segments_1.gr)))]
	length(mCG_segments_1.gr)
	# , 51765
	# Schmitz: 41844, 41805 mCG segments
	
	# Add back the mCG values for the mCG segments, so it will save properly
	values(mCG_segments_1.gr) = cbind(values(mCG_segments_1.gr),meth_by_segment(meth_CG.gr, segment_model=mCG_segments_1.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=mCG_segments_1.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_mCG_masked_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_mCG_masked_",m.sel,"_",n.sel,".tsv"))
	
	# Find the mCG segments which overlap with m_non_CG segments:
	#mCHG_mCG_segments.gr=mCG_segments_1.gr[unique(subjectHits(findOverlaps(setdiff(setdiff(genome_loci.gr, mask_loci.gr),UMRLMRsegments_CHG.gr), mCG_segments_1.gr)))]
	m_non_CG_mCG_segments.gr=mCG_segments_1.gr[unique(subjectHits(findOverlaps(m_non_CG_segments.gr, mCG_segments_1.gr)))]
	length(m_non_CG_mCG_segments.gr)
	#, 1673
	# Schmitz: 2683, 2584 mCG segments have an overlap with m_non_CG segments

	# Add back the mCG values for the mCG overlap segments, so it will save properly
	values(m_non_CG_mCG_segments.gr) = cbind(values(m_non_CG_mCG_segments.gr),meth_by_segment(meth_CG.gr, segment_model=m_non_CG_mCG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=m_non_CG_mCG_segments.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_non_CG_CG_overlap_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_non_CG_CG_overlap_",m.sel,"_",n.sel,".tsv"))

	# Take the Union of the m_non_CG best fit segments (subtracted from masked genome) and those mCG 1 0.4 segments which overlap them, to extend the m_non_CG best fit segments to the next closest CG site where appropriate:
	m_non_CG_segments.gr = reduce(union(setdiff(setdiff(genome_loci.gr, UMR_non_CG_segments.gr), mask_loci.gr), m_non_CG_mCG_segments.gr))
	length(m_non_CG_segments.gr)
	# 2364
	# Schmitz: 8492, 8423 segments
	
	# Adjust segments to remove any segments which don't overlap with any CG sites
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, m_non_CG_segments.gr)))]
	length(m_non_CG_segments.gr)
	# 2360
	# Schmitz: 8286, 8196 m_non_CG_segments

	# Add back the mCG values for the mCHG overlap segments, so it will save properly
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_CG.gr, segment_model=m_non_CG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=m_non_CG_segments.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_non_CG_segments.gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_non_CG_.tsv"))
	
	# Rebuild the non-m_non_CG fraction of the masked genome for subsequent segmentation in the CG context
	potential_CG_segments.gr = meth_CG.gr[unique(subjectHits(findOverlaps(setdiff(setdiff(genome_loci.gr, m_non_CG_segments.gr), mask_loci.gr), meth_CG.gr)))]

	### This is where we may have optimised the minimum segment length for segmenting the mCGome.  This code moved to separate script, and the numbers set here:
	
	# Schmitz numbers:
	best_n_CG = 9
	prev_best_auc = 0.936370415538498
			
			
	# Execute the preferred segmentation model, visualise and save the results
	#UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=0.012, nCpG.cutoff=9, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape.pdf"))

	
	# Here are a couple of other ideas for setting cutoffs between U and M, especially for short segments:
	#m.sel=0.025
	# Set m.sel to 99ile of C/(C+T) for "U" sites
#	m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.99)
#	cat(paste0("mCG proportion which captures 99% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 99% of sites 'U' across all samples is 0.0769230769230769.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 99% of sites 'U' across all samples is 0.0833333333333333.

	# Using pmeth, rather than median.meth with .99 threshold for unmethylated segments produces some bleed into UMR from SMR particularly.  .95 is a better cutoff:
#	m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.95)
#	cat(paste0("mCG proportion which captures 95% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 95% of sites 'U' across all samples is 0.0377358490566038.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 95% of sites 'U' across all samples is 0.0384615384615385.

	# Finally, we chose m.sel iteratively.  By building a finished model of long, medium and short segments in CG context, and assessing density distribution of mCG per segment, we identified a lot of short, lowly methylated segments, that are not credible.  We fit a mixture of Gaussians to this distribution and use their crossing point to set the mCG threshold here:
	
	m.sel=0.09309094
	
	# Segment mCG using the minimum segment which maximises AUC for ROC curve
	n.sel=best_n_CG
	
	UMRLMRsegments.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segmentation_landscape_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_CG_",m.sel,"_",n.sel,".tsv"))

 	# We observe that n.sel=9 generates reasonable segments, not breaking up the majority of GBM loci.  However, it does miss short isolated regions of mCG which look like real GBM.
	# As a strategy to address this we carry out a new segmentation with the same m.sel, but with n.sel=3 then again with n.sel=1
	# We will use these finer segmentations to split Unmethylated segments, where necessary, but will not use resulting short unmethylated segments to further divide mCG segments.

 	# Segment using minimum segment length of 3
	n.sel=3
	UMRLMRsegments_medium.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segmentation_landscape_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments_medium.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_CG_",m.sel,"_",n.sel,".tsv"))

	# Segment using minimum segment length of 1
	n.sel=1
	UMRLMRsegments_fine.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segmentation_landscape_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments_fine.gr, GRangesFilename=paste0(sample_id,"UMRsLMRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_MethylSeekR_segments_CG_",m.sel,"_",n.sel,".tsv"))

 	# We have generated Unmethylated Regions with minimum lengths of 9, 3 and 1 CG sites respectively (UMR9, UMR3 and UMR1)
	# Masked Genome (MG) = Genome - Mask
	# Long Methylated regions (LMR) = MG - (UMR9 + mCHG)
	# Medium Methylated regions (MMR) = (MG - (UMR3 + mCHG)) - LMR
	# Short methylated regions (SMR) = (MG - (UMR1 + mCHG)) - (LMR + MMR)
	# Unmethylated regions (UMR) = UMR9 - (SMR + MMR +mCHG)
	# Methylated regions (MR) = MG - UMR

	# Load parental CHG and CHH methylomes separately, ready to use for annotating segments
	meth_CHG.gr <- readMethylome(FileName=paste0(project_id,"_CHG_",sample_id,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)
	# Now the SLOW bit
	meth_CHH.gr <- readMethylome(FileName=paste0(project_id,"_CHH_",sample_id,"_TM_read_counts_across_samples.tsv"), seqLengths=sLengths)

	
	MG_segments.gr = setdiff(genome_loci.gr, mask_loci.gr)
	#LMR_segments.gr = setdiff(MG_segments.gr, UMRLMRsegments.gr)
	#MMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_medium.gr), LMR_segments.gr)
	#SMR_segments.gr = setdiff(setdiff(MG_segments.gr, UMRLMRsegments_fine.gr), c(LMR_segments.gr, MMR_segments.gr))
	#UMR_segments.gr = setdiff(UMRLMRsegments.gr, c(SMR_segments.gr, MMR_segments.gr))
	#MR_segments.gr = setdiff(MG_segments.gr, UMR_segments.gr)
	MG_non_CG_segments.gr = setdiff(setdiff(genome_loci.gr, m_non_CG_segments.gr), mask_loci.gr)
	#LMR_segments.gr = setdiff(MG_non_CG_segments.gr, UMRLMRsegments.gr)
	LMR_segments.gr = setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments.gr), m_non_CG_segments.gr), mask_loci.gr)
	#MMR_segments.gr = setdiff(setdiff(MG_non_CG_segments.gr, UMRLMRsegments_medium.gr), LMR_segments.gr)
	MMR_segments.gr = setdiff(setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments_medium.gr), m_non_CG_segments.gr), mask_loci.gr), LMR_segments.gr)
	#SMR_segments.gr = setdiff(setdiff(MG_non_CG_segments.gr, UMRLMRsegments_fine.gr), c(LMR_segments.gr, MMR_segments.gr))
	SMR_segments.gr = setdiff(setdiff(setdiff(setdiff(setdiff(genome_loci.gr, UMRLMRsegments_fine.gr), m_non_CG_segments.gr), mask_loci.gr), LMR_segments.gr), MMR_segments.gr)
	#UMR_segments.gr = setdiff(UMRLMRsegments.gr, c(SMR_segments.gr, MMR_segments.gr))
	#UMR_segments.gr = setdiff(setdiff(setdiff(setdiff(UMRLMRsegments.gr, SMR_segments.gr), MMR_segments.gr), mCHG_segments.gr), mask_loci.gr)
	#MR_segments.gr = setdiff(MG_non_CG_segments.gr, UMR_segments.gr)

	# Adjust segments to set aside any segments which don't overlap with any CG sites. We will deal with them later
	# Counts here are first, using Fishers p<0.01, Binomial p<0.01; second, using  Fishers p<0.05, Binomial p<0.005; third, masking mCHG segments prior to segmenting mCG; fourth, segmenting CHG and CHH together, then masking; fifth masking only CHRM ande CHRC
	set_aside_no_CG.gr = c(setdiff(LMR_segments.gr, LMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, LMR_segments.gr)))]), setdiff(MMR_segments.gr, MMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, MMR_segments.gr)))]), setdiff(SMR_segments.gr, SMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, SMR_segments.gr)))]))
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 5.  Schmitz: 45
	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 5.  Schmitz: 1362, 1065, 1065, 1029, 5
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 17961. Schmitz: 16500, 16909, 13659, 12401, 12419
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 4297.  Schmitz: 5294, 5315, 5898, 5053, 5069
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 9057.  Schmitz: 5635, 8385, 9949, 7955, 7648
	#MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, MR_segments.gr)))]
	#cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 27429, 30609, 29506, 25106

	
	# Annotate segments with CG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_CG.gr, segment_model=MG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_CG.gr, segment_model=LMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_CG.gr, segment_model=MMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MMR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_CG.gr, segment_model=SMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","SMR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_CG.gr, segment_model=MR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#replace values previously added: 
	values(m_non_CG_segments.gr) = cbind(meth_by_segment(meth_CG.gr, segment_model=m_non_CG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m_non_CG.sel), mMeth.classes=c("UMR","TEM")))

	# Adjust SMRs to remove any with too low mCG
	# Visual observation of the distribution of mCG density, and good sense (we are averaging over six samples) suggests SMRs with mCG<0.25 are not credible
	SMR_segments.gr = SMR_segments.gr[ as.vector((!is.na(SMR_segments.gr$pmeth)) & (SMR_segments.gr$pmeth>0.25)) ]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after removing segments with mCG<0.25.\n"))
	# 2993. Schmitz: 4619, 5255, 4570, 4562, 4564

	# Build set of UMRs corresponding to new set of LMRs, MMRs, SMRs, TEMs, set_asides:
	UMR_segments.gr = setdiff(setdiff(setdiff(setdiff(setdiff(setdiff(genome_loci.gr, SMR_segments.gr), MMR_segments.gr), LMR_segments.gr), m_non_CG_segments.gr), set_aside_no_CG.gr), mask_loci.gr)
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments.\n"))
	# 25988. Schmitz: 26495, 26255, 26178

	# Remove any UMR segments containing no CG sites, after setting them aside
	set_aside_no_CG.gr = c(set_aside_no_CG.gr, setdiff(UMR_segments.gr, UMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, UMR_segments.gr)))]))
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 25988.  Schmitz: 27188, 30383, 29026, 26389, 26137
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 5.  Schmitz: 86
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_CG.gr, segment_model=UMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	
	# Adjust segments to remove any segments which don't overlap with any CHG sites, after setting them aside
	set_aside_no_CHG.gr = c(setdiff(LMR_segments.gr, LMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, LMR_segments.gr)))]), setdiff(MMR_segments.gr, MMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, MMR_segments.gr)))]), setdiff(SMR_segments.gr, SMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, SMR_segments.gr)))]))
	cat(paste0("set_aside_no_CHG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 661. Schmitz: 1099

	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 5  Schmitz: 1344, 1047, 1047, 1012, 5
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 17922 Schmitz: 16483, 16892, 13621, 12304, 12329
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 4180  Schmitz: 5205, 5221, 5750, 4873, 4893
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 2488  Schmitz: 4677, 6672, 7964, 3797, 3731
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 25968  Schmitz: 27185, 30374, 29000, 26268, 26087
	#MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, MR_segments.gr)))]
	#cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 26365, 28785, 27335, 23169, 
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, m_non_CG_segments.gr)))]
	cat(paste0("m_non_CG_segments.gr contains ",length(m_non_CG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 2360  Schmitz:   ,  , 4437, 8193, 8117
	
	# Annotate segments with CHG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=MG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=LMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=MMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=SMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=UMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=MR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=m_non_CG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))


	# Adjust segments to remove any segments which don't overlap with any CHH sites, after setting them aside
	set_aside_no_CHH.gr = c(setdiff(LMR_segments.gr, LMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, LMR_segments.gr)))]), setdiff(MMR_segments.gr, MMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, MMR_segments.gr)))]), setdiff(SMR_segments.gr, SMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, SMR_segments.gr)))]))
	cat(paste0("set_aside_no_CHG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 661  Schmitz: 1099

	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 5  Schmitz: 1012
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 17922  Schmitz: 12328
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 4176  Schmitz: 4890
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 2471  Schmitz: 3720
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 25967  Schmitz: 26087
	#MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, MR_segments.gr)))]
	#cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 23116
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, m_non_CG_segments.gr)))]
	cat(paste0("m_non_CG_segments.gr contains ",length(m_non_CG_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 2360  Schmitz: 8117
	
	# Annotate segments with CHH methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=MG_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=LMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=MMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=SMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=UMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_.gr, segment_model=MR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=m_non_CG_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Merge the various distinct classes of segment
	segmentation_model.gr = c(LMR_segments.gr, MMR_segments.gr, SMR_segments.gr, UMR_segments.gr, m_non_CG_segments.gr)
	length(segmentation_model.gr)
	# 52896  Schmitz: 55873, 55142

	# Define a convenience function to save a segmentation model. This is adapted from MethylSeekR package
	saveSegmentationModel = function (segs, GRangesFilename = NULL, TableFilename = NULL) 
	{
		if (!is.null(GRangesFilename)) 
			saveRDS(segs, GRangesFilename)
		if (!is.null(TableFilename)) {
			D = data.frame(chr = as.character(seqnames(segs)), start = start(segs), 
				end = end(segs), type = values(segs)$type)
        write.table(D, file = TableFilename, quote = FALSE, sep = "\t", 
            row.names = FALSE)
		}
	}
	
#	saveSegmentationModel(segs=LMR_segments.gr, GRangesFilename=paste0(sample_id,"LMRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_LMR_segments_CG_",m.sel,"_",n.sel,".tsv"))
#	saveSegmentationModel(segs=MMR_segments.gr, GRangesFilename=paste0(sample_id,"MMRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MMR_segments_CG_",m.sel,"_",n.sel,".tsv"))
#	saveSegmentationModel(segs=SMR_segments.gr, GRangesFilename=paste0(sample_id,"SMRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_SMR_segments_CG_",m.sel,"_",n.sel,".tsv"))
#	saveSegmentationModel(segs=UMR_segments.gr, GRangesFilename=paste0(sample_id,"UMRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_UMR_segments_CG_",m.sel,"_",n.sel,".tsv"))
#	saveSegmentationModel(segs=MR_segments.gr, GRangesFilename=paste0(sample_id,"MRs_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MR_segments_CG_",m.sel,"_",n.sel,".tsv"))

	table(segmentation_model.gr@elementMetadata@listData$type)

    #LMR   MMR    MR   SMR   TEM   UMR 
    #17630  4068    21  2471  2352 26353 

	# 21 (Schmitz 4,1) segments ended up with type="MR" (short segments with relatively high levels of methylation) so we set these to be "SMR":
	segmentation_model.gr@elementMetadata@listData$type[(!is.na(segmentation_model.gr@elementMetadata@listData$type)) & (segmentation_model.gr@elementMetadata@listData$type=="MR")]="SMR"

	
	# Show a table of counts of segments by type

#	table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$type)
	#Schmitz data:
	#LMR   MMR   SMR   TEM   UMR 
	#11813  4850  3733  7823 27519 
	#12081  4843  3747  7565 26900 


	table(segmentation_model.gr@elementMetadata@listData$type)
    #LMR   MMR   SMR   TEM   UMR 
    #17630  4068  2492  2352 26353 	

	# We set mCG_cutoff_mean here to the value previously found, and reiterate mCHG cutoffs to be used later for classifying finished segments
	mCG_cutoff_mean = 0.09309094    # derived by fitting mixture of Gaussians to density distribution of segment mCG of whole parental methylome segmented in CG context
	m.sel = 0.09309094
	segment_mCHG_cutoff = 0.1262584 # derived by fitting mixture of Gaussians to density distribution of segment mCHG of whole parental methylome segmented in CHG context
	mCHG_cutoff_mean = 0.1262584
	mCHG_cutoff_median = 0.1262584
	segment_mCHH_cutoff = 0.1      # this is somewhat arbitrary, and is only used to classify the segments as UMR/MR in CHH.  May come back to this, and give it a better cutoff later, if this is thought to be useful
	mCHH_cutoff_mean = 0.1      
	

	# Plot mCG vs. segment size, with segments coloured by type 
	pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_genome_segment_mCG_vs_length_by_type_revised.pdf"))
	
	# Plot distribution of segment mCG and mCHG by type
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCG"))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=type), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=type), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.1)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCHG"))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHG"))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=values(segmentation_model.gr)@listData$median.meth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.2)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=data.frame(values(segmentation_model.gr)@listData)$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.2)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCHH"))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHH_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.2)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.2>mCHH_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.2)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHH"))
	
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.2>mCHH_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth.1, y=pmeth.2)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCHG") + ylab("Segment Mean mCHH"))

 
	### Cut-offs to separate UMR, GBM, TEM:
	# What threshold of mCG would capture 80% of TMR LMRs?
	#mCHG_quantile = 0.2
	#TMR_mCG_cutoff = quantile(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)) & (values(segmentation_model.gr)@listData$median.meth.1 > mCHG_cutoff) & (values(segmentation_model.gr)@listData$type=="LMR"),]$median.meth,mCHG_quantile)
	#cat(paste0(100*(1-mCHG_quantile),"% of 'transposon-like' Long Methylated Regions have mCG>",TMR_mCG_cutoff,".\n"))
	# Schmitz: 80% of 'transposon-like' Long Methylated Regions have mCG>0.754069511898385.

	# What threshold of mCG would capture 80% of GBM LMRs?
	#mCG_quantile = 0.8
	#GBM_mCG_cutoff = quantile(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)) & (values(segmentation_model.gr)@listData$median.meth.1 <= mCHG_cutoff) & (values(segmentation_model.gr)@listData$median.meth > m.sel) & (values(segmentation_model.gr)@listData$type=="LMR"),]$median.meth,mCG_quantile)
	#cat(paste0(100*mCG_quantile,"% of 'gene-body methylated-like' Long Methylated Regions have mCG<",GBM_mCG_cutoff,".\n"))
	# Schmitz: 80% of 'gene-body methylated-like' Long Methylated Regions have mCG<0.673452805081518.

	# Try to fit 3 Gaussians model to pick mCHG cutoff to separate TEMs from GBMs:
	#plot(normalmixEM(data.frame(values(segmentation_model.gr)@listData)[(data.frame(values(segmentation_model.gr)@listData)$pmeth>m.sel) & (data.frame(values(segmentation_model.gr)@listData)$pmeth.1>0.01) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)),]$pmeth.1,k=3), which=2, breaks=seq(0,1,0.01))
	
	#m.sel=0.2
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=mCHG), alpha=0.1) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>mCHG_cutoff*2,1,ifelse(median.meth>m.sel,2,3)), ncol=1) +geom_density2d() +scale_colour_gradient(low="green", high="red"))
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=type), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>mCHG_cutoff_median,paste0("1. 'Transposon-like' segments:\n median mCHG > ",round(mCHG_cutoff_median, digits=4)),ifelse(median.meth>m.sel,paste0("2. 'Gene-body-like' segments:\n median mCG > ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)),paste0("3. Unmethylated segments:\n median mCG <= ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1,mCHH=data.frame(values(segmentation_model.gr)@listData)$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=type), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	dev.off()
	
	
	# Crystalise segment methylation status as a column in the segmentation object
	# mCHG cutoff is segment_mCHG_cutoff
	# mCHH cutoff is segment_mCHH_cutoff
	# mCG cutoff is m.sel
	segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse((cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff) | (cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=data.frame(values(segmentation_model.gr)@listData)$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[,]$mCHH>segment_mCHH_cutoff),"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$pmeth>m.sel,"GBM","UMR"))

	# For Schmitz data, two of these segments have NA pmeth values
	#segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$mCHG>mCHG_cutoff_mean,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$pmeth>m.sel,"GBM","UMR"))

	# Set aside any unclassified segments
	length(segmentation_model.gr)
	# 52896  Schmitz: 55873, 55142
	
	#set_aside_no_class.gr = segmentation_model.gr[is.na(as.data.frame(values(segmentation_model.gr)@$segment.mStatus))]
	set_aside_no_class.gr = segmentation_model.gr[is.na(segmentation_model.gr@elementMetadata@listData$segment.mStatus)]
		
	segmentation_model.gr = segmentation_model.gr[!is.na(segmentation_model.gr@elementMetadata@listData$segment.mStatus)]
	length(segmentation_model.gr)
	# 52893  Schmitz: 55867, 55136

	cat(paste0("There are ",length(set_aside_no_class.gr), " segments with no methylation status classification.\n"))
	# 3  Schmitz data: 6

	# It's not exactly clear why this is not zero - if these segments don't overlap CHG or CHH sites they should have been already set aside.
	# On inspection, these segments all have mCG data (in Schmitz data set)
	# Check for overlaps with TEs or exons and characterise accordingly

	# Load TE and exon annotations so we can check for overlaps to resolve anomalies
 	annot_gff = read.delim(reference_gff, header=F, comment.char="#")	
	gff.exons = read.delim(reference_exons, header=F, comment.char="#")
	gff.transposons = annot_gff[annot_gff[,3]=="transposable_element",]
	gff.transposons$V1=mixedCaseChr(gff.transposons$V1)
	gff.exons$V1=mixedCaseChr(gff.exons$V1)
	gff.transposons$gene_ID=str_split_fixed(str_split_fixed(gff.transposons$V9, ';',3)[,1],"=",2)[,2]
	gff.exons$gene_ID=str_split_fixed(str_split_fixed(gff.exons$V9, ';',3)[,1],"=",2)[,2]
	gff.transposons=within(gff.transposons, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,3],"=",2)[,2]})
	transposon_ranges=makeGRangesFromDataFrame(df = gff.transposons, start.field = "V4", end.field = "V5",seqnames.field = "V1")
	exon_ranges=makeGRangesFromDataFrame(df = gff.exons, start.field = "V4", end.field = "V5",seqnames.field = "V1")

	# Set segments overlapping exons to "GBM" if mCG supports it
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))])@listData$pmeth>m.sel,"GBM","UMR")
	# Set segments overlapping transposons to "TEM" if mCG supports it (to override exons in case segment overlaps both)
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))])@listData$pmeth>m.sel,"TEM","UMR")
	# Presumption of "UMR" for segments which don't overlap exons or TEs
	set_aside_no_class.gr[is.na(set_aside_no_class.gr@elementMetadata@listData$segment.mStatus)]@elementMetadata@listData$segment.mStatus = "UMR"

	
	# Classify each of the other set_aside segments, where possible, and merge back into the model
	# Check first to make sure that none of the set_aside_segments.gr overlap with each other or with segments in the model (sense check)

	length(c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr))
	# 687  Schmitz: 1200
	length(reduce(c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr)))
	# 687  Schmitz: 1200
	length(findOverlaps(c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr), segmentation_model.gr))
	#0  Schmitz: 0
	
	# Split set_aside_no_*.gr into those which have data in other contexts and those that don't
	
	set_aside_no_CG_CHG.gr = setdiff(set_aside_no_CG.gr, set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, set_aside_no_CG.gr)))])
	length(set_aside_no_CG_CHG.gr)
	# 2  Schmitz: 35
	
	set_aside_no_CG_CHG_CHH.gr = setdiff(set_aside_no_CG_CHG.gr, set_aside_no_CG_CHG.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, set_aside_no_CG_CHG.gr)))])
	length(set_aside_no_CG_CHG_CHH.gr)
	# 0  Schmitz: 22  
	# These cannot be annotated - they both look like TEM by eye
	# (shouldnt have any effect as we found 0)
	set_aside_no_CG_CHG.gr = set_aside_no_CG_CHG.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, set_aside_no_CG_CHG.gr)))]
	length(set_aside_no_CG_CHG.gr)
	# 2  Schmitz: 33
	# These can be annotated using CHH
	
	set_aside_no_CG.gr = set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, set_aside_no_CG.gr)))]
	length(set_aside_no_CG.gr)
	# 3  Schmitz: 51
	
	set_aside_no_CG_CHH.gr = setdiff(set_aside_no_CG.gr, set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, set_aside_no_CG.gr)))])
	length(set_aside_no_CG_CHH.gr)
	# 1  Schmitz: 0
	set_aside_no_CG.gr = set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, set_aside_no_CG.gr)))]
	length(set_aside_no_CG.gr)
	# 2  Schmitz: 51
	# These can be annotated using CHG & CHH
	
	set_aside_no_CHG_CG.gr = setdiff(set_aside_no_CHG.gr, set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_CG.gr, set_aside_no_CHG.gr)))])
	length(set_aside_no_CHG_CG.gr)
	# 0  Schmitz: 0
	set_aside_no_CHG.gr = set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_CG.gr, set_aside_no_CHG.gr)))]
	length(set_aside_no_CHG.gr)
	# 661  Schmitz: 1099
	
	set_aside_no_CHG_CHH.gr = setdiff(set_aside_no_CHG.gr, set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, set_aside_no_CHG.gr)))])
	length(set_aside_no_CHG_CHH.gr)
	# 15  Schmitz: 12
	# These can be annotated using CG
	
	set_aside_no_CHG.gr = set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, set_aside_no_CHG.gr)))]
	length(set_aside_no_CHG.gr)
	# 646  Schmitz: 1087
	# These can be annotated using CG and CHH, although CHH is unreliable as a differentiator between GBM and TEM so also use annotation. Chances are though that any anomalies will be resolved by the merging of adjacent TEMs with short GBMs

	set_aside_no_CHH_CG.gr = setdiff(set_aside_no_CHH.gr, set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_CG.gr, set_aside_no_CHH.gr)))])
	length(set_aside_no_CHH_CG.gr)
	# 0  Schmitz: 0
	set_aside_no_CHH.gr = set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_CG.gr, set_aside_no_CHH.gr)))]
	length(set_aside_no_CHH.gr)
	# 21  Schmitz: 15
	
	set_aside_no_CHH_CHG.gr = setdiff(set_aside_no_CHH.gr, set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, set_aside_no_CHH.gr)))])
	length(set_aside_no_CHH_CHG.gr)
	# 0  Schmitz: 0
	set_aside_no_CHH.gr = set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, set_aside_no_CHH.gr)))]
	length(set_aside_no_CHH.gr)
	# 21   Schmitz: 15
	# These can be annotated using CG and CHG
	
	# Segments with no CG sites all overlap either CHH sites, or CHG and CHH sites
	# Segments with no CHG sites all overlap CG sites, and most overlap CHH sites
	# Segments with no CHH sites all overlap CG and CHG sites
	
	# Annotate set_aside_no_* with whatever from mCG, mCHG, mCHH is possible given the overlapping sites.
	# Use annotate with an available context in place of a missing one, to provide comparable data structure, but then zero out values.

	# For each of the sets of anomalous segments, set the methylation status as appropriate using annotation if necessary
	# Occasionally (2 cases in Schmitz data), classification fails even when segment overlaps some sites of a given class for unknown reason. In that case we presume "UMR" by default

	if (length(set_aside_no_CG.gr) > 0) {
		values(set_aside_no_CG.gr) = NA_by_segment(set_aside_no_CG.gr)
		values(set_aside_no_CG.gr) = cbind(values(set_aside_no_CG.gr),meth_by_segment(meth_CHG.gr, segment_model=set_aside_no_CG.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
		values(set_aside_no_CG.gr) = cbind(values(set_aside_no_CG.gr),meth_by_segment(meth_CHH.gr, segment_model=set_aside_no_CG.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
		
		# It is conceivable that a segment will only contain sites with zero coverage in one particular context, leading to pmeth=NaN
		# We try to pre-empt that by setting pmeth=0 for such segments
		set_aside_no_CG.gr@elementMetadata@listData[[12]][is.na(data.frame(set_aside_no_CG.gr@elementMetadata@listData[[12]]))]=0
		set_aside_no_CG.gr@elementMetadata@listData[[19]][is.na(data.frame(set_aside_no_CG.gr@elementMetadata@listData[[19]]))]=0
		
		# no CG sites so segments likely arise at the end of a TEM segment.  Use mCHG and mCHH to decide status
		set_aside_no_CG.gr@elementMetadata@listData$segment.mStatus = ifelse((data.frame(values(set_aside_no_CG.gr)@listData)$pmeth.1>segment_mCHG_cutoff) | (data.frame(values(set_aside_no_CG.gr)@listData)$pmeth.2>segment_mCHH_cutoff),"TEM","UMR")
	}
	
	if (length(set_aside_no_CG_CHG.gr) > 0) {
		values(set_aside_no_CG_CHG.gr) = cbind(NA_by_segment(set_aside_no_CG_CHG.gr), NA_by_segment(set_aside_no_CG_CHG.gr))
		values(set_aside_no_CG_CHG.gr) = cbind(values(set_aside_no_CG_CHG.gr), meth_by_segment(meth_CHH.gr, segment_model=set_aside_no_CG_CHG.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

		# It is conceivable that a segment will only contain sites with zero coverage in one particular context, leading to pmeth=NaN
		# We try to pre-empt that by setting pmeth=0 for such segments
		set_aside_no_CG_CHG.gr@elementMetadata@listData[[19]][is.na(data.frame(set_aside_no_CG_CHG.gr@elementMetadata@listData[[19]]))]=0

		# For this one no CG so it arises at the end of a TEM segment, most likely.  Use mCHH to decide (though may be safer to refer to annotation, as CHH in absence of CHG may be unreliable)
		set_aside_no_CG_CHG.gr@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_CG_CHG.gr)@listData)$pmeth.2>segment_mCHH_cutoff,"TEM","UMR")
	}
	
	if (length(set_aside_no_CG_CHH.gr) > 0) {
		values(set_aside_no_CG_CHH.gr) = cbind(NA_by_segment(set_aside_no_CG_CHH.gr), meth_by_segment(meth_CHG.gr, segment_model=set_aside_no_CG_CHH.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")), NA_by_segment(set_aside_no_CG_CHH.gr))

		# It is conceivable that a segment will only contain sites with zero coverage in one particular context, leading to pmeth=NaN
		# We try to pre-empt that by setting pmeth=0 for such segments
		set_aside_no_CG_CHH.gr@elementMetadata@listData[[12]][is.na(data.frame(set_aside_no_CG_CHH.gr@elementMetadata@listData[[12]]))]=0

		# For this one no CG so it arises at the end of a TEM segment, most likely.  Use mCHG to decide
		set_aside_no_CG_CHH.gr@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_CG_CHH.gr)@listData)$pmeth.1>segment_mCHG_cutoff,"TEM","UMR")
	}
	
	if (length(set_aside_no_CHG_CHH.gr) > 0) {
		values(set_aside_no_CHG_CHH.gr) = meth_by_segment(meth_CG.gr, segment_model=set_aside_no_CHG_CHH.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))
		values(set_aside_no_CHG_CHH.gr) = cbind(values(set_aside_no_CHG_CHH.gr), NA_by_segment(set_aside_no_CHG_CHH.gr), NA_by_segment(set_aside_no_CHG_CHH.gr))

		# It is conceivable that a segment will only contain sites with zero coverage in one particular context, leading to pmeth=NaN
		# We try to pre-empt that by setting pmeth=0 for such segments
		set_aside_no_CHG_CHH.gr@elementMetadata@listData[[5]][is.na(data.frame(set_aside_no_CHG_CHH.gr@elementMetadata@listData[[5]]))]=0

		# Check overlaps with TE or exons to decide whether TEM or GBM if mCG>m.sel
		#set_aside_exon_hits = unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG_CHH.gr)))
		#set_aside_transposon_hits = unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG_CHH.gr)))

		set_aside_no_CHG_CHH.gr@elementMetadata@listData$segment.mStatus = as.character(rep(NA, length(set_aside_no_CHG_CHH.gr)))
		# Set all those overlapping exons to "GBM"
		set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG_CHH.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG_CHH.gr)))])@listData)$pmeth>m.sel,"GBM","UMR")
		# Set all those overlapping transposons to "TEM" (to override exons in case overlaps both)
		set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG_CHH.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG_CHH.gr)))])@listData)$pmeth>m.sel,"TEM","UMR")
		# For those overlapping neither exons nor TEs, start with a presumption of "UMR" (why?)
		set_aside_no_CHG_CHH.gr[is.na(set_aside_no_CHG_CHH.gr@elementMetadata@listData$segment.mStatus)]@elementMetadata@listData$segment.mStatus = rep("UMR", length(set_aside_no_CHG_CHH.gr[is.na(set_aside_no_CHG_CHH.gr@elementMetadata@listData$segment.mStatus)]))
	}

	if (length(set_aside_no_CHG.gr) > 0) {
		values(set_aside_no_CHG.gr) = meth_by_segment(meth_CG.gr, segment_model=set_aside_no_CHG.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))
		values(set_aside_no_CHG.gr) = cbind(values(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr))
		values(set_aside_no_CHG.gr) = cbind(values(set_aside_no_CHG.gr),meth_by_segment(meth_CHH.gr, segment_model=set_aside_no_CHG.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

		# It is conceivable that a segment will only contain sites with zero coverage in one particular context, leading to pmeth=NaN
		# We try to pre-empt that by setting pmeth=0 for such segments
		set_aside_no_CHG.gr@elementMetadata@listData[[5]][is.na(data.frame(set_aside_no_CHG.gr@elementMetadata@listData[[5]]))]=0
		set_aside_no_CHG.gr@elementMetadata@listData[[19]][is.na(data.frame(set_aside_no_CHG.gr@elementMetadata@listData[[19]]))]=0

		# For this one no CHG so need to check CG and CHH
		set_aside_no_CHG.gr@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_CHG.gr)@listData)$pmeth.2>segment_mCHH_cutoff,"TEM",ifelse(data.frame(values(set_aside_no_CHG.gr)@listData)$pmeth>m.sel,"GBM","UMR"))
		
		# In 2 cases (Schmitz) the resulting classification failed with too little coverage in CHH
		# Set all those overlapping exons to "GBM"
		set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)][unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)])))]$segment.mStatus = ifelse(as.data.frame(values(set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)][unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)])))])@listData)$pmeth>m.sel,"GBM","UMR") 
		# Set all those overlapping transposons to "TEM" (to override exons in case overlaps both)
		set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)][unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)])))]$segment.mStatus = ifelse(as.data.frame(values(set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)][unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(values(set_aside_no_CHG.gr)@listData)$segment.mStatus)])))])@listData)$pmeth>m.sel,"GBM","UMR") 

		# For those overlapping neither exons nor TEs, start with a presumption of "UMR" (why?)
		set_aside_no_CHG.gr@elementMetadata@listData$segment.mStatus[is.na(set_aside_no_CHG.gr@elementMetadata@listData$segment.mStatus)] = "UMR"
	}
	
	if (length(set_aside_no_CHH.gr) > 0) {
		values(set_aside_no_CHH.gr) = meth_by_segment(meth_CG.gr, segment_model=set_aside_no_CHH.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))
		values(set_aside_no_CHH.gr) = cbind(values(set_aside_no_CHH.gr),meth_by_segment(meth_CHG.gr, segment_model=set_aside_no_CHH.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
		values(set_aside_no_CHH.gr) = cbind(values(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr))

		# It is conceivable that a segment will only contain sites with zero coverage in one particular context, leading to pmeth=NaN
		# We try to pre-empt that by setting pmeth=0 for such segments
		set_aside_no_CHH.gr@elementMetadata@listData[[5]][is.na(data.frame(set_aside_no_CHH.gr@elementMetadata@listData[[5]]))]=0
		set_aside_no_CHH.gr@elementMetadata@listData[[12]][is.na(data.frame(set_aside_no_CHH.gr@elementMetadata@listData[[12]]))]=0

		# For this one no CHH so need to check CG and CHG and whether overlaps exon or TE
		set_aside_no_CHH.gr@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_CHH.gr)@listData)$pmeth.1>segment_mCHG_cutoff,"TEM",ifelse(data.frame(values(set_aside_no_CHH.gr)@listData)$pmeth>m.sel,"GBM","UMR"))
	}

	if (length(set_aside_no_CG_CHG_CHH.gr) > 0) {
		# For the segments with no overlapping sites, assign some spurious values from one of the other sets to keep the data structure correct
		values(set_aside_no_CG_CHG_CHH.gr) = cbind(NA_by_segment(set_aside_no_CG_CHG_CHH.gr),NA_by_segment(set_aside_no_CG_CHG_CHH.gr),NA_by_segment(set_aside_no_CG_CHG_CHH.gr)) 
		#values(set_aside_no_CHG.gr)[1:length(set_aside_no_CG_CHG_CHH.gr),]

		# Set their methylation status to "TEM" based on visual inspection of the 2 cases in Schmitz data
		set_aside_no_CG_CHG_CHH.gr@elementMetadata@listData$segment.mStatus = as.character(rep("TEM", length(set_aside_no_CG_CHG_CHH.gr)))
	}
	
	set_aside_segments.gr = c(set_aside_no_class.gr, set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr, set_aside_no_CG_CHG.gr, set_aside_no_CG_CHH.gr, set_aside_no_CHG_CHH.gr, set_aside_no_CG_CHG_CHH.gr)
	# Check that set_aside segments are all classified now
	#length(c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr, set_aside_no_CG_CHG.gr, set_aside_no_CG_CHH.gr, set_aside_no_CHG_CHH.gr, set_aside_no_CG_CHG_CHH.gr))
	sum(is.na(set_aside_segments.gr@elementMetadata@listData$segment.mStatus))
	# 0
	
	# Now that all set aside segments have been classified, merge all the set aside segments back into the model
	segmentation_model.gr = sort(sortSeqlevels(c(segmentation_model.gr, set_aside_segments.gr)))

	# Use the reduce function on each segment type, separately, to concatenate any adjoining segments of the same type
	segmentation_model.gr = c(reduce(segmentation_model.gr[segmentation_model.gr@elementMetadata@listData$segment.mStatus == "TEM"]), reduce(segmentation_model.gr[segmentation_model.gr@elementMetadata@listData$segment.mStatus == "GBM"]), reduce(segmentation_model.gr[segmentation_model.gr@elementMetadata@listData$segment.mStatus == "UMR"]))
	length(segmentation_model.gr)
	# 51759 (should be more - missing some set_asides) Schmitz: 54028, 54341
	
	# Adjust segments to set aside any segments which don't overlap with any CG sites. We will deal with them later
	# Counts here are first masking only CHRM and CHRC
	set_aside_no_CG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 64
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 54301
	set_aside_no_CHG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CHG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 1048
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 53253
	set_aside_no_CHH.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CHH.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 15
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 53238


	# Reannotate the segmentation model (again!)
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_CG.gr, segment_model=segmentation_model.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_CHG.gr, segment_model=segmentation_model.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHG_cutoff), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_CHH.gr, segment_model=segmentation_model.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHH_cutoff), mMeth.classes=c("UMR","MR")))
	segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$pmeth>m.sel,"GBM","UMR"))

	# Check whether any unannotated segments remain in the model
	set_aside_no_class.gr = segmentation_model.gr[is.na(data.frame(values(segmentation_model.gr)@listData)$segment.mStatus)]
	length(set_aside_no_class.gr)
	# 4
	segmentation_model.gr = segmentation_model.gr[!is.na(data.frame(values(segmentation_model.gr)@listData)$segment.mStatus)]
	length(segmentation_model.gr)
	# 53234

	
	# Add back in any of the original set aside segments which don't overlap with the new model
	segmentation_model.gr = c(segmentation_model.gr, set_aside_segments.gr[!(set_aside_segments.gr %in% set_aside_segments.gr[unique(subjectHits(findOverlaps(segmentation_model.gr, set_aside_segments.gr)))])])
	length(segmentation_model.gr)
	# 54341
	
	
	
	# We have a small number of remaining 'GBM' segments which are adjacent to TEM segments, rather than being separated by UMR.  These should be merged with their neighbour.
	# We do not want to merge longer GBM segments, as they are likely real GBM, so we set a nt length cutoff, below which, any GBM segments will be merged with their neighboring TEM segment if they share a boundary:
	
	GBM_distinct_min_width=500
	
	# rebuild the segmentation model by:
	#  concatenate TEM segments with GBM segments whose segment width is below the threshold
	#  reduce the above segments to merge direct neighbours into a single segment
	#  concatenate the above merged set with the original segmentation model, and reduce this so that the new longer segments subsume the pairs of original segments
	merged_segments.gr = reduce(c(segmentation_model.gr[data.frame(values(segmentation_model.gr)@listData)$segment.mStatus == "TEM"], segmentation_model.gr[(data.frame(values(segmentation_model.gr)@listData)$segment.mStatus == "GBM") & (segmentation_model.gr@ranges@width < GBM_distinct_min_width)]))
	length(merged_segments.gr)
	# 16195, 16527

	# remove tha annotations from the segmentation model remainder, so it will merge properly with the new, unannotated, merged segments
	segmentation_model_remainder.gr = segmentation_model.gr[!((segmentation_model.gr$segment.mStatus == "TEM") | ((segmentation_model.gr$segment.mStatus == "GBM") & (segmentation_model.gr@ranges@width < GBM_distinct_min_width)))]
	length(segmentation_model_remainder.gr)
	# 36042, 35382
	segmentation_model_remainder.gr@elementMetadata@listData = list()

	#str(merged_segments.gr)
	#str(segmentation_model_remainder.gr)
		
	#Annotate the two subsets of the model with some spurious annotations - lack of anything in elementMetadata can cause sort to fail
	values(merged_segments.gr) = NA_by_segment(merged_segments.gr)
	values(segmentation_model_remainder.gr) = NA_by_segment(segmentation_model_remainder.gr)

	# put the whole model back together again by concatenating the new merged segments with all the other segments from the original model which were not considered for merging
	# May as well sort the new model as we build it, in case we need a sorted model later
	segmentation_model.gr = sort(sortSeqlevels(c(merged_segments.gr, segmentation_model_remainder.gr)))
	length(segmentation_model.gr)
	# 52237, 51909

	#Remove the spurious annotations
    segmentation_model.gr@elementMetadata@listData = list()

	
	# Adjust segments to set aside any segments which don't overlap with any CG sites. We will deal with them later
	# Counts here are first masking only CHRM and CHRC
	set_aside_no_CG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 40
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 51869
	set_aside_no_CHG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 504
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 51365
	set_aside_no_CHH.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 11
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 51354


	
	# Reannotate the segmentation model (again!)
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_CG.gr, segment_model=segmentation_model.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_CHG.gr, segment_model=segmentation_model.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHG_cutoff), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_CHH.gr, segment_model=segmentation_model.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHH_cutoff), mMeth.classes=c("UMR","MR")))
	#segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$pmeth>m.sel,"GBM","UMR"))
	segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse((data.frame(values(segmentation_model.gr)@listData)$pmeth.1>segment_mCHG_cutoff) | (data.frame(values(segmentation_model.gr)@listData)$pmeth.2>segment_mCHH_cutoff),"TEM",ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth>m.sel,"GBM","UMR"))

	# Add back in any of the original set aside segments which don't overlap with the new model, and sort for good measure
	segmentation_model.gr = sort(sortSeqlevels(c(segmentation_model.gr, set_aside_segments.gr[!(set_aside_segments.gr %in% set_aside_segments.gr[unique(subjectHits(findOverlaps(segmentation_model.gr, set_aside_segments.gr)))])])))
	length(segmentation_model.gr)
	# 51909

    # Rebuild set_aside_no_class.gr in case any new segments have arisen
	set_aside_no_class.gr = c(set_aside_no_class.gr, segmentation_model.gr[is.na(data.frame(values(segmentation_model.gr)@listData)$segment.mStatus)])
	length(set_aside_no_class.gr)
	# 1
	segmentation_model.gr = segmentation_model.gr[!is.na(data.frame(values(segmentation_model.gr)@listData)$segment.mStatus)]
	length(segmentation_model.gr)
	# 53234

	# Assign status to the set_aside_no_class.gr segments
	# Set segments overlapping exons to "GBM" if mCG supports it
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))])@listData)$pmeth>m.sel,"GBM","UMR")
	# Set segments overlapping transposons to "TEM" if mCG supports it (to override exons in case segment overlaps both)
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(data.frame(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))])@listData)$pmeth>m.sel,"TEM","UMR")
	# Presumption of "UMR" for segments which don't overlap exons or TEs
	set_aside_no_class.gr[is.na(as.data.frame(values(set_aside_no_class.gr)@listData)$segment.mStatus)]@elementMetadata@listData$segment.mStatus = "UMR"

    # Merge the set_aside_no_class.gr segments back into the model
	segmentation_model.gr = c(segmentation_model.gr, set_aside_no_class.gr)
	
	# Output a text version of the segmentation model for visualisation
	write.table(cbind("Chromosome"=data.frame(segmentation_model.gr@seqnames),"Start"=data.frame(segmentation_model.gr@ranges)$start,"End"=data.frame(segmentation_model.gr@ranges)$end,"Type"=data.frame(values(segmentation_model.gr)@listData)$segment.mStatus),file=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_draft.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Save a copy in case we want to come back to it later
	saveUMRLMRSegments(segs=segmentation_model.gr, GRangesFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_draft.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_draft.tsv"))


	table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$segment.mStatus)
	# These are the counts of segments by type, based on coverage in all 3 contexts.  There are additional segments classified on the basis of two, or one contexts, plus annotation, which are included in the second table below.
	
	#Schmitz data:
	#GBM   TEM   UMR 
	#19738  7500 26313 
	#18292  7460 25609 (segmentation_model_draft files)
	
	# After annotation-driven reannotation of GBM segments to TEM if they don't overlap a gene (segmentation_model_draft2 files):
	# GBM   TEM   UMR 
	#16629  8793 25607 

	table(as.data.frame(values(segmentation_model.gr)@listData$segment.mStatus))
	#18774  7484 25651  (segmentation_model_draft files)
	#17006  8919 25648  (segmentation_model_draft2 files)
	
	
	# Plot mCG vs. segment size, with segments coloured by segment methylation status
	pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_genome_segment_mCG_vs_length_by_mStatus.pdf"))

	# Plot distribution of segment mCG, mCHG, mCHH by status
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Median mCG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Mean mCG"))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.1)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.1)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Mean mCHG"))

	# Plot segment mCHG vs mCG
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHG"))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=values(segmentation_model.gr)@listData$median.meth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.2)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=data.frame(values(segmentation_model.gr)@listData)$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.2)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Mean mCHH"))

	# Plot segment mCHH vs mCG
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$median.meth.2>mCHH_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.2>segment_mCHH_cutoff,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHH"))
	
	# Plot segment mCHH vs mCHG
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$median.meth.2>mCHH_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth.1, y=median.meth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCHG") + ylab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.2>segment_mCHH_cutoff,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth.1, y=pmeth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCHG") + ylab("Segment Mean mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(data.frame(values(segmentation_model.gr)@listData)$pmeth.2>segment_mCHH_cutoff,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth.1, y=pmeth.2)) + geom_point(aes(colour=pmeth, size=(seg_width)/1000), alpha=0.5) + theme_minimal() + xlab("Segment Mean mCHG") + ylab("Segment Mean mCHH") + scale_colour_gradient("Segment\nMean mCG",low="green", high="red"))
	
	# Plot m* vs length by mStatus
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>mCHG_cutoff_median,paste0("1. 'Transposon-like' segments:\n median mCHG > ",round(mCHG_cutoff_median, digits=4)),ifelse(median.meth>m.sel,paste0("2. 'Gene-body-like' segments:\n median mCG > ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)),paste0("3. Unmethylated segments:\n median mCG <= ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>segment_mCHG_cutoff,paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>segment_mCHG_cutoff | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1,mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHH>segment_mCHH_cutoff,paste0("1. 'Transposon-like' segments:\n mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))

	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1,mCHH=data.frame(values(segmentation_model.gr)@listData)$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))

	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1,mCHH=data.frame(values(segmentation_model.gr)@listData)$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)),], aes(x=seg_width, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1,mCHH=data.frame(values(segmentation_model.gr)@listData)$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)),], aes(x=seg_width, y=pmeth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHH") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))


	
	dev.off()

	
	# Now the methylome is segmented irrespective of annotation.
    # Start using the annotation to refine the segmentation 
	
	# Redefine genomicranges representing the whole genome
	genome_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	levels(genome_loci$Chromosome)=c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrM","ChrC")
	for (chrom_no in 1:length(sLengths)) {
		genome_loci[chrom_no,]=list(Chromosome=names(sLengths[chrom_no]),start_site=1,end_site=sLengths[chrom_no])
	}
	genome_loci.gr=makeGRangesFromDataFrame(df=genome_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")

	gff.genes$Chromosome = mixedCaseChr(gff.genes$V1)
	gff.transposons$Chromosome = mixedCaseChr(gff.transposons$V1)

    # Mask out the part of Chr2 containing mitochondrial genes
	segment_mask_loci.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="Chr2", "start_site"=3239693, "end_site"=3505260))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")
			
	masked_segmentation_model.gr = segmentation_model.gr[unique(queryHits(findOverlaps(segmentation_model.gr, setdiff(segmentation_model.gr,segment_mask_loci.gr))))]

	
	#GBM_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	#Het_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Heterochromatic",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	#UM_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	#other_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.genes[!(gff.genes$m_class %in% c("Gene-body Methylated","Heterochromatic","Unmethylated")),c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	#GBM_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[gff.transposons$m_class=="Gene-body Methylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	#Het_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[gff.transposons$m_class=="Heterochromatic",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	#UM_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[gff.transposons$m_class=="Unmethylated",c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	#other_TE_loci.gr = reduce(makeGRangesFromDataFrame(df=gff.transposons[!(gff.transposons$m_class %in% c("Gene-body Methylated","Heterochromatic","Unmethylated")),c("Chromosome","V4","V5")], start.field="V4", end.field="V5", seqnames.field="Chromosome"))
	
	#Unannotated_loci.gr = setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr, GBM_gene_loci.gr, Het_gene_loci.gr, UM_gene_loci.gr, other_gene_loci.gr, GBM_TE_loci.gr, Het_TE_loci.gr, UM_TE_loci.gr, other_TE_loci.gr))
	Unannotated_loci.gr = setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr, gene_ranges, transposon_ranges))
	
	segmentation_model.gr = masked_segmentation_model.gr

#7316

#8719

	# At this point, segmentation_model.gr == masked_segmentation_model.gr
	
	# We will make amendments to segmentation_model.gr which we could potentially rollback by setting segmentation_model.gr=masked_segmentation_model.gr
	
	# Capture initial state of segments from model, so that x_x_x_hits will continue to work properly
	GBM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]
	TEM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="TEM"]
	UMR_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="UMR"]
	length(GBM_segments.gr)
	length(TEM_segments.gr)
	length(UMR_segments.gr)
	# 16992, Schmitz: 18730
	# 4530, Schmitz: 7447
	# 21211, Schmitz: 25648
    	
	G_G_hits=findOverlaps(gene_ranges, segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"])
	G_T_hits=findOverlaps(transposon_ranges, segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"])

	# GBM segments overlapping an annotated TE:
	G_x_T_segments = unique(subjectHits(G_T_hits))
	length(G_x_T_segments)
	# 1185, Schmitz: 1525 segments from the ones currently annotated as GBM

	# GBM segments overlapping an annotated TE, but no annotated genes:
	# (find which of the G_x_T segments overlaps with any of the gene loci ranges, and take the set difference to exclude these from the G_x_T segments)
	G_x_T_no_gene_segments.gr = setdiff(GBM_segments.gr[G_x_T_segments], GBM_segments.gr[G_x_T_segments][unique(subjectHits(findOverlaps(gene_ranges,GBM_segments.gr[G_x_T_segments])))])
	length(G_x_T_no_gene_segments.gr)
	# 437, Schmitz: 1014 segments
	
	# Set such segments to be GBM-like rather than GBM
	segmentation_model.gr[unique(subjectHits(findOverlaps(G_x_T_no_gene_segments.gr, segmentation_model.gr)))]$segment.mStatus = rep("GBM-like", length(unique(subjectHits(findOverlaps(G_x_T_no_gene_segments.gr, segmentation_model.gr)))))
	
	
	# GBM segments overlapping an annotated gene or TE:
	G_x_x_segments = unique(c(G_x_T_segments, subjectHits(G_G_hits)))
	length(G_x_x_segments)
	# 16007, Schmitz: 18015 segments from the ones currently annotated as GBM

	# GBM segments not overlapping an annotated TE, nor an annotated gene:
	# (find which non-G_x_G segments overlaps with any of the gene loci ranges, and take the set difference to exclude these from the G_x_T segments)
	G_x_x_no_annotation_segments.gr = setdiff(GBM_segments.gr, GBM_segments.gr[G_x_x_segments])
	length(G_x_x_no_annotation_segments.gr)
	# 984, Schmitz: 714 segments
	
	# Set such segments to be GBM-like rather than GBM
	segmentation_model.gr[unique(subjectHits(findOverlaps(G_x_x_no_annotation_segments.gr, segmentation_model.gr)))]$segment.mStatus = rep("GBM-like", length(unique(subjectHits(findOverlaps(G_x_x_no_annotation_segments.gr, segmentation_model.gr)))))
	
	# Recapture modified state of segments from model, so that segmentation model can be rebuilt after reduce operation
	GBM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]
	TEM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="TEM"]
	UMR_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="UMR"]
	GBM_like_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM-like"]
	length(GBM_segments.gr)
	length(TEM_segments.gr)
	length(UMR_segments.gr)
	length(GBM_like_segments.gr)
	# 15570, Schmitz: 17001
	# 4530, Schmitz: 7447
	# 21211, Schmitz: 25648
	# 1422, Schmitz: 1729

	# We still have some GBM segments which overflow past the bounds of the annotated genes.  We will trim these to end at the end of the gene body, and reannotate accordingly.  We will also annotate the remaining segments which lie outside gene bodies.
	
	G_x_G_hits = findOverlaps(reduce(gene_ranges), reduce(segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]))
	# Define set of segments representing the overlap between GBM segments and gene models
	true_GBM_segments.gr = setdiff(reduce(segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"])[subjectHits(G_x_G_hits)],setdiff(reduce(segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"])[subjectHits(G_x_G_hits)],reduce(gene_ranges)[queryHits(G_x_G_hits)]))
	length(true_GBM_segments.gr)
	# 15704, Schmitz: 17285
	# Define set of segments representing the parts of GBM segments which are external of gene models
	ext_GBM_segments.gr = setdiff(reduce(segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]), true_GBM_segments.gr) 
	length(ext_GBM_segments.gr)
	# 2306, Schmitz: 892
	
	# Read the parental methylome data in as Granges objects, to use to annotate new segments
	#meth_CG.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	#meth_CHG.gr <- readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	#meth_CHH.gr <- readMethylome(FileName=paste0(project_id,"_CHH_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)

	# First the 'real_GBM' segments
	# Set aside any which overlap no CG sites
	set_aside_no_CG.gr = setdiff(true_GBM_segments.gr, true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, true_GBM_segments.gr)))])
	cat(paste0("true_GBM_segments.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 116, Schmitz: 70
	true_GBM_segments.gr = true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, true_GBM_segments.gr)))]
	cat(paste0("true_GBM_segments.gr contains ",length(true_GBM_segments.gr)," segments after trimming segments with no CG sites.\n"))
	# 15588, Schmitz: 17215
	# Annotate the remainder with their mean mCG values
	values(true_GBM_segments.gr) = cbind(values(true_GBM_segments.gr),meth_by_segment(meth_CG.gr, segment_model=true_GBM_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any of the remainder which overlap no CHG sites 
	set_aside_no_CHG.gr = setdiff(true_GBM_segments.gr, true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, true_GBM_segments.gr)))])
	cat(paste0("true_GBM_segments.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 553, Schmitz: 389
	true_GBM_segments.gr = true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, true_GBM_segments.gr)))]
	cat(paste0("true_GBM_segments.gr contains ",length(true_GBM_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 15035, Schmitz: 16826
	# Annotate the remainder with their mean mCHG values
	values(true_GBM_segments.gr) = cbind(values(true_GBM_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=true_GBM_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any of the remainder which overlap no CHH sites 
	set_aside_no_CHH.gr = setdiff(true_GBM_segments.gr, true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, true_GBM_segments.gr)))])
	cat(paste0("true_GBM_segments.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 71, sCHMITZ: 10
	true_GBM_segments.gr = true_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, true_GBM_segments.gr)))]
	cat(paste0("true_GBM_segments.gr contains ",length(true_GBM_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 14964, sCHMITZ: 16816
	# Annotate the remainder with their mean mCHH values
	values(true_GBM_segments.gr) = cbind(values(true_GBM_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=true_GBM_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	true_GBM_segments.gr$segment.mStatus = rep("GBM", length(true_GBM_segments.gr))
			
	# Function to return NAs for cases where segmentation model doen's overlap any sites in the given context
	NA_by_segment <- function (segment_model)
	{
		segment_count = length(segment_model)
		nSites.segmentation = rep(as.integer(NA),segment_count)
		nSites = rep(as.integer(NA),segment_count)
		T = array(NA,segment_count)
		M = array(NA,segment_count)
		pmeth = array(NA,segment_count)
		median.meth = rep(as.numeric(NA),segment_count)
		type = array(NA,segment_count)

		DataFrame(nSites.segmentation, nSites, T, M, pmeth = M/T, median.meth = median.meth, type)
	}
	
	# Define a corresponding convenience function to save an annotated genomic ranges object
	# Function adapted from MethylSeekR package
	saveUMRLMRSegments <- function (segs, GRangesFilename = NULL, TableFilename = NULL) 
	{
		# Replaced nCG with nSites for generality
		if (!is.null(GRangesFilename)) 
			saveRDS(segs, GRangesFilename)
		if (!is.null(TableFilename)) {
			D = data.frame(chr = as.character(seqnames(segs)), start = start(segs), 
				end = end(segs), type = values(segs)$type, nSites.segmentation = values(segs)$nSites.segmentation, 
				nSites.seq = values(segs)$nSites, mean.meth = 100 * values(segs)$pmeth, 
				median.meth = 100 * values(segs)$median.meth)
			write.table(D, file = TableFilename, quote = FALSE, sep = "\t", 
			row.names = FALSE)
		}
	}

	# For now, we will annotate mean methylation values of the set-aside segments as NA, and annotate them all as GBM
	values(set_aside_no_CG.gr) = cbind(NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr))
	values(set_aside_no_CHG.gr) = cbind(NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr))
	values(set_aside_no_CHH.gr) = cbind(NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr))
	set_aside_segments.gr = c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr)
	set_aside_segments.gr$segment.mStatus = rep("GBM",length(set_aside_segments.gr))
	true_GBM_segments.gr = c(true_GBM_segments.gr, set_aside_segments.gr)

	
	# Now tackle the GBM_leftovers.gr segments, which are not true GBM, but may not be GBM-like either
	# Set aside any which overlap no CG sites
	set_aside_no_CG.gr = setdiff(ext_GBM_segments.gr, ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, ext_GBM_segments.gr)))])
	cat(paste0("ext_GBM_segments.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 367, sCHMITZ: 242
	ext_GBM_segments.gr = ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CG.gr, ext_GBM_segments.gr)))]
	cat(paste0("ext_GBM_segments.gr contains ",length(ext_GBM_segments.gr)," segments after trimming segments with no CG sites.\n"))
	# 1939, sCHMITZ: 650
	# Annotate the remainder with their mean mCG values
	values(ext_GBM_segments.gr) = cbind(values(ext_GBM_segments.gr),meth_by_segment(meth_CG.gr, segment_model=ext_GBM_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any which overlap no CHG sites
	set_aside_no_CHG.gr = setdiff(ext_GBM_segments.gr, ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, ext_GBM_segments.gr)))])
	cat(paste0("ext_GBM_segments.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 201, Schmitz: 47
	ext_GBM_segments.gr = ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHG.gr, ext_GBM_segments.gr)))]
	cat(paste0("ext_GBM_segments.gr contains ",length(ext_GBM_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 1738, Schmitz: 603
	# Annotate the remainder with their mean mCHG values
	values(ext_GBM_segments.gr) = cbind(values(ext_GBM_segments.gr),meth_by_segment(meth_CHG.gr, segment_model=ext_GBM_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set aside any which overlap no CHH sites
	set_aside_no_CHH.gr = setdiff(ext_GBM_segments.gr, ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, ext_GBM_segments.gr)))])
	cat(paste0("ext_GBM_segments.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 9, Schmitz: 0
	ext_GBM_segments.gr = ext_GBM_segments.gr[unique(subjectHits(findOverlaps(meth_CHH.gr, ext_GBM_segments.gr)))]
	cat(paste0("ext_GBM_segments.gr contains ",length(ext_GBM_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 1729, Schmitz: 603
	# Annotate the remainder with their mean mCHH values
	values(ext_GBM_segments.gr) = cbind(values(ext_GBM_segments.gr),meth_by_segment(meth_CHH.gr, segment_model=ext_GBM_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))

	# Set the segment class of the ext_GBM_segments.gr based on their methylation means
	ext_GBM_segments.gr@elementMetadata@listData$segment.mStatus = ifelse((cbind(data.frame(values(ext_GBM_segments.gr)@listData),mCHG=data.frame(values(ext_GBM_segments.gr)@listData)$pmeth.1, seg_width=ext_GBM_segments.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff) | (cbind(data.frame(values(ext_GBM_segments.gr)@listData),mCHH=data.frame(values(ext_GBM_segments.gr)@listData)$pmeth.2, seg_width=ext_GBM_segments.gr@ranges@width)[,]$mCHH>segment_mCHH_cutoff),"TEM",ifelse(cbind(data.frame(values(ext_GBM_segments.gr)@listData),mCHG=data.frame(values(ext_GBM_segments.gr)@listData)$pmeth.1, seg_width=ext_GBM_segments.gr@ranges@width)[,]$pmeth>m.sel,"GBM-like","UMR"))
	### We should then reduce each class of the model, since we have introduced a number of new segments with annotations that already exist in the model, and may be adjacent

	# For now, we will annotate mean methylation values of the set-aside segments as NA, and annotate them all as GBM-ext
	### A better strategy would be to allow them to take on the annotation of their neighbour, then merge them, then reannotate the neighbours
	values(set_aside_no_CG.gr) = cbind(NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr), NA_by_segment(set_aside_no_CG.gr))
	values(set_aside_no_CHG.gr) = cbind(NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr))
	values(set_aside_no_CHH.gr) = cbind(NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr))
	set_aside_segments.gr = c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr)
	
	#set_aside_segments.gr$segment.mStatus = rep("GBM-ext",length(set_aside_segments.gr))
	### Setting segment.mStatus to "UMR" is almost always the right thing to do here since a genuine GBM segment within a gene body almost invariably results in the gene having a UMR region surrounding it. In fact these set-aside segments are often an indicator that there is a larger unmethylated segment that should expand from the seed of the set-aside segment
	set_aside_segments.gr$segment.mStatus = rep("UMR",length(set_aside_segments.gr))
	
	### We would then ideally need to reduce the UMR segments in the model, and reannotate them, but the numbers are so small as to be inconsequential
	ext_GBM_segments.gr = c(ext_GBM_segments.gr, set_aside_segments.gr)
	
	# We may have now left some GBM segments without any CG sites, or having had their mCG sites shipped out to an ext_GBM_segments.gr segment.
	# Reannotate all GBM sites 
	
	GBM_like_segments.gr$ID=NULL
	GBM_like_segments.gr$CG_site_count=NULL
	GBM_like_segments.gr$variant_count=NULL
	GBM_like_segments.gr$all_M_count=NULL
	GBM_like_segments.gr$all_U_count=NULL
	GBM_like_segments.gr$m_to_u_count=NULL
	GBM_like_segments.gr$u_to_m_count=NULL
	TEM_segments.gr$ID=NULL
	TEM_segments.gr$CG_site_count=NULL
	TEM_segments.gr$variant_count=NULL
	TEM_segments.gr$all_M_count=NULL
	TEM_segments.gr$all_U_count=NULL
	TEM_segments.gr$m_to_u_count=NULL
	TEM_segments.gr$u_to_m_count=NULL
	UMR_segments.gr$ID=NULL
	UMR_segments.gr$CG_site_count=NULL
	UMR_segments.gr$variant_count=NULL
	UMR_segments.gr$all_M_count=NULL
	UMR_segments.gr$all_U_count=NULL
	UMR_segments.gr$m_to_u_count=NULL
	UMR_segments.gr$u_to_m_count=NULL
	
	
	segmentation_model.gr = sort(sortSeqlevels(c(true_GBM_segments.gr, ext_GBM_segments.gr, GBM_like_segments.gr, TEM_segments.gr, UMR_segments.gr)))
	length(segmentation_model.gr)
	# 45173, Schmitz: 53001

	#write.table(cbind("Chromosome"=as.character(as.data.frame(segmentation_model.gr)$seqnames),"Start"=as.data.frame(segmentation_model.gr)$start,"End"=as.data.frame(segmentation_model.gr)$end,"Type"=as.data.frame(values(segmentation_model.gr)@listData)$segment.mStatus),file=paste0(project_id,"_",meth_context,"_segmentation_model_draft2.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
    write.table(cbind("Chromosome"=as.character(segmentation_model.gr$seqnames),"Start"=segmentation_model.gr$start,"End"=segmentation_model.gr$end,"Type"=as.data.frame(values(segmentation_model.gr)@listData)$segment.mStatus),file=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_draft2.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Save a copy in case we want to come back to it later
	saveUMRLMRSegments(segs=segmentation_model.gr, GRangesFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_draft2.rds"), TableFilename=paste0(project_id,"_",meth_context,"_segmentation_model_detailed_draft2.tsv"))

	# Check whether there are any gaps in the segmentation
	#View(as.data.frame(setdiff(setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr)), segmentation_model.gr)))
	length(setdiff(setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr)), segmentation_model.gr))

	#Some segments ended up with NA in segment.mStatus (valid status was in type.1 and type.2):
	length(segmentation_model.gr[is.na(segmentation_model.gr$segment.mStatus)])
	#4 
	segmentation_model.gr[is.na(segmentation_model.gr$segment.mStatus)]$segment.mStatus = segmentation_model.gr[is.na(segmentation_model.gr$segment.mStatus)]$type.2
	length(segmentation_model.gr[is.na(segmentation_model.gr$segment.mStatus)])
	#0
	
	# 72, Schmitz: 50 gaps identified - on inspection, all are very short (<600bp) UMR segments, typically in heterochromatic regions, which do not contain any CHG sites.  Often they are in between a TEM and GBM-like segment, and seem enriched for variable sites.
	#Instantiate these as new UMR segments in the model:
	missing_segments.gr = setdiff(setdiff(genome_loci.gr, c(mask_loci.gr, segment_mask_loci.gr)), segmentation_model.gr)
	values(missing_segments.gr) = cbind(NA_by_segment(missing_segments.gr),NA_by_segment(missing_segments.gr),NA_by_segment(missing_segments.gr))
	missing_segments.gr$segment.mStatus = rep("UMR",length(missing_segments.gr))
	segmentation_model.gr = c(segmentation_model.gr, missing_segments.gr)
	
	# Write a copy for visualisation, add the sample name to the segment type
	write.table(cbind("Chromosome"=substr(as.character(segmentation_model.gr@seqnames),4,4),"Start"=segmentation_model.gr@ranges@start,"End"=segmentation_model.gr@ranges@start+segmentation_model.gr@ranges@width-1,"Type"=paste0(as.data.frame(values(segmentation_model.gr)@listData)$segment.mStatus,"_",sample_list[sample_list$SRA.Accession==sample_id,"Name"])),file=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_draft3.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Save a copy in case we want to come back to it later
	saveUMRLMRSegments(segs=segmentation_model.gr, GRangesFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_draft3.rds"), TableFilename=paste0(project_id,"_",meth_context,"_",sample_id,"_segmentation_model_detailed_draft3.tsv"))
	
	### THE SEGMENTATION IS FINISHED!  NOW NEED TO GO BACK ABOVE TO REGENERATE REPORTS ON OVERLAPS BETWEEN SEGMENTS AND VARIABLE SITES, SEGMENTS AND ANNOTATED MODELS, AND WHOLE GENOME PLOTS ###
	
	# What proportion of remaining GBM segments lays outside of annotated genes?
    
	G_x_G_hits=findOverlaps(reduce(gene_ranges), segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"])
	G_x_G_overlaps <- pintersect(reduce(gene_ranges)[queryHits(G_x_G_hits)], segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"][subjectHits(G_x_G_hits)])
	G_x_G_overlap_prop <- sum(width(G_x_G_overlaps)) / sum(width(segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]))
    G_x_G_overlap_prop 
    # 1

quit()


# This package brings in a load of methods from bedtools
# Tried using bed_intersect, but it seems to handle coordinates of segments differently to Granges package, and ends up extending all segments by one nt at 5'
#install.packages("valr")
#library(valr)

install.packages("VennDiagram")
library(VennDiagram)

model0 = readRDS("../../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_segmentation_model_draft3.rds") # Schmitz parental model
model1 = readRDS("../SRX445897/m1001_CG_SRX445897_segmentation_model_draft3.rds") # Dor-10
model2 = readRDS("../SRX446021/m1001_CG_SRX446021_segmentation_model_draft3.rds") # TOM 01
model3 = readRDS("../SRX446040/m1001_CG_SRX446040_segmentation_model_draft3.rds") # Col-0

model0_GBM=model0[model0$segment.mStatus=="GBM"]
model1_GBM=model1[model1$segment.mStatus=="GBM"]
model2_GBM=model2[model2$segment.mStatus=="GBM"]
model3_GBM=model3[model3$segment.mStatus=="GBM"]
model0_TEM=model0[model0$segment.mStatus=="TEM"]
model1_TEM=model1[model1$segment.mStatus=="TEM"]
model2_TEM=model2[model2$segment.mStatus=="TEM"]
model3_TEM=model3[model3$segment.mStatus=="TEM"]
sum(model0_GBM@ranges@width)
#[1] 17671482
sum(model1_GBM@ranges@width)
#[1] 27218126
sum(model2_GBM@ranges@width)
#[1] 24618643
sum(model3_GBM@ranges@width)
#[1] 17955065
sum(intersect(model0_GBM, model3_GBM)@ranges@width)
#[1] 16197790
sum(intersect(model1_GBM, model2_GBM)@ranges@width)
#[1] 21020273
sum(intersect(model1_GBM, model3_GBM)@ranges@width)
#[1] 16305608
sum(intersect(model2_GBM, model3_GBM)@ranges@width)
#[1] 16100187
# Col-0 accession shares >90% of GBM with Schmitz parental Col-0
sum(intersect(model0_GBM, model3_GBM)@ranges@width)/sum(model0_GBM@ranges@width)
#[1] 0.9166062
sum(intersect(model0_GBM, model3_GBM)@ranges@width)/sum(model3_GBM@ranges@width)
#[1] 0.9021293
sum(intersect(model1_GBM, model2_GBM)@ranges@width)/sum(model1_GBM@ranges@width)
#[1] 0.7722895
sum(intersect(model1_GBM, model2_GBM)@ranges@width)/sum(model2_GBM@ranges@width)
#[1] 0.8538356
sum(intersect(model1_GBM, model3_GBM)@ranges@width)/sum(model1_GBM@ranges@width)
#[1] 0.5990717
sum(intersect(model1_GBM, model3_GBM)@ranges@width)/sum(model3_GBM@ranges@width)
#[1] 0.9081342
sum(intersect(model2_GBM, model3_GBM)@ranges@width)/sum(model2_GBM@ranges@width)
#[1] 0.6539835
sum(intersect(model2_GBM, model3_GBM)@ranges@width)/sum(model3_GBM@ranges@width)
#[1] 0.8966933

sum(intersect(intersect(model1_GBM, model2_GBM), model3_GBM)@ranges@width)
#[1] 15354935

# Compare numbers of segments as well
length(model1_GBM)
#[1] 21167
length(model2_GBM)
#[1] 19817
length(model3_GBM)
#[1] 18566

pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_venn_between_samples_GBM.pdf"))
grid.newpage()
draw.pairwise.venn(area1=sum(model0_GBM@ranges@width), area2=sum(model3_GBM@ranges@width), cross.area=sum(intersect(model0_GBM, model3_GBM)@ranges@width),category=c("Schmitz parental consensus Col-0", "m1001 Col-0"), fill=c('red', 'blue'), cat.col=c('red', 'blue'), euler=TRUE, scaled=TRUE)
overrideTriple=TRUE  # Needed to make proportional diagrams work (approximately).  NB. areas are not exactlt proportionate as this is not necessarily possible with circles.
grid.newpage()
draw.triple.venn(area1=sum(model1_GBM@ranges@width), area2=sum(model2_GBM@ranges@width), area3=sum(model3_GBM@ranges@width), n12=sum(intersect(model1_GBM, model2_GBM)@ranges@width), n13=sum(intersect(model1_GBM, model3_GBM)@ranges@width), n23=sum(intersect(model2_GBM, model3_GBM)@ranges@width), n123=sum(intersect(model3_GBM, intersect(model1_GBM, model2_GBM))@ranges@width), category=c("Dor 10", "TOM 01", "Col-0"), fill=c('red', 'blue', 'green'), cat.col=c('red', 'blue', 'green'), euler=TRUE, scaled=TRUE)
dev.off()

sum(model0_TEM@ranges@width)
#[1] 26743518
sum(model1_TEM@ranges@width)
#[1] 27416574
sum(model2_TEM@ranges@width)
#[1] 27004117
sum(model3_TEM@ranges@width)
#[1] 26478873
sum(intersect(model0_TEM, model3_TEM)@ranges@width)
#[1] 25381894
sum(intersect(model1_TEM, model2_TEM)@ranges@width)
#[1] 24196216
sum(intersect(model1_TEM, model3_TEM)@ranges@width)
#[1] 23599089
sum(intersect(model2_TEM, model3_TEM)@ranges@width)
#[1] 23794128
sum(intersect(model0_TEM, model3_TEM)@ranges@width)/sum(model0_TEM@ranges@width)
#[1] 0.9490858
sum(intersect(model0_TEM, model3_TEM)@ranges@width)/sum(model3_TEM@ranges@width)
#[1] 0.9585715
sum(intersect(model1_TEM, model2_TEM)@ranges@width)/sum(model1_TEM@ranges@width)
#[1] 0.8825397
sum(intersect(model1_TEM, model2_TEM)@ranges@width)/sum(model2_TEM@ranges@width)
#[1] 0.8960195
sum(intersect(model1_TEM, model3_TEM)@ranges@width)/sum(model1_TEM@ranges@width)
#[1] 0.860761
sum(intersect(model1_TEM, model3_TEM)@ranges@width)/sum(model3_TEM@ranges@width)
#[1] 0.8912433
sum(intersect(model2_TEM, model3_TEM)@ranges@width)/sum(model2_TEM@ranges@width)
#[1] 0.8811296
sum(intersect(model2_TEM, model3_TEM)@ranges@width)/sum(model3_TEM@ranges@width)
#[1] 0.898608 

length(model1_TEM)
#[1] 4875
length(model2_TEM)
#[1] 5009
length(model3_TEM)
#[1] 5489

pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_venn_between_samples_TEM.pdf"))
grid.newpage()
draw.pairwise.venn(area1=sum(model0_TEM@ranges@width), area2=sum(model3_TEM@ranges@width), cross.area=sum(intersect(model0_TEM, model3_TEM)@ranges@width),category=c("Schmitz parental consensus Col-0", "m1001 Col-0"), fill=c('red', 'blue'), cat.col=c('red', 'blue'), euler=TRUE, scaled=TRUE)
grid.newpage()
draw.triple.venn(area1=sum(model1_TEM@ranges@width), area2=sum(model2_TEM@ranges@width), area3=sum(model3_TEM@ranges@width), n12=sum(intersect(model1_TEM, model2_TEM)@ranges@width), n13=sum(intersect(model1_TEM, model3_TEM)@ranges@width), n23=sum(intersect(model2_TEM, model3_TEM)@ranges@width), n123=sum(intersect(model3_TEM, intersect(model1_TEM, model2_TEM))@ranges@width), category=c("Dor 10", "TOM 01", "Col-0"), fill=c('red', 'blue', 'green'), cat.col=c('red', 'blue', 'green'), euler=TRUE, scaled=TRUE)
dev.off()
 

# Visualise length distributions of segments in each of the three models
pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_genome_segment_length_by_mStatus.pdf"))
segmentation_model.gr=model1
print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Width (nt)") +xlim(0,5000) +ylim(0,0.0015))
 
segmentation_model.gr=model2
print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Width (nt)") +xlim(0,5000) +ylim(0,0.0015))
 
segmentation_model.gr=model3
print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=data.frame(values(segmentation_model.gr)@listData)$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Width (nt)") +xlim(0,5000) +ylim(0,0.0015))
dev.off() 

 # Use the GenometriCorr package to develop some statistics about overlaps between the models
install.packages("devtools")
library(devtools)
install_github('favorov/GenometriCorr')
 
# alternatively:
#install.packages('GenometriCorr',repos='http://genometricorr.sourceforge.net/R/',type='source')
 
 
library("GenometriCorr")

VisualiseTwoGRanges(model1_GBM, model2_GBM)
model1_to_model2 = GenometriCorrelation(model1_GBM, model2_GBM, keep.distributions = TRUE, showProgressBar = FALSE)

GenometriCorrelation(setdiff(model1_GBM,intersect(model1_GBM, model3_GBM)), setdiff(model2_GBM,intersect(model2_GBM, model3_GBM)), keep.distributions = TRUE, showProgressBar = FALSE)

# Limit hypothesis testing to gene_ranges only:
map_space=RangedData(ranges=gene_ranges)
mapped_model1_GBM = MapRangesToGenomicIntervals(what.to.map=model1_GBM, where.to.map=gene_ranges)
mapped_model2_GBM = MapRangesToGenomicIntervals(what.to.map=model2_GBM, where.to.map=gene_ranges)
mapped_result=GenometriCorrelation(mapped_model1_GBM, mapped_model2_GBM, keep.distributions=TRUE, permut.number=permut.number)
 
 

# break data set up into batches of X samples, and write a separate script to download each batch
batch_size = 10

cat("Done.\n")


# Read the sample metadata in from a text file in the source directory
#sample_metadata=read.delim(paste0(source_dir,project_id,"_metadata.tsv"), header=TRUE, sep="\t")
#rownames(sample_metadata)=sample_metadata$Identifier

#samples = rownames(sample_metadata)
#no_samples=length(samples)

# Read the Schmitz et al, 2011 DMRs data, so this can be used in masking where needed
#Schmitz_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_CG_DMRs.txt")
#Schmitz_non_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_non-CG_DMRs.txt")
# read in sperm data
#x=read.table(file=paste0(raw_data_dir,"/GSM2323848_WT.spm.CG.tair10 (1).gff")))



no_samples = 1
sample_info=data.frame(Cov_C=integer(no_samples), Cov_T=integer(no_samples), CT_ratio=numeric(no_samples), ISNA_nonC=integer(no_samples), ISNA_C=integer(no_samples), Cprop_cutoff=numeric(no_samples))
rownames(sample_info)=sample_id

# quantile cutoff for C/T ratio in Chloroplast meth_context sites
Cprop_quantile=0.99

coverage_data = data.frame("Sample"=as.factor(NULL),"Chromosome"=as.factor(NULL),"Locus"=as.integer(NULL),"Strand"=as.factor(NULL),"Cov_C"=as.integer(NULL),"Cov_T"=as.integer(NULL))


for(sample_name in rownames(sample_info)) {
	if(opt$verbose) {cat(paste0("Summarising: ",sample_name,"\n"))}
	
		# Summarise the data for the sample
		sample_info[sample_name,"Cov_C"]=sum(all_samples_ct[,paste0(sample_name,"_C"),with=FALSE],na.rm=TRUE)
		sample_info[sample_name,"Cov_T"]=sum(all_samples_ct[,paste0(sample_name,"_T"),with=FALSE],na.rm=TRUE)
		sample_info[sample_name,"CT_ratio"]=sample_info[sample_name,"Cov_C"]/sample_info[sample_name,"Cov_T"]
		sample_info[sample_name,"ISNA_nonC"]=sum(is.na(all_samples_ct[all_samples_ct$Chromosome!="CHRC",paste0(sample_name,"_C"),with=FALSE]))
		sample_info[sample_name,"ISNA_C"]=sum(is.na(all_samples_ct[all_samples_ct$Chromosome=="CHRC",paste0(sample_name,"_C"),with=FALSE]))
		
		
		# In the methylation calling stage later, the C/(C+T) ratio in the chloroplast is used as a proxy estimate of the combined technical errors introduced by sequencing and bisulphite conversion.
		# A threshold C/(C+T) ratio is calculated here for each sample's chloroplast reads, against which to test all the other sites for whether they contain significantly more C reads than an unmethylated site, and/or significantly fewer C reads than a methylated site.
		# Two different methods of calculating the C/(C+T) threshold are proposed.
		# The first method uses the .99 quantile of chloroplast C/(C+T), on the basis that this should be a reasonable upper bound for unmethylated sites. One drawback of this method is that, for low coverage samples, a larger number of sites will fail to be called due to not being significantly more methylated than a U site while neither being significantly less methylkated than an M site.  Such sites ultimately receive an I classification.  A second drawback is that it can potentially increase the number of sites called as U which are in reality methylated at a low level.
		# The second method uses the mean of chloroplast C/(C+T), so makes the test for a U site much more stringent, but at the same time, allows a call to be made for more sites with lower coverage.  This might be justified in cases with low coverage, but where replicates are available to increase the likelihood that sites called incorrectly due to low coverage will be contradicted by the call in their replicate.

		# The third method uses a fixed value for the threshold.  In Lizzie's data, typical observed values are 0.2-0.4%
		
		# Method 1
		# Find C/(C+T) cutoff thresholds for the sample from 99% quantile of C/(C+T) distribution on Chloroplast	
		sample_info[sample_name,"Cprop_cutoff"]=quantile(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="ChrC"),paste0(sample_name,"_C"),with=FALSE])[,1]/(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="ChrC"),paste0(sample_name,"_C"),with=FALSE])[,1]+as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="ChrC"),paste0(sample_name,"_T"),with=FALSE])[,1]),Cprop_quantile,na.rm=TRUE)
		
		# Method 2
		# Find C/(C+T) ratio for chloroplast as a whole (estimation of conversion failure rate + sequencing error)
	
		# Method 3
		#sample_info[sample_name,"Cprop_cutoff"]=0.005
				
		#sample_info[sample_name,"Cprop_cutoff"]=sum(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_C"),with=FALSE])[,1],na.rm=TRUE)/(sum(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_C"),with=FALSE])[,1],na.rm=TRUE)+sum(as.data.frame(all_samples_ct[(all_samples_ct$Chromosome=="CHRC"),paste0(sample_name,"_T"),with=FALSE])[,1],na.rm=TRUE))
		
		# Add the sample coverage data to the long table
		sample_coverage_data=as.data.frame(all_samples_ct[,c("Chromosome","Locus","Strand","Cov_C"=paste0(sample_name,"_C"),"Cov_T"=paste0(sample_name,"_T")),with=FALSE])
		sample_coverage_data=cbind(rep(sample_name,nrow(sample_coverage_data)),sample_coverage_data)
		colnames(sample_coverage_data)=c("Sample","Chromosome","Locus","Strand","Cov_C","Cov_T")	
		coverage_data=rbind(coverage_data,sample_coverage_data)
		rm(sample_coverage_data)
}

	
# Dump the sample info to a file for convenience
write.table(sample_info,file=paste0(project_id,"_",meth_context,"_",sample_id,"_sample_info.tsv"),sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)
saveRDS(sample_info, file=paste0(project_id,"_",meth_context,"_",sample_id,"_sample_info.rds"))
	
# Dump the coverage data table to file for a convenience cache
saveRDS(coverage_data, file=paste0(project_id,"_",meth_context,"_",sample_id,"_coverage_data_detail.rds"))
	
coverage_data=data.table(coverage_data)


# Plot histogram of coverage
pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_coverage_histogram_CG.pdf"))
		
# plot data for all Cytosines
print(ggplot(coverage_data) + geom_histogram(aes(x=coverage_data$Cov_C+coverage_data$Cov_T), binwidth=1) + xlab("Coverage, all CG sites") +xlim(0,300))
dev.off()


# Plot histograms of methylation proportion
pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_meth_prop_histogram_CG.pdf"))
		
# plot data for all Cytosines
print(ggplot(coverage_data) + geom_histogram(aes(x=coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T)), binwidth=0.02) + xlab("C/(C+T) all CG sites"))
# plot data for all Cytosines with methylation between 0.01 and 1.01
print(ggplot(coverage_data) + geom_histogram(aes(x=coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T)), binwidth=0.02) + xlab("C/(C+T) all CG sites") +xlim(0.01,1.01))
# plot data for all Cytosines with methylation between 0.25 and 0.75
#print(ggplot(coverage_data) + geom_histogram(aes(x=coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T)), binwidth=0.02) + xlab("C/(C+T) all CG sites") +xlim(0.25,0.75))

# plot data for Cytosines with coverage above threshold
min_cov=20

print(ggplot(coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]) + geom_histogram(aes(x=coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]$Cov_C/(coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]$Cov_C+coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]$Cov_T)), binwidth=0.02) + xlab(paste0("C/(C+T) CG sites with coverage>",min_cov)))
print(ggplot(coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]) + geom_histogram(aes(x=coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]$Cov_C/(coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]$Cov_C+coverage_data[(coverage_data$Cov_C+coverage_data$Cov_T)>min_cov,]$Cov_T)), binwidth=0.02) + xlab(paste0("C/(C+T) CG sites with coverage>",min_cov) ) +xlim(0.1,1.01))


dev.off()


	fisherTestMulti <- function(input_matrix, maxN = 100, Mprop_cutoff=0, Uprop_cutoff = 0, bh = TRUE) {
	
		# This function accepts a matrix of pairs of #C and #T values, and a cutoff for errors
		# For each pair it estimates the p-values/FDR of the pair being distinct from those expected from a methylated and from an unmethylated locus
		
		# Confidence level for Fisher test
		fisher_conf_level = 0.95
		
		# Set up two matrices, one for Fisher test results to rule out Unmethylated, and one to rule out Methylated 
        fmatU <- matrix(NA, maxN+1, maxN+1)
        fmatM <- matrix(NA, maxN+1, maxN+1)
		
		# Precompute the test results for C and T values between 0 and 50 and store in matrix
		# Matrix slots will have offset of 1 from relevant value: e.g. values 0-50, slots #1-#51
		
        for (i in 0:maxN) {
            for (j in 0:maxN) {
				Cprop=i/(i+j)
		
				# Minimum number of C calls expected if methylated with errors
				# Truncate to remain on the conservative side (control FP)
				expMC=trunc((i+j)*(1-Mprop_cutoff))
				# Maximum number of T calls expected if methylated with errors
				expMT=i+j-expMC
		
				# Minimum number of T calls expected if unmethylated with errors
				# Truncate to remain on the conservative side (control FP)
				expUT=trunc((i+j)*(1-Uprop_cutoff))
				# Maximum number of C calls expected if unmethylated with errors
				expUC=i+j-expUT

				fmatU[i+1,j+1] = fisher.test(matrix(c(i, j, expUC,expUT), byrow = TRUE, 2, 2), or = 1, alternative = "greater", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
				fmatM[i+1,j+1] = fisher.test(matrix(c(i, j, expMC, expMT), byrow = TRUE, 2, 2), or = 1, alternative = "less", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
			}
		}
        pVfun <- function(x) {
			#is.na test added to cope with NA values in input matrix
			if (is.na(x[2])) {
				# Can't differentiate from either M or U
				return(c(1,1))
			} else if ((x[1] < (maxN + 1)) & (x[2] < (maxN + 1))) { 
				# C and T counts are both small so we can use precomputed p-values
				p_valueU = fmatU[x[1]+1, x[2]+1]
				p_valueM = fmatM[x[1]+1, x[2]+1]
			}  else {
				# We need to calculate the p-values
				Cprop=x[1]/(x[1]+x[2])
		
				# Minimum number of C calls expected if methylated with errors
				# Truncate to remain on the conservative side (control FP)
				expMC=trunc((x[1]+x[2])*(1-Mprop_cutoff))
				# Maximum number of T calls expected if methylated with errors
				expMT=x[1]+x[2]-expMC
		
				# Minimum number of T calls expected if unmethylated with errors
				# Truncate to remain on the conservative side (control FP)
				expUT=trunc((x[1]+x[2])*(1-Uprop_cutoff))
				# Maximum number of C calls expected if unmethylated with errors
				expUC=x[1]+x[2]-expUT

				# p-value of pair being distinct from those expected from a U site
				p_valueU = fisher.test(matrix(c(x[1], x[2], expUC,expUT), byrow = TRUE, 2, 2), or = 1, alternative = "greater", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
				# p-value of pair being distinct from those expected from a M site
				p_valueM = fisher.test(matrix(c(x[1], x[2], expMC, expMT), byrow = TRUE, 2, 2), or = 1, alternative = "less", conf.int = TRUE, conf.level = fisher_conf_level)$p.value
			}
			return(c(p_valueM, p_valueU))
		}
		
		pValues <- apply(input_matrix, 1, pVfun)
        if (bh) {
            pValues[1,] <- p.adjust(pValues[1,], method = "BH")
            pValues[2,] <- p.adjust(pValues[2,], method = "BH")
		}
		#pValues rescaling removed - I think they included this	step to turn p-value into a positive integer for convenience of storage
		#pValues <- -round(10 * log10(pValues))
        return(pValues)
    }
 	
	#sample_Cprop_cutoff is derived from the proportion of C calls in Chloroplast meth_context sites and represents an estimate of the bounds of error
	#Chloroplast should not be methylated, so any C calls (as opposed to T calls) result from one of the various sources of error
	#the sample_Cprop_cutoff is calculated by looking at the Cprop_cutoff quantile of the proportion of C calls among all Chloroplst meth_context sites in the given sample

	meth_data_fisher=matrix(, nrow=nrow(coverage_data), ncol=2)
		
	for(sample_name in rownames(sample_info)) {
		cat(paste0("Fisher's exact tests for sample ",sample_name,"\n"))
		meth_data_fisher[coverage_data$Sample==sample_name,] = t(fisherTestMulti(as.matrix(cbind(coverage_data[coverage_data$Sample==sample_name,]$Cov_C,coverage_data[coverage_data$Sample==sample_name,]$Cov_T)), Mprop_cutoff=sample_info[rownames(sample_info)==sample_name,]$Cprop_cutoff, Uprop_cutoff=sample_info[rownames(sample_info)==sample_name,]$Cprop_cutoff))
	}	

	# We now have a matrix of pairs of FDR-adjusted p-values, one pair per site, indicating whether the site is significantly distinct from the expected #C and #T coverage of a Methylated or from an Unmethylated site with errors
	
	# Save the Fisher's test results to cache file
	saveRDS(meth_data_fisher, file=paste0(project_id,"_",meth_context,"_",sample_id,"_meth_data_fisher_2.rds"))

	# We have used the Fisher's Exact test to identify sites which can't reliably be distinguished from U or M status
	# Now we use the binomial test to identify sites that can be reliably inferred to be Unmethylated
	

	# Function to calculate a matrix of Binomial p-values given an input matrix of C-coverage and C+T coverage values for loci
	# Original function taken from methylPipe$binomTestMulti.  Modifications are commented below.
    binomTestMulti <- function(mat, maxN = 100, p = bc, bh = TRUE) {
        bmat <- matrix(NA, maxN, maxN)
        for (i in 1:maxN) {
            for (j in i:maxN) bmat[i, j] <- binom.test(x = i, 
                n = j, p = p, alternative = "greater")$p.value
        }
        pVfun <- function(x) {
			#is.na test added to cope with NA values in input matrix
			if (is.na(x[2])) {
				return(NA)
			#x[1] == 0 test added after discussion with methylPipe author: https://support.bioconductor.org/p/102762/#102765
			} else if (x[1] == 0) {
				return(1)
			} else if (x[2] < (maxN + 1)) { 
                return(bmat[x[1], x[2]])
			}
            else {
				return(binom.test(x = x[1], n = x[2], p = p, 
                alternative = "greater")$p.value)
			}
		}
        pValues <- apply(mat, 1, pVfun)
        if (bh) {
            pValues <- p.adjust(pValues, method = "BH")
		}
		#pValues rescaling removed - they included this step to turn p-value into a positive integer for convenience of storage
		#pValues <- -round(10 * log10(pValues))
        return(pValues)
    }
 	
	meth_data_binom=matrix(, nrow=nrow(coverage_data), ncol=1)
	
	for(sample_name in rownames(sample_info)) {
		cat(paste0("Binomial tests for sample ",sample_name,"\n"))
		meth_data_binom[coverage_data$Sample==sample_name,] = binomTestMulti(as.matrix(cbind(coverage_data[coverage_data$Sample==sample_name,]$Cov_C,coverage_data[coverage_data$Sample==sample_name,]$Cov_C+coverage_data[coverage_data$Sample==sample_name,]$Cov_T)), p=sample_info[rownames(sample_info)==sample_name,]$Cprop_cutoff)
	}	

	# Save the Binomial test results to cache file
	saveRDS(meth_data_binom, file=paste0(project_id,"_",meth_context,"_",sample_id,"_meth_data_binom.rds"))


	#View(cbind(coverage_data,meth_data_fisher,meth_data_binom)[1:100,])

	
	# p-value/FDR for rejecting hypotheses of 'Methylated' or 'Unmethylated'
	# Making this limit more strict (smaller p-value) has the effect of excluding more sites due to insufficient coverage
	# 0.01 may be a better threshold to use when no replicates are available.  Where replicates are available, errors made in this step will be compensated for by comparison between replicates and discard of sites where replicates disagree
	#fisher_p_value_cutoff = 0.01
	# 0.01 may be too strict in removing low-coverage sites.  We use 0.05 and allow disagreement between replicates to control for errors made in this call.
	#fisher_p_value_cutoff = 0.05
	# Tighter definition of 'I' to remove more errors:
	fisher_p_value_cutoff = 0.005 

	# p-value/FDR for testing whether sites have too much methylation to be considered Unmethylated
	# Making this limit more strict (smaller p-value) has the effect of calling more sites (with a few methylated reads) as Unmethylated, rather than Methylated/Partial
	# Lister et al, 2009 used binomial test with B-H adjusted p-values for FDR<0.01
	# In their work, anything that passed this test was classified as having some methylation, but we will instead use failure to pass this test as a criterion to call sites Unmethylated
	# 0.01 still leaves a mode of mCG ~ 0.1 in the density distribution of 'partially' methylated sites.  0.005 removes this mode, so is likely to return more true positive 'U' sites.
	#binomial_p_value_threshold=0.01	
	#binomial_p_value_threshold=0.005	
	# Tighter definition of 'U' to be more certain of Us:
	binomial_p_value_threshold=0.05	
	#binomial_p_value_threshold=0.001	

	# We also impose a strict cutoff on minimum coverage before a site will be considered for calling
	# 3 reads is a very generous minimum, allowing the most sites to be considered.  10 reads would be a conservative minimum, if the experiment contains enough data.  There is literature support for thresholds in this range:
	# Lister et al, 2009 >=3 reads coverage
	# Ziller et al, 2013 >=5 reads coverage
	# Kunde-Ramamoorthy et al, 2013 >=10 reads coverage
	# Cokus et al, 2008 >=5 reads coverage.
	#binomial_coverage_cutoff=3
	# Tighter coverage constraint to remove uncertainty:
	binomial_coverage_cutoff=10
	# Second constraint, to rule out 1:k sites from being called as U if k<24 (called I instead)
	binomial_coverage_cutoff2 = binomial_coverage_cutoff * 2.5

	
	# Using either the Fisher or binomial method alone, a large number of CG dinucleotides are identified as Methylated at one nucleotide and Unmethylated at the other.  This typically happens in cases where one of the sites in the pair has high coverage, and a clear 'U' or 'M' call, and the other has low coverage and a borderline call. For this reason, we use both methods to make a final call - we use Fisher's test to identify sites with ambiguous calls, which either pass or fail both tests (Unmethylated and Methylated).  We use the binomial test to identify sites clearly Unmethylated.  For the remainder (sites not Unmethylated, and with sufficient but not too high coverage), we consider them Methylated if more than a cutoff fraction of reads are methylated, otherwise we consider them partially methylated.  At the moment this fraction is chosen arbitrarily.  Ideally we could optimise this cutoff, however, while reducing it to maximise coverage and minimise sites lost as 'Partially methylated', we would likely also increase the number of adjacent CG dinucleotides showing inconsistent methylation calls.	
	
	# We impose a cutoff for the minimum proportion #C/(#C+#T) we accept to call the locus methylated. This differs from Lister et al, 2009 approach where they allowed the binomial test alone to call site methylation status.
	# According to Gardiner et al, 2015, Cokus et al, 2008 had identified that it is generally possible to impose a far higher proportion as the threshold in the CG context (they used 80%) than in other contexts (they used 25% and 10% for CHG and CHH respectively)
	### We should be able to work out this threshold empirically for CG context using ROC analysis on the basis that methylation status in adjacent CG pairs should be concordant
	### We might also consider whether different thresholds are appropriate dependant on whether each C is considered in isolation, or whether CG dinucleotides are considered as a single site
	
	# New approach to this: we try a range of values from 0.1 to 0.9 for partial_meth_ctuff, and instead of classifying sites below the line as 'P' we classify them as 'U'.  This increases the concordance among reps.  We maximise the mean proportion of sites with concordance among reps in the script 5-methylation_calls_optimise_cutoff_2018-02-23.R.  This comes out with a threshold of mCG=0.45 for the Schmitz data set.  We then use this cutoff, along with a more stringent binomial test threshold to distinguish unmethylated from methylated sites.
	
	partial_meth_cutoff=0
	if (meth_context=="CG") {
		#partial_meth_cutoff = 0.45
		# New narrower definition of 'P' to allow more heterozygous sites to be called 'M'
		partial_meth_cutoff = 0.25
	} else if (meth_context=="CHG") {
		partial_meth_cutoff = 0.25
	} else if (meth_context=="CHH") {
		partial_meth_cutoff = 0.1
	}
	
	# Create a matrix ready to hold results, for efficiency
	#meth_data=matrix(, nrow=nrow(coverage_data), ncol=2)
	# Leave out p-values to save memory
	meth_data=matrix(, nrow=nrow(coverage_data), ncol=1)
	
	# Classify each site dependent on which of the two Fisher tests passes the FDR threshold
	fisher_class=ifelse((meth_data_fisher[,1]<fisher_p_value_cutoff),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),"P","U"),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),"M","I"))
	# Classify each site dependent on whether it passes the binomial test

	#binom_class=ifelse((coverage_data$Cov_C+coverage_data$Cov_T)<binomial_coverage_cutoff,"I",ifelse(meth_data_binom<binomial_p_value_threshold,"M","U"))
	#NEW VERSION OF BINOM_CLASS THAT HAS A MORE STRINGENT COVERAGE CUT OFF WHEN 1 SITE IS C
	# The basis of this extra constraint is that if a site has only a single read signalling methylated, and coverage overall is not very high, we discard the site as unreliable to call using the binomial test (it will probably be classified as U, due to only one read signalling M, but will have a relatively high %mCG due to low coverage)
	binom_class=ifelse((coverage_data$Cov_C+coverage_data$Cov_T)<binomial_coverage_cutoff,"I"
                   ,ifelse(coverage_data$Cov_C==1 & (coverage_data$Cov_C+coverage_data$Cov_T)<binomial_coverage_cutoff2, "I"
                           ,ifelse(meth_data_binom<binomial_p_value_threshold,"M","U")))

	# Classify each site based on the relevant hybrid assessment
	# Added new test here for binom_class=="I" so that the minimum coverage cutoff is enforced properly
	#meth_data[,1] = ifelse(binom_class=="I","I",ifelse(fisher_class=="I","I",ifelse(binom_class=="U","U",ifelse((coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T))>=partial_meth_cutoff,"M","P"))))
	# version without p-values to save memory:
	meth_data = ifelse(binom_class=="I","I",ifelse(fisher_class=="I","I",ifelse(binom_class=="U","U",ifelse((coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T))>=partial_meth_cutoff,"M","P"))))
	
	#Alternative approach: consider mC<partial_meth_cutoff calls to be U, rather than P
	#meth_data[,1] = ifelse(fisher_class=="I","I",ifelse(binom_class=="U","U",ifelse((coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T))>=partial_meth_cutoff,"M","U")))

	rm(fisher_class, binom_class)
	
	# Choose the relevant p-value to retain (from binomial or fishers test), dependent on which status was chosen
	#### This needs more work, but not done now, because at the moment we do not make further use of this composite p-value
	#### Fisher's method (using a chisquared test) could be used to combine multiple p-values addressing the same general hypothesis.  Implemented in chiCombP function in the methylPipe package.
	# Removed to save memory by not retaining p-values:
	#meth_data[,2] = ifelse((meth_data_fisher[,1]<fisher_p_value_cutoff),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),pmax(meth_data_fisher[,1],meth_data_fisher[,2]),meth_data_fisher[,1]),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),meth_data_fisher[,2],pmin(meth_data_fisher[,1],meth_data_fisher[,2])))

	# Discard status and replace with "I" if coverage is too high to consider trusting the call (5 sigma)
	# This assumes that coverage is roughly Gaussian distributed (although it typically is long-tailed)
	for(sample_name in rownames(sample_info)) {
		coverage_upper_cutoff = sd(coverage_data[Sample==sample_name,Cov_C+Cov_T],na.rm=TRUE)*5
		#meth_data[(coverage_data$Sample==sample_name) & (!is.na(coverage_data$Cov_C)) & ((coverage_data$Cov_C+coverage_data$Cov_T)>coverage_upper_cutoff),1] = "I"
		# revised version with no p-values:
		meth_data[(coverage_data$Sample==sample_name) & (!is.na(coverage_data$Cov_C)) & ((coverage_data$Cov_C+coverage_data$Cov_T)>coverage_upper_cutoff)] = "I"
	}
	

	# Organise meth_data from the large matrix to a data frame
	# no need if not storing p-values, but we are going to anyway so it has a colname for use in the ggplotting
	meth_data=data.frame(meth_data)
	colnames(meth_data) = c("meth_status")
	#meth_data$p_value=as.numeric(levels(meth_data$p_value)[meth_data$p_value])

	
	# Dump the methylation calls data frame to file for a convenience cache
	saveRDS(meth_data, file=paste0(project_id,"_",meth_context,"_",sample_id,"_meth_data_",fisher_p_value_cutoff,"_",binomial_p_value_threshold,"_",binomial_coverage_cutoff,"_",binomial_coverage_cutoff2,".rds"))
	
  	
	# Plot the distribution of C proportions at sites with #C more or less than half coverage 
	### scaffolds mess up these plots
	#pdf(paste0(project_id,"_",meth_context,"_site_C_proportion_ge_0.5_by_sample.pdf"))
	#print(ggplot(cbind(coverage_data,meth_data)[(meth_status!="I") & (Cov_C>=Cov_T),], aes(x=Cov_C/(Cov_C+Cov_T),colour=Chromosome)) + geom_freqpoly(binwidth = .01) + facet_wrap(~Sample))
	#dev.off()
	#pdf(paste0(project_id,"_",meth_context,"_site_C_proportion_lt_0.5_by_sample.pdf"))
	#print(ggplot(cbind(coverage_data,meth_data)[(meth_status!="I") & (Cov_C<Cov_T),], aes(x=Cov_C/(Cov_C+Cov_T),colour=Chromosome)) + geom_freqpoly(binwidth = .01) + facet_wrap(~Sample))
	#dev.off()
	
	# Plot C/T coverage distribution overall, and by methylation status
	plot_coverage_limit=100
	pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_C_T_coverage_distribution_smooth.pdf"))
	print(smoothScatter(x=coverage_data[(coverage_data$Cov_T<=plot_coverage_limit) & (coverage_data$Cov_C<=plot_coverage_limit),]$Cov_T, y=coverage_data[(coverage_data$Cov_T<=plot_coverage_limit) & (coverage_data$Cov_C<=plot_coverage_limit),]$Cov_C))
	dev.off()
	

	# This version can't manage to allocate enough memory for the cbind
	#pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_point.pdf"))
	#print(ggplot(cbind(coverage_data,meth_data)[sample(nrow(coverage_data),900000),][!is.na(Cov_C),], aes(x=Cov_T, y=Cov_C)) + geom_point(aes(colour = meth_status), alpha=0.05, size=2) + #xlim(0,plot_coverage_limit) + ylim(0,plot_coverage_limit) + theme_minimal())
	#dev.off()
	z=sample(nrow(coverage_data),90000)
	pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_C_T_coverage_distribution_point.pdf"))
	print(ggplot(cbind(coverage_data[z,],meth_data[z,])[!is.na(Cov_C),], aes(x=Cov_T, y=Cov_C)) + geom_point(aes(colour = V2), alpha=0.05, size=2) + xlim(0,plot_coverage_limit) + ylim(0,plot_coverage_limit) + theme_minimal())
	dev.off()

	# Plot C/T coverage distribution as classified by Binomial test
	# This version can't manage to allocate enough memory for the cbind
	#pdf(paste0(project_id,"_",meth_context,"_C_T_coverage_distribution_binomial_sig.pdf"))
	#print(ggplot(cbind(coverage_data,meth_data,meth_data_binom,binom_class=ifelse(is.na(meth_data_binom),"I",ifelse((coverage_data$Cov_C+coverage_data$Cov_T)<binomial_coverage_cutoff,"I",ifelse(meth_data_binom<binomial_p_value_threshold,"M","U"))))[sample(nrow(coverage_data),300000),][(!is.na(Cov_C)),], aes(x=Cov_T, y=Cov_C)) + geom_point(aes(colour = binom_class.V1), alpha=0.05, size=2) + xlim(0,plot_coverage_limit) + ylim(0,plot_coverage_limit) + theme_minimal())
	#dev.off()
	pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_C_T_coverage_distribution_binomial_sig.pdf"))
	print(ggplot(cbind(coverage_data[z,],meth_data[z,],meth_data_binom[z],binom_class=ifelse(is.na(meth_data_binom[z]),"I",ifelse((coverage_data[z,]$Cov_C+coverage_data[z,]$Cov_T)<binomial_coverage_cutoff,"I",ifelse(meth_data_binom[z]<binomial_p_value_threshold,"M","U"))))[(!is.na(Cov_C)),], aes(x=Cov_T, y=Cov_C)) + geom_point(aes(colour = binom_class), alpha=0.05, size=2) + xlim(0,plot_coverage_limit) + ylim(0,plot_coverage_limit) + theme_minimal())
	dev.off()

	
	# Generate the 'all_samples_meth_status' format that is used by subsequent analysis steps
	meth_data=data.table(meth_data)
	  
	  #### Build per-site tables of methylation calls across all samples
	  
	  # - characters in sample IDs are problematic as column names.  Replace them with _ characters
	  rownames(sample_info)= gsub("-","_",rownames(sample_info))
	  sample_info$Sample = rownames(sample_info)
	  sample_info$Identifier = rownames(sample_info)
	  coverage_data$Sample = gsub("-","_",coverage_data$Sample)
	  #colnames(all_samples_ct) = gsub("-","_",colnames(all_samples_ct))
	  
	  # Create empty tables to accumulate each of the 'all samples' columns data
	  all_samples_meth_status=NULL
	  all_samples_p_value=NULL
	  
	  sample_no=0
	  for(sample_name in rownames(sample_info)) {
	    sample_no = sample_no+1
	    if(opt$verbose) {cat(paste0("Adding methylation calls to all-samples per-site table: ",sample_name,"\n"))}
	    
	    sample_meth_status=NULL
	    sample_p_value=NULL
	    
	    # Make separate tables for the sample's meth_call, p_value and meth_call_basis, keyed on locus
	    sample_meth_status=cbind(coverage_data,meth_data)[coverage_data$Sample==sample_name,c(2,3,4,7)]
	    colnames(sample_meth_status)[4]=sample_name
	    # remopve this if not doing p-values:
	    #sample_p_value=cbind(coverage_data,meth_data)[coverage_data$Sample==sample_name,c(2,3,4,8)]
	    #colnames(sample_p_value)[4]=sample_name
	    
	    # Merge each of the sample's 3 tables into the relevant accumulation table
	    if (sample_no==1) {
	      all_samples_meth_status=sample_meth_status
	      #all_samples_p_value=sample_p_value
	    } else {
	      all_samples_meth_status=merge(all_samples_meth_status, sample_meth_status, by=c("Chromosome","Locus","Strand"), all=TRUE)
	      #all_samples_p_value=merge(all_samples_p_value, sample_p_value, by=c("Chromosome","Locus","Strand"), all=TRUE)
	    }
	    
	    # Tidy up the temporary objects
	    rm(sample_meth_status, sample_p_value)
	  } # end for each sample
	  
	  #Add a column to identify the number of samples with missing data
	  #######CG_pair_meth_status$m_count <- apply(CG_pair_meth_status, 1, function(x) sum(is.na(x))/2)
	  
	  if(opt$verbose) {
	    ### Numbers of consistent calls by replicate set:
	  }
	  
	  # Sort the data frames sensibly (locus has received an alpha rather than numerical sort at this point)
	  # Dump the site all-samples meth status data table to file for a convenience cache
	  all_samples_meth_status=all_samples_meth_status[order(all_samples_meth_status$Chromosome,all_samples_meth_status$Locus),]
	  #saveRDS(all_samples_meth_status, file=paste0(project_id,"_",meth_context,"_all_samples_meth_status.rds"))
	  saveRDS(all_samples_meth_status, file=paste0(project_id,"_",meth_context,"_",sample_id,"_all_samples_meth_status_",fisher_p_value_cutoff,"_",binomial_p_value_threshold,"_",binomial_coverage_cutoff,"_",binomial_coverage_cutoff2,".rds"))
	  #all_samples_p_value=all_samples_p_value[order(all_samples_p_value$Chromosome,all_samples_p_value$Locus),]
	
  
	
