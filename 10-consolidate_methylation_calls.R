#!/usr/bin/env Rscript

# 10-consolidate_methylation_calls.R
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

not_processed_samples = c()
processed_samples = c()
all_samples_meth_status = NA
no_samples = 0
for (this_sample in sample_list$SRA.Accession) {
  this_file = paste0(project_id,"_",meth_context,"_",this_sample,"_all_samples_meth_status_0.005_0.05_10_25.rds")
  if (file.exists(this_file)) {
	if(opt$verbose) {cat(paste0("Loading ",this_sample,"\n"))}
    no_samples = no_samples + 1
    processed_samples=c(processed_samples, this_sample)
	this_sample_meth_status = data.frame(readRDS(this_file))
    if (no_samples == 1) {
	  all_samples_meth_status = this_sample_meth_status
	} else {
	  all_samples_meth_status = merge(all_samples_meth_status, this_sample_meth_status, by=c("Chromosome","Locus","Strand"), all=TRUE)
    }
  } else {
    not_processed_samples=c(not_processed_samples, this_sample)
  }
}

# Dump the data table to file for a convenience cache
saveRDS(all_samples_meth_status, file=paste0(project_id,"_",meth_context,"_all_samples_meth_status_0.005_0.05_10_25.rds"))
# saved for batches 1-29 so far
#all_samples_meth_status = readRDS(file=paste0(project_id,"_",meth_context,"_all_samples_meth_status_0.005_0.05_10_25.rds"))


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


# Count X calls across all samples
no_mutations = matrix(0, nrow=length(valid_samples), ncol=2)
for(sample_no in 1:length(valid_samples)) {
  sample_name = valid_samples[sample_no]
  no_mutations[sample_no,1] = sample_name
  no_mutations[sample_no,2] = sum(ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse(all_samples_meth_status[,sample_name]=="X",1,0)))
}

# This version looks for M/U calls across all samples, independent of rep structure
### This produces FALSE almost everywhere due to missing data etc
sample_no=0
first_sample = TRUE
for(sample_name in valid_samples) {
  #sample_no = sample_no+1
  if(opt$verbose) {cat(paste0("Checking methylation calls in all-samples per-site table: ",sample_name,"\n"))}
  if (first_sample) {
    #all_m_or_u_regardless=ifelse(((all_samples_meth_status[,sample_no+3]=="U") | (all_samples_meth_status[,sample_no+3]=="M")),TRUE,FALSE)
    all_m=ifelse(is.na(all_samples_meth_status[,sample_name]),TRUE,ifelse(all_samples_meth_status[,sample_name]=="M",TRUE,FALSE))
    all_u=ifelse(is.na(all_samples_meth_status[,sample_name]),TRUE,ifelse(all_samples_meth_status[,sample_name]=="U",TRUE,FALSE))
    first_sample=FALSE
  } else {
    #all_m_or_u_regardless=all_m_or_u_regardless & ifelse(((all_samples_meth_status[,sample_no+3]=="U") | (all_samples_meth_status[,sample_no+3]=="M")),TRUE,FALSE)
    all_m=all_m & ifelse(is.na(all_samples_meth_status[,sample_name]),TRUE,ifelse(all_samples_meth_status[,sample_name]=="M",TRUE,FALSE))
    all_u=all_u & ifelse(is.na(all_samples_meth_status[,sample_name]),TRUE,ifelse(all_samples_meth_status[,sample_name]=="U",TRUE,FALSE))
  }
}

# Identify parental consensus calls
parental_consensus = ifelse((average_methylation[,1]/average_methylation[,2])>0.5, "M", "U")


# Identify sites where at least one offspring line differs from parental consensus, and also, count how many offspring differ at the locus
line_gen_no=0
first_line_gen = TRUE
site_varies_from_parentals = NA
num_lines_varying_from_parentals = NA
num_lines_m_to_u_from_parentals = NA
num_lines_u_to_m_from_parentals = NA
num_lines_m_or_u_excl_parentals = NA
for(sample_name in valid_samples) {
#for(line_gen in gen3_30) {
#  line_gen_no = line_gen_no + 1
  if(opt$verbose) {cat(paste0("Checking variation in methylation calls from parental consensus: ",sample_name,"\n"))}
  if (first_line_gen) {
    #all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    site_varies_from_parentals=ifelse((all_samples_meth_status[,sample_name] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_samples_meth_status[,sample_name]!=parental_consensus),TRUE,FALSE)
    num_lines_varying_from_parentals=ifelse((all_samples_meth_status[,sample_name] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_samples_meth_status[,sample_name]!=parental_consensus),1,0)
    num_lines_m_to_u_from_parentals=ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse((all_samples_meth_status[,sample_name] == "U") & (parental_consensus == "M"),1,0))
    num_lines_u_to_m_from_parentals=ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse((all_samples_meth_status[,sample_name] == "M") & (parental_consensus == "U"),1,0))
    num_lines_m_or_u_excl_parentals=ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse((all_samples_meth_status[,sample_name] == "M") | (all_samples_meth_status[,sample_name] == "U"),1,0))
    first_line_gen=FALSE
  } else {
    #all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    site_varies_from_parentals=site_varies_from_parentals | ifelse((all_samples_meth_status[,sample_name] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_samples_meth_status[,sample_name]!=parental_consensus),TRUE,FALSE)
    num_lines_varying_from_parentals=num_lines_varying_from_parentals + ifelse((all_samples_meth_status[,sample_name] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_samples_meth_status[,sample_name]!=parental_consensus),1,0)
    num_lines_m_to_u_from_parentals=num_lines_m_to_u_from_parentals + ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse((all_samples_meth_status[,sample_name] == "U") & (parental_consensus == "M"),1,0))
    num_lines_u_to_m_from_parentals=num_lines_u_to_m_from_parentals + ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse((all_samples_meth_status[,sample_name] == "M") & (parental_consensus == "U"),1,0))
    num_lines_m_or_u_excl_parentals=num_lines_m_or_u_excl_parentals + ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse((all_samples_meth_status[,sample_name] == "M") | (all_samples_meth_status[,sample_name] == "U"),1,0))
    #first_line_gen=FALSE
  }
}


### This is mostly bobbins
if (opt$verbose) {
  cat(paste("No of",meth_context,"loci:",nrow(all_samples_meth_status),"\n"))
  #cat(paste("No of",meth_context,"loci with all calls M or U:",nrow(all_samples_meth_status[all_m_or_u==TRUE,]),"\n"))
  cat(paste("No of",meth_context,"loci with all calls M:",nrow(all_samples_meth_status[all_m==TRUE,]),"\n"))
  cat(paste("No of",meth_context,"loci with all calls U:",nrow(all_samples_meth_status[all_u==TRUE,]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls (all parentals agree, all lines have agreement between reps):",nrow(all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3:",nrow(all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3==TRUE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3:",nrow(all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3==TRUE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3_30==TRUE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3_30==TRUE)),]),"\n"))
  cat(paste("No of",meth_context,"loci with any variation in M/U calls (at least 2 parentals agree, at least one progeny line disagrees and has agreement between reps):",nrow(all_samples_meth_status[site_varies_from_parentals,]),"\n"))
  #cat(paste("No of",meth_context,"loci with any variation in M/U calls, M at generation 3:",nrow(all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"\n"))
  #cat(paste("No of",meth_context,"loci with any variation in M/U calls, U at generation 3:",nrow(all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"\n"))
  cat(paste("No of",meth_context,"loci with a clear call in parental consensus and all offspring:",nrow(cbind.data.frame(parental_consensus, site_varies_from_parentals, num_lines_m_to_u_from_parentals, num_lines_u_to_m_from_parentals, num_lines_m_or_u_excl_parentals)[(num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),]),"\n"))
  # Commented out reporting of generation 31 changes for now.  We should redo this analysis in a more sensitive way - check each generation 31 line for changes compared to the corresponding generation 30 line, in concert with reps, but ignoring coverage and calls in generation 3 and in all other generation 30/31 lines.
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3_30==TRUE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3_30==TRUE)),]),"\n"))
}

# Generate tables of sites with various patterns of calls
across_samples_meth_status=cbind(all_samples_meth_status[,1:3],"average_methylation"=average_methylation[,1]/average_methylation[,2])[average_methylation[,2]>0,]
#variant_calls=all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE)),]
variant_calls=all_samples_meth_status[site_varies_from_parentals,]

# removed these lot for space reasons and because they never get used again:

#M_to_U_calls=all_samples_meth_status[site_varies_from_parentals & (num_lines_m_to_u_from_parentals>0),]
#U_to_M_calls=all_samples_meth_status[site_varies_from_parentals & (num_lines_u_to_m_from_parentals>0),]
#all_M_calls=all_samples_meth_status[all_m==TRUE,]
#all_U_calls=all_samples_meth_status[all_u==TRUE,]
#all_clean_calls=all_samples_meth_status[(num_lines_m_or_u_excl_parentals==1001) & (parental_consensus %in% c("M","U")),]

#M_parent_calls = all_samples_meth_status[parental_consensus=="M",]
#U_parent_calls = all_samples_meth_status[parental_consensus=="U",]



# Make loci 2nt long for CG sites
if (meth_context=="CG") {
  across_samples_meth_status=within(across_samples_meth_status, {Locus2=Locus+1})
  variant_calls=within(variant_calls, {Locus2=Locus+1})
  #M_to_U_calls=within(M_to_U_calls, {Locus2=Locus+1})
  #U_to_M_calls=within(U_to_M_calls, {Locus2=Locus+1})
  #all_M_calls=within(all_M_calls, {Locus2=Locus+1})
  #all_U_calls=within(all_U_calls, {Locus2=Locus+1})
  #all_clean_calls=within(all_clean_calls, {Locus2=Locus+1})
  
  #M_parent_calls=within(M_parent_calls, {Locus2=Locus+1})
  #U_parent_calls=within(U_parent_calls, {Locus2=Locus+1})
} else {
  across_samples_meth_status=within(across_samples_meth_status, {Locus2=Locus})
  variant_calls=within(variant_calls, {Locus2=Locus})
  #M_to_U_calls=within(M_to_U_calls, {Locus2=Locus})
  #U_to_M_calls=within(M_to_U_calls, {Locus2=Locus})
  #all_M_calls=within(all_M_calls, {Locus2=Locus})
  #all_U_calls=within(all_U_calls, {Locus2=Locus})
  #all_clean_calls=within(all_clean_calls, {Locus2=Locus})
  
  #M_parent_calls=within(M_parent_calls, {Locus2=Locus})
  #U_parent_calls=within(U_parent_calls, {Locus2=Locus})
}

# Get the genes annotation and start carving up the CG sites by gene

## Load the annotation
annot_gff = read.delim(reference_gff, header=F, comment.char="#")	
gff.exons = read.delim(reference_exons, header=F, comment.char="#")
gff.introns = read.delim(reference_introns, header=F, comment.char="#")
gff.5UTR = read.delim(reference_5UTR, header=F, comment.char="#")
gff.3UTR = read.delim(reference_3UTR, header=F, comment.char="#")

# Grab the portion relating to genes
gff.genes = annot_gff[annot_gff[,3]=="gene",]
# Grab the portion relating to transposons
gff.transposons = annot_gff[annot_gff[,3]=="transposable_element",]

rm(annot_gff)

# Convert chromosome names to uppercase to match previous objects
### All previous objects are NOT in capitals
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
gff.exons$gene_ID=str_split_fixed(str_split_fixed(gff.exons$V9, ';',3)[,1],"=",2)[,2]
gff.introns$gene_ID=str_split_fixed(str_split_fixed(gff.introns$V9, ';',3)[,1],"=",2)[,2]
gff.5UTR$gene_ID=str_split_fixed(str_split_fixed(gff.5UTR$V9, ';',3)[,1],"=",2)[,2]
gff.3UTR$gene_ID=str_split_fixed(str_split_fixed(gff.3UTR$V9, ';',3)[,1],"=",2)[,2]

#Gene names are tricky - need to find the element containing "symbol="
gff.genes=within(gff.genes, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,2],"=",2)[,2]})
gff.transposons=within(gff.transposons, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,3],"=",2)[,2]})

# Fix the gene names in exon, intron and UTRs
gff.exons$exon_ID = gff.exons$gene_ID
gff.exons$gene_ID = substr(gff.exons$gene_ID,1,9)
gff.introns$intron_ID = gff.introns$gene_ID
gff.introns$gene_ID = substr(gff.introns$gene_ID,1,9)
gff.5UTR$UTR_ID = gff.5UTR$gene_ID
gff.5UTR$gene_ID = substr(gff.5UTR$gene_ID,1,9)
gff.3UTR$UTR_ID = gff.3UTR$gene_ID
gff.3UTR$gene_ID = substr(gff.3UTR$gene_ID,1,9)

# Add the exon and intron numbers
gff.exons$exon_no = as.numeric(str_split_fixed(gff.exons$exon_ID, ':',3)[,3])
gff.introns$exon_no = as.numeric(str_split_fixed(gff.introns$intron_ID, ':',3)[,3])

# Find out how many exons each gene has
#install.packages("sqldf")
#library(sqldf)
gff.genes=merge(gff.genes, sqldf('SELECT genes.gene_id, MAX(exon_no) AS no_exons FROM [gff.exons] exons, [gff.genes] genes WHERE genes.gene_ID=exons.gene_ID GROUP BY genes.gene_id'), by="gene_ID", all=TRUE)

library(GenomicRanges)

# Make a GRanges for gene space
gene_ranges=makeGRangesFromDataFrame(df = gff.genes, start.field = "V4", end.field = "V5",seqnames.field = "V1")
# Make a GRanges for transposon space
transposon_ranges=makeGRangesFromDataFrame(df = gff.transposons, start.field = "V4", end.field = "V5",seqnames.field = "V1")
# Make a GRanges for exon space
exon_ranges=makeGRangesFromDataFrame(df = gff.exons, start.field = "V4", end.field = "V5",seqnames.field = "V1")
# Make a GRanges for intron space
intron_ranges=makeGRangesFromDataFrame(df = gff.introns, start.field = "V4", end.field = "V5",seqnames.field = "V1")
# Make a GRanges for 5' UTR space
UTR5_ranges=makeGRangesFromDataFrame(df = gff.5UTR, start.field = "V4", end.field = "V5",seqnames.field = "V1")
# Make a GRanges for 3' UTR space
UTR3_ranges=makeGRangesFromDataFrame(df = gff.3UTR, start.field = "V4", end.field = "V5",seqnames.field = "V1")


#  Make a GRanges for all CG sites
CG_site_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for all_M sites
all_M_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[all_m==TRUE,], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for all_U sites
all_U_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[all_u==TRUE,], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for variable sites
variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals,], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for M->U sites
m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (num_lines_m_to_u_from_parentals>0),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for U->M sites
u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (num_lines_u_to_m_from_parentals>0),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for clean call sites
### Commented out as not relevant here
#clean_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_m_or_u_excl_parentals==1001) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make 'clean' versions of the site ranges, where data is available for all lines
# No need to make for all_M and all_U as they were already 'clean'
#clean_all_M_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(all_m==TRUE) & (num_lines_m_or_u_excl_parentals==4),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
#clean_all_U_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(all_u==TRUE) & (num_lines_m_or_u_excl_parentals==4),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
#clean_variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
#clean_m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
#clean_u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make single change clean versions of the variable site ranges - sites where data is available for all lines, and where a change occurs but only in one offspring line
#clean_single_variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
#clean_single_m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (no_m_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
#clean_single_u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (no_u_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")


### These fail due to too many NAs:
M_parent_call_ranges = makeGRangesFromDataFrame(df = all_samples_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
U_parent_call_ranges = makeGRangesFromDataFrame(df = all_samples_meth_status[(parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")

### Need new set of DMRs relevant to 1001 methylomes paper
DMRs.gr = reduce(c(makeGRangesFromDataFrame(df=read.table(file=Schmitz_CG_DMRs, header=FALSE, sep="\t"), start.field = "V2", end.field = "V3",seqnames.field = "V1"),makeGRangesFromDataFrame(df=read.table(file=Schmitz_non_CG_DMRs, header=FALSE, sep="\t"), start.field = "V2", end.field = "V3",seqnames.field = "V1")))
levels(DMRs.gr@seqnames@values) = toupper(as.character(levels(DMRs.gr@seqnames@values)))
DMRs.gr@seqinfo@seqnames = levels(DMRs.gr@seqnames@values)

variant_call_ranges = setdiff(variant_call_ranges, DMRs.gr, ignore.strand=TRUE)
m_to_u_call_ranges = setdiff(m_to_u_call_ranges, DMRs.gr, ignore.strand=TRUE)
u_to_m_call_ranges = setdiff(u_to_m_call_ranges, DMRs.gr, ignore.strand=TRUE)
clean_call_ranges = setdiff(clean_call_ranges, DMRs.gr, ignore.strand = TRUE)
clean_variant_call_ranges = setdiff(clean_variant_call_ranges, DMRs.gr, ignore.strand = TRUE)
clean_m_to_u_call_ranges = setdiff(clean_m_to_u_call_ranges, DMRs.gr, ignore.strand = TRUE)
clean_u_to_m_call_ranges = setdiff(clean_u_to_m_call_ranges, DMRs.gr, ignore.strand = TRUE)
clean_single_variant_call_ranges = setdiff(clean_variant_call_ranges, DMRs.gr, ignore.strand = TRUE)
clean_single_m_to_u_call_ranges = setdiff(clean_m_to_u_call_ranges, DMRs.gr, ignore.strand = TRUE)
clean_single_u_to_m_call_ranges = setdiff(clean_u_to_m_call_ranges, DMRs.gr, ignore.strand = TRUE)
M_parent_call_ranges = setdiff(M_parent_call_ranges, DMRs.gr, ignore.strand = TRUE)
U_parent_call_ranges = setdiff(U_parent_call_ranges, DMRs.gr, ignore.strand = TRUE)


# Find the overlaps
olaps_CG_sites = findOverlaps(CG_site_ranges, gene_ranges)
olaps_all_M = findOverlaps(all_M_ranges, gene_ranges)
olaps_all_U = findOverlaps(all_U_ranges, gene_ranges)
olaps_variable = findOverlaps(variant_call_ranges, gene_ranges)
olaps_m_to_u = findOverlaps(m_to_u_call_ranges, gene_ranges)
olaps_u_to_m = findOverlaps(u_to_m_call_ranges, gene_ranges)

tolaps_CG_sites = findOverlaps(CG_site_ranges, transposon_ranges)
tolaps_all_M = findOverlaps(all_M_ranges, transposon_ranges)
tolaps_all_U = findOverlaps(all_U_ranges, transposon_ranges)
tolaps_variable = findOverlaps(variant_call_ranges, transposon_ranges)
tolaps_m_to_u = findOverlaps(m_to_u_call_ranges, transposon_ranges)
tolaps_u_to_m = findOverlaps(u_to_m_call_ranges, transposon_ranges)

### exon intron and UTR olaps fail due to capitalisation or something
eolaps_CG_sites = findOverlaps(CG_site_ranges, exon_ranges)
eolaps_all_M = findOverlaps(all_M_ranges, exon_ranges)
eolaps_all_U = findOverlaps(all_U_ranges, exon_ranges)
eolaps_variable = findOverlaps(variant_call_ranges, exon_ranges)
eolaps_m_to_u = findOverlaps(m_to_u_call_ranges, exon_ranges)
eolaps_u_to_m = findOverlaps(u_to_m_call_ranges, exon_ranges)

iolaps_CG_sites = findOverlaps(CG_site_ranges, intron_ranges)
iolaps_all_M = findOverlaps(all_M_ranges, intron_ranges)
iolaps_all_U = findOverlaps(all_U_ranges, intron_ranges)
iolaps_variable = findOverlaps(variant_call_ranges, intron_ranges)
iolaps_m_to_u = findOverlaps(m_to_u_call_ranges, intron_ranges)
iolaps_u_to_m = findOverlaps(u_to_m_call_ranges, intron_ranges)

olaps5_CG_sites = findOverlaps(CG_site_ranges, UTR5_ranges)
olaps5_all_M = findOverlaps(all_M_ranges, UTR5_ranges)
olaps5_all_U = findOverlaps(all_U_ranges, UTR5_ranges)
olaps5_variable = findOverlaps(variant_call_ranges, UTR5_ranges)
olaps5_m_to_u = findOverlaps(m_to_u_call_ranges, UTR5_ranges)
olaps5_u_to_m = findOverlaps(u_to_m_call_ranges, UTR5_ranges)

olaps3_CG_sites = findOverlaps(CG_site_ranges, UTR3_ranges)
olaps3_all_M = findOverlaps(all_M_ranges, UTR3_ranges)
olaps3_all_U = findOverlaps(all_U_ranges, UTR3_ranges)
olaps3_variable = findOverlaps(variant_call_ranges, UTR3_ranges)
olaps3_m_to_u = findOverlaps(m_to_u_call_ranges, UTR3_ranges)
olaps3_u_to_m = findOverlaps(u_to_m_call_ranges, UTR3_ranges)


# Generate a count of CG sites per gene
gff.genes$CG_sites_count=table(gff.genes[subjectHits(olaps_CG_sites),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$CG_sites_count=ifelse(is.na(gff.genes$CG_sites_count),0,gff.genes$CG_sites_count)

# Generate a count of all_M sites per gene
gff.genes$all_M_count=table(gff.genes[subjectHits(olaps_all_M),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$all_M_count=ifelse(is.na(gff.genes$all_M_count),0,gff.genes$all_M_count)

# Generate a count of all_U sites per gene
gff.genes$all_U_count=table(gff.genes[subjectHits(olaps_all_U),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$all_U_count=ifelse(is.na(gff.genes$all_U_count),0,gff.genes$all_U_count)

# Generate a count of variable sites per gene
gff.genes$variable_count=table(gff.genes[subjectHits(olaps_variable),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$variable_count=ifelse(is.na(gff.genes$variable_count),0,gff.genes$variable_count)

# Generate a count of m->u sites per gene
gff.genes$m_to_u_count=table(gff.genes[subjectHits(olaps_m_to_u),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$m_to_u_count=ifelse(is.na(gff.genes$m_to_u_count),0,gff.genes$m_to_u_count)

# Generate a count of u->m sites per gene
gff.genes$u_to_m_count=table(gff.genes[subjectHits(olaps_u_to_m),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$u_to_m_count=ifelse(is.na(gff.genes$u_to_m_count),0,gff.genes$u_to_m_count)


# Generate a count of CG sites per transposon
gff.transposons$CG_sites_count=table(gff.transposons[subjectHits(tolaps_CG_sites),]$gene_ID)[gff.transposons$gene_ID]
# Replace NA values in variant_count with 0s
gff.transposons$CG_sites_count=ifelse(is.na(gff.transposons$CG_sites_count),0,gff.transposons$CG_sites_count)

# Generate a count of all_M sites per transposon
gff.transposons$all_M_count=table(gff.transposons[subjectHits(tolaps_all_M),]$gene_ID)[gff.transposons$gene_ID]
# Replace NA values in variant_count with 0s
gff.transposons$all_M_count=ifelse(is.na(gff.transposons$all_M_count),0,gff.transposons$all_M_count)

# Generate a count of all_U sites per transposon
gff.transposons$all_U_count=table(gff.transposons[subjectHits(tolaps_all_U),]$gene_ID)[gff.transposons$gene_ID]
# Replace NA values in variant_count with 0s
gff.transposons$all_U_count=ifelse(is.na(gff.transposons$all_U_count),0,gff.transposons$all_U_count)

# Generate a count of variable sites per transposon
gff.transposons$variable_count=table(gff.transposons[subjectHits(tolaps_variable),]$gene_ID)[gff.transposons$gene_ID]
# Replace NA values in variant_count with 0s
gff.transposons$variable_count=ifelse(is.na(gff.transposons$variable_count),0,gff.transposons$variable_count)

# Generate a count of m->u sites per transposon
gff.transposons$m_to_u_count=table(gff.transposons[subjectHits(tolaps_m_to_u),]$gene_ID)[gff.transposons$gene_ID]
# Replace NA values in variant_count with 0s
gff.transposons$m_to_u_count=ifelse(is.na(gff.transposons$m_to_u_count),0,gff.transposons$m_to_u_count)

# Generate a count of u->m sites per transposon
gff.transposons$u_to_m_count=table(gff.transposons[subjectHits(tolaps_u_to_m),]$gene_ID)[gff.transposons$gene_ID]
# Replace NA values in variant_count with 0s
gff.transposons$u_to_m_count=ifelse(is.na(gff.transposons$u_to_m_count),0,gff.transposons$u_to_m_count)


# Find out which genes are really genes, and which are 'heterochromatic genes' as per Zemach et al, 2013 Figure 6A.
### This is just putting total nos sites in the average_methylation column
gff.genes$average_methylation=table(gff.genes[subjectHits(olaps_CG_sites),]$gene_ID)[gff.genes$gene_ID]

# Find out which transposons are really transposons, and which are 'euchromatic transposons' as per my own method.
### This needs more
gff.transposons$average_methylation=table(gff.transposons[subjectHits(tolaps_CG_sites),]$gene_ID)[gff.transposons$gene_ID]


# Create a data frame with gene ID and average methylation for each overlap between a site and a gene model
olaps_average_methylation=data.frame(cbind("gene_ID"=gff.genes[subjectHits(olaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(olaps_CG_sites)]))
tolaps_average_methylation=data.frame(cbind("gene_ID"=gff.transposons[subjectHits(tolaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(tolaps_CG_sites)]))
# Convert average values to numerics
olaps_average_methylation$site_average_methylation=as.numeric(levels(olaps_average_methylation$site_average_methylation)[olaps_average_methylation$site_average_methylation])
tolaps_average_methylation$site_average_methylation=as.numeric(levels(tolaps_average_methylation$site_average_methylation)[tolaps_average_methylation$site_average_methylation])

# Find the mean of sites average methylation levels across each gene, and merge this with the genes table
#x = merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], mean), by="gene_ID", all=TRUE)
#gff.genes$average_methylation = x$site_average_methylation[match(gff.genes$gene_ID, x$gene_ID)]
gff.genes$average_methylation = merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], mean), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.genes$gene_ID, merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], mean), by="gene_ID", all=TRUE)$gene_ID)]
gff.transposons$average_methylation = merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], mean), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.transposons$gene_ID, merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], mean), by="gene_ID", all=TRUE)$gene_ID)]

# Find the variance of sites average methylation levels across each gene, and merge this with the genes table
gff.genes$variance_methylation = merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], var), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.genes$gene_ID, merge(gff.genes, aggregate(. ~ gene_ID, olaps_average_methylation[], var), by="gene_ID", all=TRUE)$gene_ID)]
gff.transposons$variance_methylation = merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], var), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.transposons$gene_ID, merge(gff.transposons, aggregate(. ~ gene_ID, tolaps_average_methylation[], var), by="gene_ID", all=TRUE)$gene_ID)]

# For a more hands-off way to set the threshold, fit a mixture of Gaussians to the dispersion index density, and identify the intersection of the fit models
library(mixtools)
library(rootSolve)
# Fit a mixture of 2 Gaussians to the dispersion index density
#mixmdl = normalmixEM(gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$variance_methylation/gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$average_methylation)

# Fit a mixture of 3 Gaussians to the dispersion index density
mixmdl = normalmixEM(gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$variance_methylation/gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$average_methylation, k=ifelse(meth_context=="CG",3,ifelse(meth_context=="CHG",2,2)))

# Mixture of 2 Gaussians fits better once 'P' CG sites are removed from the picture
#mixmdl = normalmixEM(gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$variance_methylation/gff.genes[!is.na(gff.genes$variance_methylation/gff.genes$average_methylation),]$average_methylation, k=ifelse(meth_context=="CG",2,ifelse(meth_context=="CHG",2,2)))

# Print the characteristics of the fit curves
cat(paste0("Fit lambda values (mixture proportions):\n"))
mixmdl$lambda
# Schmitz data: 0.2853651 0.1365236 0.5781114.  After converting P to U/M with mCG=0.45 cutoff: 0.2381327 0.7618673. Using Fishers p<0.05, Binomial p<0.005: 0.2413691 0.1474132 0.6112177
# Becker data:  0.08638043 0.17313572 0.74048385
# 1001 methylomes:   0.1665908 0.5831132 0.2502960

cat(paste0("Fit mu values (means):\n"))
mixmdl$mu
# Schmitz data: 0.01777965 0.10733632 0.64071733.  After converting P to U/M with mCG=0.45 cutoff: 0.08997885 0.65362097. Using Fishers p<0.05, Binomial p<0.005: 0.01632469 0.10906799 0.64090168
# Becker data:  0.003837026 0.070503558 0.635220579
# 1001 methylomes:   0.001042594 0.295110614 0.016815323

cat(paste0("Fit sigma values (st.devs):\n"))
mixmdl$sigma
# Schmitz data: 0.01006439 0.05852205 0.14541845.  After converting P to U/M with mCG=0.45 cutoff: 0.07725693 0.14292752. Using Fishers p<0.05, Binomial p<0.005: 0.01027012 0.05664957 0.14650136
# Becker data:  0.004954326 0.045214098 0.165311575
# 1001 methylomes:   0.0006754618 0.1620353312 0.0132676875

# Find the intersection of the two larger components in the mixture - here we cut off
mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
  dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }
low_model=ifelse(meth_context=="CG",2,ifelse(meth_context=="CHG",1,1))
#low_model=1
high_model=low_model+1
heterochromatic_gene_dispersion_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=mixmdl$mu[low_model], sd1=mixmdl$sigma[low_model], m2=mixmdl$mu[high_model], sd2=mixmdl$sigma[high_model], p1=mixmdl$lambda[low_model], p2=mixmdl$lambda[high_model])
heterochromatic_gene_dispersion_cutoff
# Schmitz CG data: 0.2517447, Becker CG data: 0.1896873
# Schmitz CG data after 'P' sites removed: 0.2768782, Becker:
# Schmitz data with Fishers p<0.05, Binomial p<0.005: 0.2499017
# Becker CHG data: 0.9766602
# 1001 methylomes CG data: 0.04848119
plot(mixmdl, whichplots=2)

# Plot gene methylation level vs. gene length, colour by chromosome
heterochromatic_gene_coverage_cutoff=0.6
heterochromatic_gene_coverage_minimum=0.2

gff.genes$m_class = ifelse(gff.genes$CG_sites_count<=1,"CG-poor",ifelse(gff.genes$all_M_count+gff.genes$all_U_count+gff.genes$variable_count==0,"Unknown",ifelse(gff.genes$average_methylation==0,"Unmethylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation>heterochromatic_gene_dispersion_cutoff),"Gene-body Methylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.genes$average_methylation<heterochromatic_gene_coverage_minimum),"Unmethylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.genes$average_methylation>=heterochromatic_gene_coverage_minimum),"Heterochromatic","Really Unknown"))))))

gff.transposons$m_class = ifelse(gff.transposons$CG_sites_count<=1,"CG-poor",ifelse(gff.transposons$all_M_count+gff.transposons$all_U_count+gff.transposons$variable_count==0,"Unknown",ifelse(gff.transposons$average_methylation==0,"Unmethylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation>heterochromatic_gene_dispersion_cutoff),"Gene-body Methylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.transposons$average_methylation<heterochromatic_gene_coverage_minimum),"Unmethylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.transposons$average_methylation>=heterochromatic_gene_coverage_minimum),"Heterochromatic","Really Unknown"))))))



# Make a data frame with each gene being a row and columns showing the methylation level of that gene in each sample
methylation_per_gene_per_sample = data.frame(gff.genes$gene_ID)
colnames(methylation_per_gene_per_sample)="gene_ID"
sample_no = 0
for(sample_name in valid_samples) {
  sample_no = sample_no + 1
  if(opt$verbose) {cat(paste0("Summing methylation levels in sample ",sample_name,"\n"))}
  sample_average_methylation = matrix(0, nrow=nrow(all_samples_meth_status), ncol=2)
  sample_average_methylation[,1]=ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse(all_samples_meth_status[,sample_name]=="U",0,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",0.25,0))))
  sample_average_methylation[,2]=ifelse(is.na(all_samples_meth_status[,sample_name]),0,ifelse(all_samples_meth_status[,sample_name]=="U",1,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",1,0))))

  sample_olaps_average_methylation=data.frame(cbind("gene_ID"=gff.genes[subjectHits(olaps_CG_sites),"gene_ID"],"site_average_methylation"=(sample_average_methylation[,1]/sample_average_methylation[,2])[queryHits(olaps_CG_sites)]))
  sample_olaps_average_methylation$site_average_methylation = as.numeric(levels(sample_olaps_average_methylation$site_average_methylation)[sample_olaps_average_methylation$site_average_methylation])
  #methylation_per_gene_per_sample = merge(methylation_per_gene_per_sample, sample_name=merge(gff.genes, aggregate(. ~ gene_ID, sample_olaps_average_methylation[], mean), by="gene_ID", all=TRUE)$site_average_methylation[match(gff.genes$gene_ID, merge(gff.genes, aggregate(. ~ gene_ID, sample_olaps_average_methylation[], mean), by="gene_ID", all=TRUE)$gene_ID)], by="gene_ID", all=TRUE)
  methylation_per_gene_per_sample = merge(methylation_per_gene_per_sample, aggregate(. ~ gene_ID, sample_olaps_average_methylation[], mean), by="gene_ID", all=TRUE)
  colnames(methylation_per_gene_per_sample)[sample_no + 1] = sample_name
}
write.table(methylation_per_gene_per_sample, file="mean_mCG_per_gene_per_sample.txt", sep="\t", quote=FALSE, rownames=FALSE, colnames=FALSE)
# methylation_per_gene_per_sample = read.delim(file="mean_mCG_per_gene_per_sample.txt")

row.names(methylation_per_gene_per_sample)=methylation_per_gene_per_sample[,1]

#library(gplots)
#heatmap.2(data.matrix(methylation_per_gene_per_sample[1:100,2:101]), col=redgreen(75))

library(pheatmap)
library("RColorBrewer")
#pheatmap(data.matrix(na.omit(methylation_per_gene_per_sample)[1:1010,2:ncol(methylation_per_gene_per_sample)]), color=brewer.pal(9,"Blues"))
pheatmap(data.matrix(methylation_per_gene_per_sample[as.logical((rowSums(is.na(methylation_per_gene_per_sample))-ncol(methylation_per_gene_per_sample)+1)),][1:500,2:ncol(methylation_per_gene_per_sample)]), color=brewer.pal(9,"Blues"))

pca = prcomp(na.omit(data.matrix(methylation_per_gene_per_sample[methylation_per_gene_per_sample$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,][as.logical((rowSums(is.na(methylation_per_gene_per_sample[methylation_per_gene_per_sample$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,]))-ncol(methylation_per_gene_per_sample)+1)),][,2:ncol(methylation_per_gene_per_sample)])), center=TRUE, scale=TRUE)
library(scatterplot3d)
scatterplot3d(pca$x[,1:3], pch=20, color="blue")





gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))
gene_body_loci.gr=makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)

unmethylated_gene_loci.gr=makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)




segmentation_model.gr = readRDS(file = "../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_segmentation_model_draft3.rds")

# Capitalise chromosome names in segmentation_model.gr
### Commented out - we are not using capitals any more
#levels(segmentation_model.gr@seqnames)=toupper(levels(segmentation_model.gr@seqnames))
#segmentation_model.gr@seqinfo@seqnames=levels(segmentation_model.gr@seqnames)

#levels(variant_call_ranges@seqnames)=toupper(levels(variant_call_ranges@seqnames))
#variant_call_ranges@seqinfo@seqnames=levels(variant_call_ranges@seqnames)


GBM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]
TEM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="TEM"]
UMR_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="UMR"]
length(GBM_segments.gr)
length(TEM_segments.gr)
length(UMR_segments.gr)




# Three different versions of target_ranges:

# Identify GenomicRegions covering GBM genes which don't have a TE segment overlapping
target_ranges = setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])

### should we exclude gene_body_loci.gr if they only have 1 or 2 GBM sites?

# Identify GenomicRegions covering 'body' of such GBM genes (excluding 2kb 5' and 1kb 3')
#target_ranges = 

# Identify GenomicRegions covering GBM range of such GBM genes (from first to last methylated base in gen 3)
#target_ranges = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])))]



  
# Find U sites in target_ranges 
target_U_sites = U_parent_call_ranges[queryHits(findOverlaps(U_parent_call_ranges, target_ranges))]

# Find M sites in target_ranges
target_M_sites = M_parent_call_ranges[queryHits(findOverlaps(M_parent_call_ranges, target_ranges))]

#z=cbind(as.data.frame(target_U_sites), dist_to_nearest_M = rep(NA,length(target_U_sites)), u_to_m_count = rep(NA,length(target_U_sites)))
zU=cbind(as.data.frame(target_U_sites), dist_to_nearest_M = rep(NA,length(target_U_sites)), sites_to_nearest_M = rep(NA,length(target_U_sites)), dist_to_nearest_U = rep(NA,length(target_U_sites)), sites_to_nearest_U = rep(NA,length(target_U_sites)), u_to_m_count = rep(NA,length(target_U_sites)))



# Find sites in target_ranges
target_U_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_ranges))]



num_metrics = 8
methylation_summary = matrix(, nrow=length(valid_samples), ncol=num_metrics + 2)
methylation_summary[,1] = valid_samples
methylation_summary[,2] = as.character(sample_list[match(valid_samples, sample_list$SRA.Accession), "Title"])
	
for(sample_no in 1:length(valid_samples)) {
  cat(paste0("Processing ", valid_samples[sample_no],"\n"))

  this_sample = all_samples_meth_status[, valid_samples[sample_no]]
  
  # work out whole genome methylation for sample
  methylation_summary[sample_no, 3] = sum(this_sample %in% c("M"))/sum(this_sample %in% c("M","U"))

  # work out whole genome coverage for sample
  methylation_summary[sample_no, 4] = sum(this_sample %in% c("M","U"))/sum(this_sample %in% c("M","U","P","I"))
  
  # work out proportion of missing data due to mutations
  methylation_summary[sample_no, 9] = sum(this_sample %in% c("X"))/sum(this_sample %in% c("M","U","P","I","X"))
  # proportion of sites with identified mutations, ignoring NA (zero coverage?)

  # work out proportion of missing data due to NAs (zero coverage? missing due to long indels?)
  methylation_summary[sample_no, 10] = sum(is.na(this_sample))/length(this_sample)
  # proportion of sites with identified mutations, ignoring NA (zero coverage?)


  # work out mean GBM gene methylation for sample
  target_ranges = gene_body_loci.gr
  target_sample = this_sample[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
  methylation_summary[sample_no, 5] = sum(target_sample %in% c("M"))/sum(target_sample %in% c("M","U"))
  
  # work out coverage of GBM genes
  #target_ranges = gene_body_loci.gr
  #target_sample = this_sample[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
  methylation_summary[sample_no, 6] = sum(target_sample %in% c("M","U"))/sum(target_sample %in% c("M","U","P","I"))
  # proportion of clean calls vs. all sites with any call except "X" (mutation) or NA (zero coverage?)
   
  # work out mean GBM in Col-0 GBM segments per sample
  target_ranges = GBM_segments.gr
  target_sample = this_sample[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
  methylation_summary[sample_no, 7] = sum(target_sample %in% c("M"))/sum(target_sample %in% c("M","U"))
  
  # work out coverage of Col-0 GBM segments
  #target_ranges = GBM_segments.gr
  #target_sample = this_sample[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
  methylation_summary[sample_no, 8] = sum(target_sample %in% c("M","U"))/sum(target_sample %in% c("M","U","P","I"))
  # proportion of clean calls vs. all sites with any call except "X" (mutation) or NA (zero coverage?)
   
   
   
   
  # work out no. genes with GBM for sample

  # work out mean gbM as average of gbM at individual CG sites

}

methylation_summary = as.data.frame(methylation_summary)

colnames(methylation_summary) = c("SRA.Accession", "Title", "genome_mCG", "genome_coverage", "gbm_gene_mCG", "gbm_gene_coverage", "Col_0_gbM_segment_mCG", "Col_0_gbM_segment_coverage", "missing_mutation", "missing_NA")
colnames(no_mutations) = c("SRA.Accession", "CG_site_mutations")
# Merge in the number of mutated sites
methylation_summary = merge(methylation_summary, no_mutations, by="SRA.Accession", all=TRUE)

#colnames(methylation_summary) = c("SRA.Accession", "Title", "genome_mCG", "genome_coverage", "gbm_gene_mCG", "gbm_gene_coverage", "Col_0_gbM_segment_mCG", "Col_0_gbM_segment_coverage", "missing_mutation", "missing_NA", "CG_site_mutations")
for(i in 3:(num_metrics + 3)) {
  methylation_summary[,i] = as.numeric(levels(methylation_summary[,i])[methylation_summary[,i]])
}

write.table(merge(methylation_summary, sample_info_merged[sample_info_merged$SRA.Accession %in% valid_samples,], by="SRA.Accession", all=TRUE), file="1001_methylomes_summary.txt", sep="\t")

pdf(paste0(project_id,"_",meth_context,"_",sample_id,"_summary_plots_CG.pdf"))

# plot missing data:
print(ggplot(methylation_summary[,]) + geom_point(aes(x=missing_NA, y=missing_mutation)))
# missing_mutation rises exponentially with missing_NA
# missing_NA>0.5 would seem a reasonable cutoff.  missing_mutation rises very rapidly above 0.125 (missing_NA=0.875)

# Compare mCG with missing data
print(ggplot(methylation_summary[,]) + geom_point(aes(x=missing_NA, y=genome_mCG)))
# All three accessions with global mCG>0.4 have missing_NA>0.875

# plot mCG vs coverage in each context
print(ggplot(methylation_summary[,]) + geom_point(aes(x=genome_coverage, y=genome_mCG)))
print(ggplot(methylation_summary[,]) + geom_point(aes(x=gbm_gene_coverage, y=gbm_gene_mCG)))
print(ggplot(methylation_summary[,]) + geom_point(aes(x=Col_0_gbM_segment_coverage, y=Col_0_gbM_segment_mCG)))
# these look uncorrelated for the first 450 samples
# The sites with the highest genome_mCG (and genic) also have high coverage, so are not artificially being reported as hyper-methylated due to lack of coverage
# The sites with the lowest gbM segment mCG also have very low coverage

# plot coverage in genes vs. whole genome
print(ggplot(methylation_summary[,]) + geom_point(aes(x=genome_coverage, y=gbm_gene_coverage)))
# Not correlated, some outliers in both directions with high coverage in genome but low in genes and vice versa
 
# plot coverage in segments vs. genes
print(ggplot(methylation_summary[,]) + geom_point(aes(x=gbm_gene_coverage, y=Col_0_gbM_segment_coverage)))
# very strongly correlated

# compare gene coverage with missing data
print(ggplot(methylation_summary[,]) + geom_point(aes(x=missing_NA, y=gbm_gene_coverage)))
# accessions with missing data are largely separate from accessions with lack of gbm gene coverage
# former is likely due to mutation, latter maybe more likely coverage
# Probably need filters on both: missing_NA<0.5 & gbm_gene_coverage>0.5:

nrow(methylation_summary[methylation_summary$gbm_gene_coverage>0.5 & methylation_summary$missing_NA<0.5,])
#[1] 417
nrow(methylation_summary)
#[1] 453

minimum_gbm_gene_coverage = 0.5
maximum_missing_NA = 0.5
# compare mCG between contexts
#print(ggplot(methylation_summary[,]) + geom_point(aes(x=genome_mCG, y=gbm_gene_mCG)))
print(ggplot(cbind(methylation_summary[,], "coverage_gt_0.5"=methylation_summary$gbm_gene_coverage>=minimum_gbm_gene_coverage & methylation_summary$missing_NA<=maximum_missing_NA)) + geom_point(aes(x=genome_mCG, y=gbm_gene_mCG, colour=coverage_gt_0.5)))

#print(ggplot(methylation_summary[,]) + geom_point(aes(x=gbm_gene_mCG, y=Col_0_gbM_segment_mCG)))
print(ggplot(cbind(methylation_summary[,], "coverage_gt_0.5"=methylation_summary$gbm_gene_coverage>=minimum_gbm_gene_coverage & methylation_summary$missing_NA<=maximum_missing_NA)) + geom_point(aes(x=gbm_gene_mCG, y=Col_0_gbM_segment_mCG, colour=coverage_gt_0.5)))

minimum_gbm_gene_coverage = 0.8
maximum_missing_NA = 0.2
# compare mCG between contexts
#print(ggplot(methylation_summary[,]) + geom_point(aes(x=genome_mCG, y=gbm_gene_mCG)))
print(ggplot(cbind(methylation_summary[,], "coverage_gt_0.8"=methylation_summary$gbm_gene_coverage>=minimum_gbm_gene_coverage & methylation_summary$missing_NA<=maximum_missing_NA)) + geom_point(aes(x=genome_mCG, y=gbm_gene_mCG, colour=coverage_gt_0.8)))

#print(ggplot(methylation_summary[,]) + geom_point(aes(x=gbm_gene_mCG, y=Col_0_gbM_segment_mCG)))
print(ggplot(cbind(methylation_summary[,], "coverage_gt_0.8"=methylation_summary$gbm_gene_coverage>=minimum_gbm_gene_coverage & methylation_summary$missing_NA<=maximum_missing_NA)) + geom_point(aes(x=gbm_gene_mCG, y=Col_0_gbM_segment_mCG, colour=coverage_gt_0.8)))

# For accessions with enough coverage (let's say 0.8 for now) do the ones with abnormally high mCG in Col-0 segments show higher alignment rates due to being closer to Col-0?
print(ggplot(cbind(methylation_summary[methylation_summary$gbm_gene_coverage>=minimum_gbm_gene_coverage & methylation_summary$missing_NA<=maximum_missing_NA,], "like_Col_0"=methylation_summary[methylation_summary$gbm_gene_coverage>=minimum_gbm_gene_coverage & methylation_summary$missing_NA<=maximum_missing_NA,]$missing_NA<0.03)) + geom_point(aes(x=gbm_gene_mCG, y=Col_0_gbM_segment_mCG, colour=like_Col_0)))
# Accessions with <3% lost coverage due to genetic distance tend to show significantly higher Col0_gbM_segment_mCG relative to gbm_gene_mCG as expected since they are closer to Col-0


dev.off()


# bodge to incorporate Schmitz data as a check
all_samples_meth_status_col0 = readRDS(file="../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_all_samples_meth_status.rds")# Define a convenience function to convert chromosome names to mixed-case (e.g. Chr9, ChrC)
mixedCaseChr <- function(s, strict = FALSE) {
    paste0(toupper(substr(s,1,1)),tolower(substr(s,2,3)),toupper(substr(s,4,4)))
}
all_samples_meth_status_col0$Chromosome = mixedCaseChr(all_samples_meth_status_col0$Chromosome)
x=merge(all_samples_meth_status, all_samples_meth_status_col0, by=c("Chromosome", "Locus", "Strand"), all.x=TRUE, sort=TRUE)
CG_site_ranges1 =makeGRangesFromDataFrame(df = x, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")

rm(all_samples_meth_status_col0)

for(sample_no in 477:499) {
  cat(paste0("Processing ", sample_no,"\n"))

  this_sample = x[, sample_no]
  
  # work out whole genome methylation for sample
  cat(sum(this_sample %in% c("M"))/sum(this_sample %in% c("M","U")))
  
  # work out mean GBM gene methylation for sample
  target_ranges = gene_body_loci.gr
  target_sample = this_sample[queryHits(findOverlaps(CG_site_ranges1, target_ranges))]
  cat(sum(target_sample %in% c("M"))/sum(target_sample %in% c("M","U")))
  
  # work out mean GBM in Col-0 GBM segments per sample
  target_ranges = GBM_segments.gr
  target_sample = this_sample[queryHits(findOverlaps(CG_site_ranges1, target_ranges))]
  cat(sum(target_sample %in% c("M"))/sum(target_sample %in% c("M","U")))
  
  # work out no. genes with GBM for sample

}

rm(x)
rm(CG_site_ranges1)
  




'


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



