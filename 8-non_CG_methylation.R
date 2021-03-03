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

library(data.table)
library(stringr)
library(ggplot2)
library(plyr)
#library(ggpubr)


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
sample_list1 = read.table(paste0("../0-raw_data/","GSE43857.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=43857
fastq_list1 = read.table(paste0("../0-raw_data/","ENA_SRA065807.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/SRA065807

# Swedish Arabidopsis accessions from separate study (284 of)
sample_list2 = read.table(paste0("../0-raw_data/","GSE54292.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=54292
fastq_list2 = read.table(paste0("../0-raw_data/","ENA_SRP035593.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/PRJNA236110

sample_list = rbind.data.frame(sample_list1, sample_list2)
fastq_list = rbind.data.frame(fastq_list1, fastq_list2)

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

quit()

# get the list of runs for the sample
# concatenate the methylation data for the runs
all_samples_ct=NA
run_no = 0
for (this_run in fastq_list[fastq_list$experiment_accession==sample_id,"run_accession"]) {
	run_no = run_no + 1
	if(opt$verbose) {cat(paste0("Loading ",this_run,"\n"))}
	x=read.table(file=paste0(source_dir, this_run, "_TAIR10.cg.methratio2"), header=TRUE)

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

	all_samples_ct[,paste0(sample_id,"_C")] = all_samples_ct[,paste0(sample_id,"_C")] + all_samples_ct[,paste0(this_run,"_C")]
	all_samples_ct[,paste0(sample_id,"_T")] = all_samples_ct[,paste0(sample_id,"_T")] + all_samples_ct[,paste0(this_run,"_T")]
}

colnames(all_samples_ct)[1]="Chromosome"
colnames(all_samples_ct)[2]="Locus"
colnames(all_samples_ct)[3]="Strand"

all_samples_ct = all_samples_ct[,c(1, 2, 3, 4 + run_no*2, 5 + run_no*2)]
all_samples_ct = data.table(all_samples_ct,key=c("Chromosome", "Locus","Strand"))

	
# Dump the data table to file for a convenience cache
saveRDS(all_samples_ct, file=paste0(project_id,"_",meth_context,"_",sample_id,"_all_samples_ct.rds"))

# call methylation

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
#x=read.table(file="../../0-raw_data/GSM2323848_WT.spm.CG.tair10 (1).gff")



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
	
  
	
