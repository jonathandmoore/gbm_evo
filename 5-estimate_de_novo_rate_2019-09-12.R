####EARLY SETUP STUFF COPIED FROM 5-methylation_calls_2018-06-26.R

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
  #project_id="test"
  project_id="SRA035939"  # Schmitz et al, 2011
  #project_id="PRJEB2678"  # Becker et al, 2011
} else {
  project_id=opt$project
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


# Platform-specific stuff:
#  Base of path to project
#  Where to find bioconductor (packages need to be pre-installed in R using singularity if running on the cluster)
on_cluster = FALSE
pathroot = ""
if (.Platform$OS.type=="windows") {
  pathroot="D:"
  #source("http://bioconductor.org/biocLite.R")
} else if (.Platform$OS.type=="unix") {
  if(substr(getwd(),1,20)=="/jic/scratch/groups/") {
    # assume we are running on cluster
    pathroot="/jic/scratch/groups/Daniel-Zilberman"
	on_cluster = TRUE
} else {
    # assume we are running in Virtualbox with shared folder /media/sf_D_DRIVE
    pathroot="/media/sf_D_DRIVE"
    #source("http://bioconductor.org/biocLite.R")
  }
}
pathroot="X:/Daniel-Zilberman/"


if (!on_cluster) {
  # Assume that if we are on the cluster, then we are running in a VM with all relevant packages pre-installed
  install.packages("BiocManager")
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
#if(!require(ggpubr)){
#  install.packages("ggpubr")
#  library(ggpubr)
#}

#### This part sets up libraries and contextual information and metadata about the project

#library(data.table)
#library(stringr)
#library(ggplot2)
library(plyr)
#library(ggpubr)


# Source directory containing alignments from bs_sequel
source_dir = paste0(pathroot,"/Projects/Jay-SMP_variation/3-alignments/",project_id,"/")

# Working directory where outputs will be directed to
setwd(dir = paste0(pathroot,"/Projects/Jay-SMP_variation/5-analysis/",project_id,"/"))

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
sample_metadata=read.delim(paste0(source_dir,project_id,"_metadata.tsv"), header=TRUE, sep="\t")
rownames(sample_metadata)=sample_metadata$Identifier

# Read the Schmitz et al, 2011 DMRs data, so this can be used in masking where needed
Schmitz_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_CG_DMRs.txt")
Schmitz_non_CG_DMRs = paste0(pathroot,"/Literature/Schmitz_et_al_2011_non-CG_DMRs.txt")


#if (action=="analyse") {
  #### Reload the data objects developed in earlier actions
  #all_samples_ct <- readRDS(paste0(project_id,"_",meth_context,"_all_samples_ct.rds")) # Commented out as not needed any more?
  sample_info <- readRDS(paste0(project_id,"_",meth_context,"_sample_info.rds"))
  coverage_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_coverage_data.rds")))
  meth_data <- data.table(readRDS(paste0(project_id,"_",meth_context,"_meth_data.rds")))
  meth_data_binom <- data.table(readRDS(paste0(project_id,"_",meth_context,"_meth_data_binom.rds")))
  meth_data_fisher <- data.table(readRDS(paste0(project_id,"_",meth_context,"_meth_data_fisher_2.rds")))
  # We reload the following as data tables, as the logic tests for all_m, all_u etc fail with data.tables
  all_samples_meth_status <- data.frame(readRDS(paste0(project_id,"_",meth_context,"_all_samples_meth_status.rds")))
  all_samples_p_value <- data.frame(readRDS(paste0(project_id,"_",meth_context,"_all_samples_p_value.rds")))
  if (meth_context=="CG") {
    # nothing special so far?
  }
#} else {
#  all_samples_meth_status = data.frame(all_samples_meth_status)
#  all_samples_p_value = data.frame(all_samples_p_value)
#  if (meth_context=="CG") {
#    # nothing special so far?
#  }
#}


# What is the highest p-value from the binomial test for a site we identified as "Partial"?
cat(paste0("Max binomial p-value for any 'P' site:",max(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="P",3]),"\n"))
# What is the highest p-value from the binomial test for a site we identified as "Methylated"?
cat(paste0("Max binomial p-value for any 'M' site:",max(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="M",3]),"\n"))
# What is the lowest p-value from the binomial test for a site we identified as "Indeterminate"?
cat(paste0("Min binomial p-value for any 'I' site:",min(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="I",3],na.rm=TRUE),"\n"))

# What is the lowest p-value from the binomial test for a site we identified as "Indeterminate" where overall coverage >=10?
cat(paste0("Min binomial p-value for any 'I' site with coverage>=10:",min(cbind(coverage_data,meth_data,meth_data_binom)[meth_data$meth_status=="I" & ((coverage_data$Cov_C+coverage_data$Cov_T)>=10),9],na.rm=TRUE),"\n"))

# What is the lowest p-value from the binomial test for a site we identified as "Unmethylated"?
cat(paste0("Min binomial p-value for any 'U' site:",min(cbind(meth_data,meth_data_binom)[meth_data$meth_status=="U",3],na.rm=TRUE),"\n"))

# What is the lowest p-value from the binomial test for a site we identified as "Unmethylated" where overall coverage >=10?
cat(paste0("Min binomial p-value for any 'U' site with coverage>=10:",min(cbind(coverage_data,meth_data,meth_data_binom)[meth_data$meth_status=="U" & ((coverage_data$Cov_C+coverage_data$Cov_T)>=10),9],na.rm=TRUE),"\n"))


# Merge the sample metadata with the sample information
sample_info=cbind(sample_info,sample_metadata)

# Identify site with either an "M" or "U" call in all samples
# Identify sites among these without identical calls in all samples

# Some samples may need to be excluded from this analysis
# For Becker and Schmitz data sets we want to exclude Line 69 as it is a methylation mutant
excluded_samples = c()
valid_samples=c()
if (project_id=="SRA035939") {
  excluded_samples = c("SRR342380","SRR342390")
} else if (project_id=="PRJEB2678") {
  excluded_samples = c("ERR046563","ERR046564")
}
valid_samples=rownames(sample_info)[!(rownames(sample_info) %in% excluded_samples)]
valid_samples_meth_status=all_samples_meth_status[,!names(all_samples_meth_status) %in% excluded_samples]


# Find out which sites are generally always methylated, and which unmethylated, across the genome

# This part sums up methylation calls across lines to give an 'average' methylation level per site irrespective of line
# Partial calls are assigned an arbitrary methylation level of 0.25
average_methylation = matrix(0, nrow=nrow(all_samples_meth_status), ncol=2)
for(sample_name in valid_samples) {
  if(opt$verbose) {cat(paste0("Summing methylation levels in sample ",sample_name,"\n"))}
  average_methylation[,1]=average_methylation[,1]+ifelse(all_samples_meth_status[,sample_name]=="U",0,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",0.25,0)))
  average_methylation[,2]=average_methylation[,2]+ifelse(all_samples_meth_status[,sample_name]=="U",1,ifelse(all_samples_meth_status[,sample_name]=="M",1,ifelse(all_samples_meth_status[,sample_name]=="P",1,0)))
}


# Accumulate methylation calls which are concordant among replicates per line_generation

line_generation_no = 0
for(this_line in unique(sample_info$Line)) {
  for(this_generation in unique(sample_info[sample_info$Line==this_line,]$Generation)) {
    sample_count=0
    reps_concordant_meth_status=NULL
    for(sample_name in rownames(sample_info[(!rownames(sample_info) %in% excluded_samples) & sample_info$Line==this_line & sample_info$Generation==this_generation,])) {
      cat(paste(this_line, this_generation, sample_name, "\n"))
      sample_count=sample_count+1
      if (sample_count == 1) {
        # This is the first replicate of replicate set
        line_generation_no = line_generation_no + 1
        # Cache the methylation status and move on
        reps_concordant_meth_status=all_samples_meth_status[,sample_name]
      } else {
        # We are doing a second or higher rep - check for concordance between this sample and the consensus of the previous reps_concordant_meth_status
        # If the current sample disagrees with the previous ones, record a status of "D" for disagreement/discordant for the locus, unless status is opposite (M/U) or "O" already found, in which case "O"
        reps_concordant_meth_status=as.factor(ifelse(levels(reps_concordant_meth_status)[reps_concordant_meth_status]==all_samples_meth_status[,sample_name],levels(reps_concordant_meth_status)[reps_concordant_meth_status],ifelse(((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="M") & (all_samples_meth_status[,sample_name]=="U")) | ((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="U") & (all_samples_meth_status[,sample_name]=="M")) | ((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="O") & ((all_samples_meth_status[,sample_name]=="M") | (all_samples_meth_status[,sample_name]=="U"))),"O","D")))
        ### Could do something here with the Fisher method to combine p-values from the concordant reps
      }
    }
    # All reps done for this line and generation, report and cache the results
    if (sample_count>0) {
      if (line_generation_no==1) {
        all_reps_meth_status=cbind(all_samples_meth_status[,c("Chromosome","Locus","Strand")],reps_concordant_meth_status)
        #all_reps_p_value=sample_p_value
      } else {
        all_reps_meth_status=cbind(all_reps_meth_status, reps_concordant_meth_status)
        #all_reps_p_value=cbind(all_samples_p_value, sample_p_value, by=c("Chromosome","Locus","Strand"), all=TRUE)
      }
      colnames(all_reps_meth_status)[line_generation_no+3]=paste0("Line",this_line,"Gen",this_generation)
      cat(paste0(sum(all_reps_meth_status[,line_generation_no+3]=="O")," sites have opposite calls (M/U) among reps.\n"))
      # Tidy up the temporary objects
      rm(reps_concordant_meth_status, sample_count)
    }
  }
}

saveRDS(all_reps_meth_status, file=paste0(project_id,"_",meth_context,"_all_reps_meth_status.rds"))

# Summarise the results of the rep concordance analysis
line_gen_no = 0
rep_concordance_summary=matrix(, nrow=ncol(all_reps_meth_status)-3, ncol=6)
rownames(rep_concordance_summary)=colnames(all_reps_meth_status[,4:ncol(all_reps_meth_status)])
for(line_gen in rownames(rep_concordance_summary)) {
  line_gen_no = line_gen_no + 1
  cat(paste0("Summarising methylation call basis for ",line_gen,"\n"))
  status_counts = count(all_reps_meth_status[,line_gen])
  # Needs a hack for when "P" is zero (nrow should be 6 otherwise):
  if(nrow(status_counts)==5) {
    status_counts$x = as.character(status_counts$x)
    status_counts[6,]=status_counts[5,]
    status_counts[5,]=c("P",0)
    status_counts$x = as.factor(status_counts$x)
    status_counts$freq=as.numeric(status_counts$freq)
  }
  rep_concordance_summary[line_gen_no,]=status_counts$freq
  colnames(rep_concordance_summary)=status_counts$x
}
# output the summary of numbers of concordant and discordant sites
rep_concordance_summary
# output the proportion of sites in each line with a 'clean' M/U call
(rep_concordance_summary[,3]+rep_concordance_summary[,6])/sum(rep_concordance_summary[,1:6])*nrow(rep_concordance_summary)

# This part looks for M/U calls across all samples, cognisant of rep structure
line_gen_no=0
first_line_gen = TRUE
for(line_gen in rownames(rep_concordance_summary)) {
  line_gen_no = line_gen_no + 1
  if(opt$verbose) {cat(paste0("Checking variation in methylation calls concordant between reps: ",line_gen,"\n"))}
  if (first_line_gen) {
    all_m_or_u=ifelse(((all_reps_meth_status[,line_gen_no+3]=="U") | (all_reps_meth_status[,line_gen_no+3]=="M")),TRUE,FALSE)
    all_m=ifelse(all_reps_meth_status[,line_gen_no+3]=="M",TRUE,FALSE)
    all_u=ifelse(all_reps_meth_status[,line_gen_no+3]=="U",TRUE,FALSE)
    first_line_gen=FALSE
  } else {
    all_m_or_u=all_m_or_u & ifelse(((all_reps_meth_status[,line_gen_no+3]=="U") | (all_reps_meth_status[,line_gen_no+3]=="M")),TRUE,FALSE)
    all_m=all_m & ifelse(all_reps_meth_status[,line_gen_no+3]=="M",TRUE,FALSE)
    all_u=all_u & ifelse(all_reps_meth_status[,line_gen_no+3]=="U",TRUE,FALSE)
  }
}

# This version looks for M/U calls across all samples, independent of rep structure
sample_no=0
first_sample = TRUE
for(sample_name in valid_samples) {
  sample_no = sample_no+1
  if(opt$verbose) {cat(paste0("Checking methylation calls in all-samples per-site table: ",sample_name,"\n"))}
  if (first_sample) {
    all_m_or_u_regardless=ifelse(((all_samples_meth_status[,sample_no+3]=="U") | (all_samples_meth_status[,sample_no+3]=="M")),TRUE,FALSE)
    #all_m=ifelse(all_samples_meth_status[,sample_no+3]=="M",TRUE,FALSE)
    #all_u=ifelse(all_samples_meth_status[,sample_no+3]=="U",TRUE,FALSE)
    first_sample=FALSE
  } else {
    all_m_or_u_regardless=all_m_or_u_regardless & ifelse(((all_samples_meth_status[,sample_no+3]=="U") | (all_samples_meth_status[,sample_no+3]=="M")),TRUE,FALSE)
    #all_m=all_m & ifelse(all_samples_meth_status[,sample_no+3]=="M",TRUE,FALSE)
    #all_u=all_u & ifelse(all_samples_meth_status[,sample_no+3]=="U",TRUE,FALSE)
  }
}

### What proportion of sites which vary between reps also have a clean M/U call in all samples?
# This is our strict version of identifying variable sites.  It requires that all samples have sufficient coverage at the site, that there is agreement among reps in all lines, and that there are no differences between lines at generation 3.  One effect of this strictness is to filter out all highly labile sites.

# This part looks for variable calls which are invariable at generation 3
gen3=c()
if (project_id=="SRA035939") {
  gen3=c("Line1Gen3","Line12Gen3","Line19Gen3")
  gen3_30=c("Line29Gen30","Line49Gen30","Line59Gen30","Line119Gen30")
} else if(project_id=="PRJEB2678") {
  gen3=c("Line4Gen3","Line8Gen3")
}

line_gen_no=0
first_line_gen = TRUE
for(line_gen in gen3) {
  line_gen_no = line_gen_no + 1
  if(opt$verbose) {cat(paste0("Checking variation in methylation calls concordant within Gen 3: ",line_gen,"\n"))}
  if (first_line_gen) {
    #all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    all_m_3=ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
    all_u_3=ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
    first_line_gen=FALSE
  } else {
    #all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    all_m_3=all_m_3 & ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
    all_u_3=all_u_3 & ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
  }
}

# This part looks for variable calls which are invariable at generation 3 and 30/1 
### We should redo this analysis in a more sensitive way - check each generation 31 line for changes compared to the corresponding generation 30 line, in concert with reps, but ignoring coverage and calls in generation 3 and in all other generation 30/31 lines.  The current analysis does not really add anything useful.

gen3_30=c()
if (project_id=="SRA035939") {
  # analysis makes no sense for this project, as there is no generation 31
} else if(project_id=="PRJEB2678") {
  gen3_30=c("Line4Gen3","Line8Gen3","Line29Gen31","Line39Gen31","Line49Gen31","Line59Gen31","Line79Gen31","Line89Gen31","Line99Gen31","Line109Gen31","Line119Gen31")
}

line_gen_no=0
first_line_gen = TRUE
for(line_gen in gen3_30) {
  line_gen_no = line_gen_no + 1
  if(opt$verbose) {cat(paste0("Checking variation in methylation calls concordant within Gen 3 and 30: ",line_gen,"\n"))}
  if (first_line_gen) {
    #all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    all_m_3_30=ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
    all_u_3_30=ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
    first_line_gen=FALSE
  } else {
    #all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    all_m_3_30=all_m_3_30 & ifelse(all_reps_meth_status[,line_gen]=="M",TRUE,FALSE)
    all_u_3_30=all_u_3_30 & ifelse(all_reps_meth_status[,line_gen]=="U",TRUE,FALSE)
  }
}

# New more sensitive approach to calling variable sites:
# Find a consensus call for parental lines: if at least two parental lines have calls concordant with reps, and they agree with each other, this is the consensus parental call

gen3=c()
gen3_30=c()
if (project_id=="SRA035939") {
  gen3=c("Line1Gen3","Line12Gen3","Line19Gen3")
  gen3_30=c("Line29Gen30","Line49Gen30","Line59Gen30","Line119Gen30")
} else if(project_id=="PRJEB2678") {
  gen3=c("Line4Gen3","Line8Gen3")
  gen3_30=c("Line29Gen31","Line39Gen31","Line49Gen31","Line59Gen31","Line79Gen31","Line89Gen31","Line99Gen31","Line109Gen31","Line119Gen31")
}

# Identify parental consensus calls
line_gen_no=0
first_line_gen = TRUE
for(line_gen in gen3) {
  line_gen_no = line_gen_no + 1
  if(opt$verbose) {cat(paste0("Finding consensus call at Gen 3: ",line_gen,"\n"))}
  if (first_line_gen) {
    #all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    no_m_3=ifelse(all_reps_meth_status[,line_gen]=="M",1,0)
    no_u_3=ifelse(all_reps_meth_status[,line_gen]=="U",1,0)
    first_line_gen=FALSE
  } else {
    #all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    no_m_3=no_m_3 + ifelse(all_reps_meth_status[,line_gen]=="M",1,0)
    no_u_3=no_u_3 + ifelse(all_reps_meth_status[,line_gen]=="U",1,0)
  }
}
parental_consensus = ifelse(no_m_3>=2,"M",ifelse(no_u_3>=2,"U","I"))

# Identify sites where at least one offspring line differs from parental consensus, and also, count how many offspring differ at the locus
line_gen_no=0
first_line_gen = TRUE
site_varies_from_parentals = NA
line_varies_from_parentals = list()
num_lines_varying_from_parentals = NA
num_lines_m_to_u_from_parentals = NA
num_lines_u_to_m_from_parentals = NA
num_lines_m_or_u_excl_parentals = NA

for(line_gen in gen3_30) {
  line_gen_no = line_gen_no + 1
  if(opt$verbose) {cat(paste0("Checking variation in methylation calls from parental consensus at Gen 30: ",line_gen,"\n"))}
 
  line_varies_from_parentals[[line_gen]]=ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),TRUE,FALSE)
  if (first_line_gen) {
    #all_m_or_u_3=ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    site_varies_from_parentals=ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),TRUE,FALSE)
    num_lines_varying_from_parentals=ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),1,0)
    num_lines_m_to_u_from_parentals=ifelse((all_reps_meth_status[,line_gen] == "U") & (parental_consensus == "M"),1,0)
    num_lines_u_to_m_from_parentals=ifelse((all_reps_meth_status[,line_gen] == "M") & (parental_consensus == "U"),1,0)
    num_lines_m_or_u_excl_parentals=ifelse((all_reps_meth_status[,line_gen] == "M") | (all_reps_meth_status[,line_gen] == "U"),1,0)
    first_line_gen=FALSE
  } else {
    #all_m_or_u_3=all_m_or_u_3 & ifelse(((all_reps_meth_status[,line_gen]=="U") | (all_reps_meth_status[,line_gen]=="M")),TRUE,FALSE)
    site_varies_from_parentals=site_varies_from_parentals | ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),TRUE,FALSE)
    num_lines_varying_from_parentals=num_lines_varying_from_parentals + ifelse((all_reps_meth_status[,line_gen] %in% c("M","U")) & (parental_consensus %in% c("M","U")) & (all_reps_meth_status[,line_gen]!=parental_consensus),1,0)
    num_lines_m_to_u_from_parentals=num_lines_m_to_u_from_parentals + ifelse((all_reps_meth_status[,line_gen] == "U") & (parental_consensus == "M"),1,0)
    num_lines_u_to_m_from_parentals=num_lines_u_to_m_from_parentals + ifelse((all_reps_meth_status[,line_gen] == "M") & (parental_consensus == "U"),1,0)
    num_lines_m_or_u_excl_parentals=num_lines_m_or_u_excl_parentals + ifelse((all_reps_meth_status[,line_gen] == "M") | (all_reps_meth_status[,line_gen] == "U"),1,0)
  }
}

# Plot histograms of numbers of times sites change in independent lines
qq=num_lines_m_to_u_from_parentals[num_lines_m_to_u_from_parentals>0]
ggplot(data.frame(qq)) + geom_histogram(aes(x=qq), binwidth=1) +xlab("No. lines site changes M->U")
qq=num_lines_u_to_m_from_parentals[num_lines_u_to_m_from_parentals>0]
ggplot(data.frame(qq)) + geom_histogram(aes(x=qq), binwidth=1) +xlab("No. lines site changes U->M")
hist(num_lines_m_to_u_from_parentals)$counts
# [1] 2766341       0       0       0   28868       0       0       0       0    5532       0       0       0       0    1223       0 
#[17]       0       0       0     229
hist(num_lines_u_to_m_from_parentals)$counts
# [1] 2767277       0       0       0   26902       0       0       0       0    5426       0       0       0       0    1967       0
#[17]       0       0       0     621


if (opt$verbose) {
  cat(paste("No of",meth_context,"loci:",nrow(all_reps_meth_status),"\n"))
  cat(paste("No of",meth_context,"loci with all calls M or U:",nrow(all_reps_meth_status[all_m_or_u==TRUE,]),"\n"))
  cat(paste("No of",meth_context,"loci with all calls M:",nrow(all_reps_meth_status[all_m==TRUE,]),"\n"))
  cat(paste("No of",meth_context,"loci with all calls U:",nrow(all_reps_meth_status[all_u==TRUE,]),"\n"))
  cat(paste("No of",meth_context,"loci with strict variation in M/U calls (all parentals agree, all lines have agreement between reps):",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE)),]),"\n"))
  cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3==TRUE)),]),"\n"))
  cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3==TRUE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3_30==TRUE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3_30==TRUE)),]),"\n"))
  cat(paste("No of",meth_context,"loci with any variation in M/U calls (at least 2 parentals agree, at least one progeny line disagrees and has agreement between reps):",nrow(all_reps_meth_status[site_varies_from_parentals,]),"\n"))
  cat(paste("No of",meth_context,"loci with any variation in M/U calls, M at generation 3:",nrow(all_reps_meth_status[site_varies_from_parentals & (no_m_3>=2),]),"\n"))
  cat(paste("No of",meth_context,"loci with any variation in M/U calls, U at generation 3:",nrow(all_reps_meth_status[site_varies_from_parentals & (no_u_3>=2),]),"\n"))
  cat(paste("No of",meth_context,"loci with a clear call in parental consensus and all offspring:",nrow(cbind.data.frame(parental_consensus, site_varies_from_parentals, num_lines_m_to_u_from_parentals, num_lines_u_to_m_from_parentals, num_lines_m_or_u_excl_parentals)[(num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),]),"\n"))
  # Commented out reporting of generation 31 changes for now.  We should redo this analysis in a more sensitive way - check each generation 31 line for changes compared to the corresponding generation 30 line, in concert with reps, but ignoring coverage and calls in generation 3 and in all other generation 30/31 lines.
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, M at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3_30==TRUE)),]),"\n"))
  #cat(paste("No of",meth_context,"loci with strict variation in M/U calls, U at generation 3 and 30:",nrow(all_reps_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3_30==TRUE)),]),"\n"))
}

# Generate tables of sites with various patterns of calls
across_samples_meth_status=cbind(all_samples_meth_status[,1:3],"average_methylation"=average_methylation[,1]/average_methylation[,2])[average_methylation[,2]>0,]
#variant_calls=all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE)),]
variant_calls=all_samples_meth_status[site_varies_from_parentals,]
#M_to_U_calls=all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_m_3==TRUE)),]
#U_to_M_calls=all_samples_meth_status[((all_m_or_u==TRUE) & (all_u==FALSE) & (all_m==FALSE) & (all_u_3==TRUE)),]
M_to_U_calls=all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),]
U_to_M_calls=all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),]
all_M_calls=all_samples_meth_status[all_m==TRUE,]
all_U_calls=all_samples_meth_status[all_u==TRUE,]
all_clean_calls=all_samples_meth_status[(num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),]

M_parent_calls = all_samples_meth_status[parental_consensus=="M",]
U_parent_calls = all_samples_meth_status[parental_consensus=="U",]

# Make loci 2nt long for CG sites
if (meth_context=="CG") {
  across_samples_meth_status=within(across_samples_meth_status, {Locus2=Locus+1})
  variant_calls=within(variant_calls, {Locus2=Locus+1})
  M_to_U_calls=within(M_to_U_calls, {Locus2=Locus+1})
  U_to_M_calls=within(U_to_M_calls, {Locus2=Locus+1})
  all_M_calls=within(all_M_calls, {Locus2=Locus+1})
  all_U_calls=within(all_U_calls, {Locus2=Locus+1})
  all_clean_calls=within(all_clean_calls, {Locus2=Locus+1})
  
  M_parent_calls=within(M_parent_calls, {Locus2=Locus+1})
  U_parent_calls=within(U_parent_calls, {Locus2=Locus+1})
} else {
  across_samples_meth_status=within(across_samples_meth_status, {Locus2=Locus})
  variant_calls=within(variant_calls, {Locus2=Locus})
  M_to_U_calls=within(M_to_U_calls, {Locus2=Locus})
  U_to_M_calls=within(M_to_U_calls, {Locus2=Locus})
  all_M_calls=within(all_M_calls, {Locus2=Locus})
  all_U_calls=within(all_U_calls, {Locus2=Locus})
  all_clean_calls=within(all_clean_calls, {Locus2=Locus})
  
  M_parent_calls=within(M_parent_calls, {Locus2=Locus})
  U_parent_calls=within(U_parent_calls, {Locus2=Locus})
}

# Load the annotation
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
gff.genes$V1=toupper(gff.genes$V1)
gff.transposons$V1=toupper(gff.transposons$V1)
gff.exons$V1=toupper(gff.exons$V1)
gff.introns$V1=toupper(gff.introns$V1)
gff.5UTR$V1=toupper(gff.5UTR$V1)
gff.3UTR$V1=toupper(gff.3UTR$V1)

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
# commented out as sqldf is not installed on cluster R_zg-1.0.0
#library(sqldf)
#gff.genes=merge(gff.genes, sqldf('SELECT genes.gene_id, MAX(exon_no) AS no_exons FROM [gff.exons] exons, [gff.genes] genes WHERE genes.gene_ID=exons.gene_ID GROUP BY genes.gene_id'), by="gene_ID", all=TRUE)


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
m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for U->M sites
u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make GRanges for each individual line's changes from parentals
line_m_to_u_call_ranges = list()
line_u_to_m_call_ranges = list()
for (this_line in gen3_30) {
	line_m_to_u_call_ranges[[this_line]] = makeGRangesFromDataFrame(df = all_samples_meth_status[line_varies_from_parentals[[this_line]] & (no_m_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
	line_u_to_m_call_ranges[[this_line]] = makeGRangesFromDataFrame(df = all_samples_meth_status[line_varies_from_parentals[[this_line]] & (no_u_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
}

# Make a GRanges for clean call sites
clean_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make 'clean' versions of the site ranges, where data is available for all lines
# No need to make for all_M and all_U as they were already 'clean'
#clean_all_M_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(all_m==TRUE) & (num_lines_m_or_u_excl_parentals==4),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
#clean_all_U_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(all_u==TRUE) & (num_lines_m_or_u_excl_parentals==4),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
clean_variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
clean_m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
clean_u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make single change clean versions of the variable site ranges - sites where data is available for all lines, and where a change occurs but only in one offspring line
clean_single_variant_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus %in% c("M","U")),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
clean_single_m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (no_m_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
clean_single_u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[(num_lines_varying_from_parentals == 1) & (no_u_3>=2) & (num_lines_m_or_u_excl_parentals==4) & (parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")

M_parent_call_ranges = makeGRangesFromDataFrame(df = all_samples_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
U_parent_call_ranges = makeGRangesFromDataFrame(df = all_samples_meth_status[(parental_consensus == "U"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")

# DMRs are differentially methylated regions as assessed by the Schmitz paper
# We exclude these from further consideration
# We also add in the sites whose methylation Lizzie has already found to be highly labile
DMRs.gr = reduce(c(makeGRangesFromDataFrame(df=read.table(file=Schmitz_CG_DMRs, header=FALSE, sep="\t"), start.field = "V2", end.field = "V3",seqnames.field = "V1"),makeGRangesFromDataFrame(df=read.table(file=Schmitz_non_CG_DMRs, header=FALSE, sep="\t"), start.field = "V2", end.field = "V3",seqnames.field = "V1")))
levels(DMRs.gr@seqnames@values) = toupper(as.character(levels(DMRs.gr@seqnames@values)))
DMRs.gr@seqinfo@seqnames = levels(DMRs.gr@seqnames@values)
labile_sites = readRDS("../SeqRunGenAll_CG_p10_2cutoff_unique_triple_overlap_sites_withoutmywt_allsitesnotjustvariables.rds")
DMRs.gr = reduce(c(DMRs.gr, makeGRangesFromDataFrame(df=labile_sites, start.field="Locus", end.field="Locus", seqnames.field="Chromosome")))

# Remove DMRs.gr from consideration
variant_call_ranges = setdiff(variant_call_ranges, DMRs.gr, ignore.strand=TRUE)
m_to_u_call_ranges = setdiff(m_to_u_call_ranges, DMRs.gr, ignore.strand=TRUE)
u_to_m_call_ranges = setdiff(u_to_m_call_ranges, DMRs.gr, ignore.strand=TRUE)
for (this_line in gen3_30) {
	line_m_to_u_call_ranges[[this_line]] = setdiff(line_m_to_u_call_ranges[[this_line]], DMRs.gr, ignore.strand=TRUE)
	line_u_to_m_call_ranges[[this_line]] = setdiff(line_u_to_m_call_ranges[[this_line]], DMRs.gr, ignore.strand=TRUE)
}
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
olaps_parent_M = findOverlaps(M_parent_call_ranges, gene_ranges)
olaps_parent_U = findOverlaps(U_parent_call_ranges, gene_ranges)

tolaps_CG_sites = findOverlaps(CG_site_ranges, transposon_ranges)
tolaps_all_M = findOverlaps(all_M_ranges, transposon_ranges)
tolaps_all_U = findOverlaps(all_U_ranges, transposon_ranges)
tolaps_variable = findOverlaps(variant_call_ranges, transposon_ranges)
tolaps_m_to_u = findOverlaps(m_to_u_call_ranges, transposon_ranges)
tolaps_u_to_m = findOverlaps(u_to_m_call_ranges, transposon_ranges)
tolaps_parent_M = findOverlaps(M_parent_call_ranges, transposon_ranges)
tolaps_parent_U = findOverlaps(U_parent_call_ranges, transposon_ranges)

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

# Generate a count of parental_M sites per gene
gff.genes$parental_M_count=table(gff.genes[subjectHits(olaps_parent_M),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$parental_M_count=ifelse(is.na(gff.genes$parental_M_count),0,gff.genes$parental_M_count)

# Generate a count of all_U sites per gene
gff.genes$all_U_count=table(gff.genes[subjectHits(olaps_all_U),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$all_U_count=ifelse(is.na(gff.genes$all_U_count),0,gff.genes$all_U_count)

# Generate a count of parental_U sites per gene
gff.genes$parental_U_count=table(gff.genes[subjectHits(olaps_parent_U),]$gene_ID)[gff.genes$gene_ID]
# Replace NA values in variant_count with 0s
gff.genes$parental_U_count=ifelse(is.na(gff.genes$parental_U_count),0,gff.genes$parental_U_count)

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
gff.genes$average_methylation=table(gff.genes[subjectHits(olaps_CG_sites),]$gene_ID)[gff.genes$gene_ID]

# Find out which transposons are really transposons, and which are 'euchromatic transposons' as per my own method.
gff.transposons$average_methylation=table(gff.transposons[subjectHits(tolaps_CG_sites),]$gene_ID)[gff.transposons$gene_ID]


# Create a data frame with gene ID and average methylation for each overlap between a site and a gene model
olaps_average_methylation=data.frame(cbind("gene_ID"=gff.genes[subjectHits(olaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(olaps_CG_sites)]))
tolaps_average_methylation=data.frame(cbind("gene_ID"=gff.transposons[subjectHits(tolaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(tolaps_CG_sites)]))
# Convert average values to numerics
#olaps_average_methylation$site_average_methylation=as.numeric(levels(olaps_average_methylation$site_average_methylation)[olaps_average_methylation$site_average_methylation])
#tolaps_average_methylation$site_average_methylation=as.numeric(levels(tolaps_average_methylation$site_average_methylation)[tolaps_average_methylation$site_average_methylation])
olaps_average_methylation$site_average_methylation=as.numeric(olaps_average_methylation$site_average_methylation)
tolaps_average_methylation$site_average_methylation=as.numeric(tolaps_average_methylation$site_average_methylation)

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

cat(paste0("Fit mu values (means):\n"))
mixmdl$mu
# Schmitz data: 0.01777965 0.10733632 0.64071733.  After converting P to U/M with mCG=0.45 cutoff: 0.08997885 0.65362097. Using Fishers p<0.05, Binomial p<0.005: 0.01632469 0.10906799 0.64090168
# Becker data:  0.003837026 0.070503558 0.635220579

cat(paste0("Fit sigma values (st.devs):\n"))
mixmdl$sigma
# Schmitz data: 0.01006439 0.05852205 0.14541845.  After converting P to U/M with mCG=0.45 cutoff: 0.07725693 0.14292752. Using Fishers p<0.05, Binomial p<0.005: 0.01027012 0.05664957 0.14650136
# Becker data:  0.004954326 0.045214098 0.165311575

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


# Plot gene methylation level vs. gene length, colour by chromosome
heterochromatic_gene_coverage_cutoff=0.6
heterochromatic_gene_coverage_minimum=0.2

gff.genes$m_class = ifelse(gff.genes$CG_sites_count<=1,"CG-poor",ifelse(gff.genes$all_M_count+gff.genes$all_U_count+gff.genes$variable_count==0,"Unknown",ifelse(gff.genes$average_methylation==0,"Unmethylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation>heterochromatic_gene_dispersion_cutoff),"Gene-body Methylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.genes$average_methylation<heterochromatic_gene_coverage_minimum),"Unmethylated",ifelse((gff.genes$variance_methylation/gff.genes$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.genes$average_methylation>=heterochromatic_gene_coverage_minimum),"Heterochromatic","Really Unknown"))))))

gff.transposons$m_class = ifelse(gff.transposons$CG_sites_count<=1,"CG-poor",ifelse(gff.transposons$all_M_count+gff.transposons$all_U_count+gff.transposons$variable_count==0,"Unknown",ifelse(gff.transposons$average_methylation==0,"Unmethylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation>heterochromatic_gene_dispersion_cutoff),"Gene-body Methylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.transposons$average_methylation<heterochromatic_gene_coverage_minimum),"Unmethylated",ifelse((gff.transposons$variance_methylation/gff.transposons$average_methylation<=heterochromatic_gene_dispersion_cutoff) & (gff.transposons$average_methylation>=heterochromatic_gene_coverage_minimum),"Heterochromatic","Really Unknown"))))))








#gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))
gene_body_loci.gr=makeGRangesFromDataFrame(df=na.omit(gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)




segmentation_model.gr = readRDS(file = "SRA035939_CG_segmentation_model_draft3.rds")

# Capitalise chromosome names in segmentation_model.gr
levels(segmentation_model.gr@seqnames)=toupper(levels(segmentation_model.gr@seqnames))
segmentation_model.gr@seqinfo@seqnames=levels(segmentation_model.gr@seqnames)

#levels(variant_call_ranges@seqnames)=toupper(levels(variant_call_ranges@seqnames))
#variant_call_ranges@seqinfo@seqnames=levels(variant_call_ranges@seqnames)


GBM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]
TEM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="TEM"]
UMR_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="UMR"]
length(GBM_segments.gr)
length(TEM_segments.gr)
length(UMR_segments.gr)

library(data.table) 

# count how many M sites each gene has in parental lines:
# removed this line because it is now done earlier as gff.genes$parental_M_count
#gff.genes$M_parent_calls = table(gff.genes[subjectHits(findOverlaps(M_parent_call_ranges, gene_ranges)),]$gene_ID)[gff.genes$gene_ID]

#  Identify de novo methylated sites in previously unmethylated regions

#define size of zones at each end of the unmethylated segment, to exclude cooperative effects from neighbouring methylated segments

#define size of extended 5' zone to accommodate region which is never methylated and has high density of unmethylated CG sites
extended_L_5_prime_length_bp = 2000
extended_L_3_prime_length_bp = 1000
#extended_L_5_prime_length = 40
#extended_L_3_prime_length = 20
extended_L_5_prime_length = 20
extended_L_3_prime_length = 10

# Make a set of genomic ranges to represent the masked ends of genes
end_ranges_5_prime = makeGRangesFromDataFrame(df = as.data.frame(cbind("chromosome" = gff.genes$V1, "start" = ifelse(gff.genes$V7 == "+", gff.genes$V4, ifelse((gff.genes$V5-extended_L_5_prime_length_bp)<gff.genes$V4, gff.genes$V4, gff.genes$V5-extended_L_5_prime_length_bp)), "end" = ifelse(gff.genes$V7 == "+", ifelse((gff.genes$V4+extended_L_5_prime_length_bp)>gff.genes$V5,gff.genes$V5,gff.genes$V4+extended_L_5_prime_length_bp), gff.genes$V5))))#, start_field = "start", end_field = "end")
end_ranges_3_prime = makeGRangesFromDataFrame(df = as.data.frame(cbind("chromosome" = gff.genes$V1, "start" = ifelse(gff.genes$V7 == "+", ifelse((gff.genes$V5-extended_L_3_prime_length_bp)<gff.genes$V4,gff.genes$V4,gff.genes$V5-extended_L_3_prime_length_bp), gff.genes$V4), "end" = ifelse(gff.genes$V7 == "+", gff.genes$V5, ifelse((gff.genes$V4+extended_L_3_prime_length_bp)>gff.genes$V5,gff.genes$V5,gff.genes$V4+extended_L_3_prime_length_bp)))))#, start_field = "start", end_field = "end")
	
masked_ranges = setdiff(gene_ranges, union(end_ranges_5_prime, end_ranges_3_prime, ignore.strand = TRUE), ignore.strand = TRUE)
# 7064 genes have some part of gene body remaining after masking ends
	
pdf(paste0(project_id,"_length_distribution_of_masked_genes.pdf"))
ggplot(as.data.frame(masked_ranges)) + geom_histogram(aes(x=width), binwidth=100) + xlab("width of masked gene body (2kb 5', 1kb 3')")
dev.off()

# Make segments to represent unmethylated genes
unmethylated_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=na.omit(gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))

meth_quartiles = c(0,quantile(gff.genes[gff.genes$m_class=="Gene-body Methylated","average_methylation"],c(0.25, 0.5, 0.75, type=1), na.rm=TRUE))
meth_quartiles
#                  25%        50%        75%       100% 
#Schmitz:  0.00000000 0.08759728 0.18859130 0.32126244 0.77777778 
#Becker:   0.00000000 0.07471366 0.17397265 0.30336126 0.66666667 


# Identify genes containing only a single segment of a given length

# genes which overlap one or more GBM singleton segments
#961 genes overlap 1 or more singleton GBM segments
#895 of these overlap no TEM segments
#324 of these only overlap a single GBM segment
#Of which 84 are valid with 40,20 buffer, 195 valid with 20,10

# genes which overlap a single 2 site segment
#564 genes

# genes which overlap a single 3 site segment
#508 genes

# genes which overlap a single 4 site segment
#353 genes

# genes which overlap a single 5 site segment
#273 genes

# genes which overlap a single 6 site segment
#230 genes





### Redo the analysis using genes with a single gbM segment of specified number of parental methylated sites
### Plot the genes with 1, 2, 3 parental M sites in phylogenetic context

# Define a convenience function to convert chromosome names to mixed-case (e.g. Chr9, ChrC)
mixedCaseChr <- function(s, strict = FALSE) {
	paste0(toupper(substr(s,1,1)),tolower(substr(s,2,3)),toupper(substr(s,4,4)))
}

# Read in the matrix of meth status values
# genes only:
#tree_samples_meth_status <- readRDS(paste0("4-summary/AncestralAnalysis_CG_gene_CG_genmeth.rds"))

# all sites in genome:
tree_samples_meth_status <- readRDS(paste0("../AncestralAnalysis_CG_all_samples_meth_status.rds"))

# all sites in genome in runs including root taxa:
tree_samples_meth_status2 <- readRDS(paste0("../SeqRunGenAll_CG_all_samples_meth_status.rds"))

# Merge the three sets of 'root taxa' data with the other samples
tree_samples_meth_status = cbind(tree_samples_meth_status, tree_samples_meth_status2[,4:6])
rm(tree_samples_meth_status2)



# previously calculated genome-wide set:
#tree_samples_meth_status = valid_samples_meth_status

tree_samples_meth_status$Chromosome = mixedCaseChr(tree_samples_meth_status$Chromosome)
tree_samples_meth_status$Chromosome = toupper(tree_samples_meth_status$Chromosome)

# Read in the segmentation model in case we need it later
#segmentation_model.gr = readRDS(paste0("4-summary/SRA035939_CG_segmentation_model_draft3.rds"))
levels(segmentation_model.gr@seqnames@values) = mixedCaseChr(as.character(levels(segmentation_model.gr@seqnames@values)))
levels(segmentation_model.gr@seqnames@values) = toupper(as.character(levels(segmentation_model.gr@seqnames@values)))
segmentation_model.gr@seqinfo@seqnames = levels(segmentation_model.gr@seqnames@values)

# Generate a masked segmentation model
### this part might fail due to capitalisation of chromosome ID
#segment_mask_loci.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="Chr2", "start_site"=3239693, "end_site"=3505260))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")		
segment_mask_loci.gr = makeGRangesFromDataFrame(df=as.data.frame(rbind(c("Chromosome"="CHR2", "start_site"=3239693, "end_site"=3505260))), start.field="start_site", end.field="end_site", seqnames.field="Chromosome")		
masked_segmentation_model.gr = segmentation_model.gr[unique(queryHits(findOverlaps(segmentation_model.gr, setdiff(segmentation_model.gr,segment_mask_loci.gr))))]



# Split GBM segments by segment length (no CG sites)
# Calculate probability of loss / site in cluster
#no_losses_observed/(no_segments*length_of_segment) by segment length group
previous_segment_length = 0 
GBM_segments.gr = masked_segmentation_model.gr[(masked_segmentation_model.gr$segment.mStatus=="GBM") & (!is.na(masked_segmentation_model.gr$nSites))]
segment_length_groups = c(1,2,3,4,5,6,7,8,9,10,11, 12, 13, 14, 15,20,25,30,35,40,45,50,100,500)
for (this_segment_length in segment_length_groups) {
  this_segment_group = GBM_segments.gr[(GBM_segments.gr$nSites>previous_segment_length) & (GBM_segments.gr$nSites<=this_segment_length)]
  cat(paste(this_segment_length, length(this_segment_group), sum(this_segment_group$nSites), length(findOverlaps(clean_m_to_u_call_ranges, this_segment_group)), length(findOverlaps(clean_m_to_u_call_ranges, this_segment_group))/sum(this_segment_group$nSites)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)), length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[,], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group)))
  
  # above figures are related to length of segment.  instead let's consider no. M sites in segment:
  #length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), this_segment_group))
  
  for (this_line_no in 1:length(gen3_30)) {
    this_line = gen3_30[this_line_no]
    
    # find the subset of m_to_u and u_to_m which are in the target zone - should speed up the overlap searches later
    #line_u_to_m_target_ranges = line_u_to_m_call_ranges[[this_line]][queryHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
    #line_m_to_u_target_ranges = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
    
    line_losses = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], this_segment_group))]
    line_prop_losses = length(findOverlaps(makeGRangesFromDataFrame(df=all_reps_meth_status[(parental_consensus == "M"),], start.field = "Locus", end.field = "Locus", seqnames.field="Chromosome"), line_losses))
    cat(paste0(" ",length(line_losses), " ", line_prop_losses))
    # line_losses and line_prop_losses are always coming out the same, but included just in case there is a problem
  }
  cat("\n")
  previous_segment_length=this_segment_length
}




# Read in the table of genes with their characteristics
#gff.genes =  readRDS(paste0("4-summary/SRA035939_CG_gff.genes.rds"))

# This object never gets used so should probably comment this section out:
# Make a set of genomic ranges corresponding to GBM genes
GBM_gene_loci=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")]
colnames(GBM_gene_loci)=c("Chromosome", "start_site", "end_site")
GBM_gene_loci$Chromosome=mixedCaseChr(GBM_gene_loci$Chromosome)
GBM_gene_loci.gr=reduce(makeGRangesFromDataFrame(df=GBM_gene_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome"))

# Make a set of genomic ranges corresponding to the sites in the matrix of meth status values
meth_status_sites.gr = makeGRangesFromDataFrame(df=tree_samples_meth_status, start.field="Locus", end.field="Locus", seqnames.field="Chromosome")



library(phangorn)

# Read RAxML rooted tree
#labelled_tree=read.tree(file="5-analysis/RAxML_RootedTree.draft_tree_rooted")
#labelled_tree=read.tree(file="../RAxML_RootedTree.draft_tree_rooted")
labelled_tree=read.tree(file="rerooted_tree.nwk")
plot(labelled_tree, cex=0.5)
nodelabels(labelled_tree$node.label, cex=0.5, frame="none")

# labelled tree has JAY_1 listed as JAY_1_0 for some reason due to the rerooting I think.  Need to fix this so names match up:
labelled_tree$tip.label[which(labelled_tree$tip.label=="JAY_1_0")] = "JAY_1"

# Alternative plotting using ggtree package

library(ggtree)

#MRCA(labelled_tree, tip=c("SRR342385", "SRR342381"))
#105
# Plot the basic tree
tree_plot = ggtree(labelled_tree,) 
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
tree_plot = tree_plot + geom_cladelabel(node=163, label="MA lines", color="green") + geom_hilight(node=163, fill="green")
# Add colour and label for SALK lines clade
tree_plot = tree_plot + geom_cladelabel(node=112, label="SALK lines", color="pink") + geom_hilight(node=112, fill="pink")
tree_plot = tree_plot + geom_cladelabel(node=195, label="SALK lines", color="pink") + geom_hilight(node=195, fill="pink")
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


# iterate over buffer_length and meth_quartile
no_stat_entries = 0				
segments_finished = 0

for (buffer_length in c(10)) {

	for (gbm_segment_length in 0:10) {
		cat("buffer length:",buffer_length,"\n")
		# meth_quartile 0 means do the unmethylated genes; 1-k means do the relevant subset of gene-body methylated genes by length of GBM segment

		
		pdf(paste0("m_pattern_gbM_segments_length_",gbm_segment_length,".pdf"))
		plot.new()


		
		# identify segments of interest
		#target_segments = unmethylated_segments∩gene-body_methylated_genes
		if (gbm_segment_length==0) {
			# Do unmethylated genes
			# Make segments to represent relevant target genes
			#unmethylated_gene_loci.gr = reduce(makeGRangesFromDataFrame(df=na.omit(gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))
			target_segments.gr = reduce(makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Unmethylated") ,c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))[1:10]
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
			target_segments.gr = reduce(makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$gene_ID %in% names(table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[!is.na(segmentation_model.gr$nSites)&(segmentation_model.gr$nSites==gbm_segment_length)&(segmentation_model.gr$segment.mStatus=="GBM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))]))) & !(gff.genes$gene_ID %in% names(table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[(segmentation_model.gr$segment.mStatus=="TEM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))]))) & (gff.genes$gene_ID %in% names(table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[(segmentation_model.gr$segment.mStatus=="GBM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))]))[table(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),"gene_ID"][subjectHits(findOverlaps(segmentation_model.gr[(segmentation_model.gr$segment.mStatus=="GBM"),],makeGRangesFromDataFrame(df=na.omit(gff.genes[(gff.genes$m_class=="Gene-body Methylated"),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)))])==1]),c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))
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

#segment_change_stats$core_U_M = as.numeric(levels(segment_change_stats$core_U_M)[segment_change_stats$core_U_M])
#segment_change_stats$core_M_U = as.numeric(levels(segment_change_stats$core_M_U)[segment_change_stats$core_M_U])
#segment_change_stats$core_length = as.numeric(levels(segment_change_stats$core_length)[segment_change_stats$core_length])
#segment_change_stats$core_length_bp = as.numeric(levels(segment_change_stats$core_length_bp)[segment_change_stats$core_length_bp])
#segment_change_stats$buffer_length = as.numeric(levels(segment_change_stats$buffer_length)[segment_change_stats$buffer_length])
segment_change_stats$core_U_M = as.numeric(segment_change_stats$core_U_M)
segment_change_stats$core_M_U = as.numeric(segment_change_stats$core_M_U)
segment_change_stats$core_length = as.numeric(segment_change_stats$core_length)
segment_change_stats$core_length_bp = as.numeric(segment_change_stats$core_length_bp)
segment_change_stats$buffer_length = as.numeric(segment_change_stats$buffer_length)

saveRDS(segment_change_stats, file=paste0(project_id,"_",meth_context,"_segment_change_stats_by_gbM_segment_length.rds"))

pdf(paste0(project_id,"_distribution_of_M-U_core_transitions_by_gbM_segment_length_P.pdf"))
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
			#changed_site_counts = hist(segment_change_stats[(segment_change_stats$buffer_length==this_buffer_length) & (segment_change_stats$gbm_segment_length==this_gbm_segment_length) & (segment_change_stats$line == this_line) & (!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),"core_U_M"], breaks=seq(from=-0.5, to=10.5, by=1), plot=plot_histograms)$counts
		
			# Calculate proportion of sites losing methylation
			changed_site_counts = hist(segment_change_stats[(segment_change_stats$buffer_length==this_buffer_length) & (segment_change_stats$gbm_segment_length==this_gbm_segment_length) & (segment_change_stats$line == this_line) & (!is.na(segment_change_stats$segment_CG_sites)) & (segment_change_stats$L_5_prime_U_M==0) & (segment_change_stats$L_3_prime_U_M==0),"core_M_U"], breaks=seq(from=-0.5, to=10.5, by=1), plot=plot_histograms)$counts
		
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



target_ranges = setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])

### should we exclude gene_body_loci.gr if they only have 1 or 2 GBM sites? no because then we'd lose cooperative maintenance effects

# Identify GenomicRegions covering 'body' of such GBM genes (excluding 2kb 5' and 1kb 3')
#target_ranges = 

# Identify GenomicRegions covering GBM range of such GBM genes (from first to last methylated base in gen 3)
#target_ranges = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])))]

# Identify gbM segments
#target_ranges = GBM_segments.gr


no_steps = 20
first_step = 1
step_size = 50
no_lines = length(gen3_30)
loss_rate_stats = matrix(0, nrow = no_steps+2, ncol = no_lines)
gain_rate_stats = matrix(0, nrow = no_steps+2, ncol = no_lines)
all_rate_stats = data.frame()
for (this_shift in 1:no_steps) { 
  #cat(paste0(this_shift*step_size," "))
  target_ranges = c(shift(flank(GBM_segments.gr, step_size), -step_size*this_shift), shift(flank(GBM_segments.gr, step_size, start=FALSE), step_size*this_shift))
  
  # ensure none of this overlaps with neighbouring segment flanks
  target_ranges = setdiff(target_ranges, c(GBM_segments.gr, TEM_segments.gr, flank(GBM_segments.gr, step_size*(this_shift - 1)), flank(GBM_segments.gr, step_size*(1-this_shift), start=FALSE)))

  #for (this_line_no in 1:length(gen3_30)) {
  for (this_line_no in 1:no_lines) {
    this_line = gen3_30[this_line_no]
    #cat(paste0(this_line," "))
    
    #cat(length(target_ranges))
    target_stats = merge(table(subjectHits(findOverlaps(U_parent_call_ranges, target_ranges))), table(subjectHits(findOverlaps(M_parent_call_ranges, target_ranges))), by="Var1", all=TRUE)
    colnames(target_stats) = c("Var1", "U_count", "M_count")
  
    # find the subset of m_to_u and u_to_m which are in the target zone - should speed up the overlap searches later
    #line_u_to_m_target_ranges = line_u_to_m_call_ranges[[this_line]][queryHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
    #line_m_to_u_target_ranges = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
  
    if(length(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_ranges)) >0) {
      target_stats = merge(target_stats, table(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_ranges))), by="Var1", all=TRUE)
    } else {target_stats = cbind(target_stats, rep(NA, nrow(target_stats)))}
    colnames(target_stats)[4] = "m_to_u_count"
    
    if(length(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_ranges)) >0) {
      target_stats = merge(target_stats,  table(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_ranges))), by="Var1", all=TRUE)
    } else {target_stats = cbind(target_stats, rep(NA, nrow(target_stats)))}
    colnames(target_stats)[5] = "u_to_m_count"
    
    target_stats[is.na(target_stats$M_count),"M_count"] = 0
    target_stats[is.na(target_stats$U_count),"U_count"] = 0
    target_stats[is.na(target_stats$m_to_u_count),"m_to_u_count"] = 0
    target_stats[is.na(target_stats$u_to_m_count),"u_to_m_count"] = 0
 
    #library(ggplot2)
    #pdf(paste0(project_id,"_",this_line,"_density_of_loss_and_gain_rate_across_gbM_segments.pdf"))

    #print(ggplot(target_stats) + geom_histogram(aes(x=m_to_u_count/M_count), bins=100))

    #print(ggplot(target_stats) + geom_histogram(aes(x=u_to_m_count/U_count), bins=100))
  
    #dev.off() 

    loss_rates = target_stats$m_to_u_count/target_stats$M_count
    gain_rates = target_stats$u_to_m_count/target_stats$U_count

    # These are mean rates across any/all lines
    #cat(paste0("Mean loss rate: ",mean(loss_rates, na.rm=TRUE), "\n"))
    #[1] 0.2318511

    #cat(paste0("Mean gain rate: ",mean(gain_rates, na.rm=TRUE), "\n"))
    #[1] 0.1749956

    no_gen=30
    cc_per_gen=30
    #no_lines = 4  #
    #cat(paste0("Mean loss rate: ",mean(loss_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)))
    #[1] 6.440308e-05

    #cat(paste0("Mean gain rate: ",mean(gain_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)))
    #[1] 4.86099e-05

    # Overall average loss and gain rates
    #sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))
    #[1] 0.1487853

    #sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))
    #[1] 0.1412831

    
    cat(paste0(this_shift*step_size," "))
    cat(paste0(this_line_no," "))
    # translated to per cell cycle rates:
    #cat(paste0("Loss rate: ",sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen*no_lines)))
    #[1] 4.132925e-05
    # we don't divide by no_lines any more, as we are doign line by line
    cat(paste0("Loss rate: ",sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)))
    
    #cat(paste0("Gain rate: ",sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen*no_lines)))
    #[1] 3.924532e-05
    cat(paste0("Gain rate: ",sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)))
    cat(paste0("\n"))
    loss_rate_stats[this_shift+1-first_step, this_line_no] = sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)
    gain_rate_stats[this_shift+1-first_step, this_line_no] = sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)
    if (nrow(all_rate_stats) == 0) {
      all_rate_stats = c("context"="GBM segments", "shift"=this_shift*step_size-step_size/2, "line_no"=this_line_no, "type"="lossrate", "rate"=sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen))
    } else {
      all_rate_stats = rbind(all_rate_stats, c("context"="GBM segments", "shift"=this_shift*step_size-step_size/2, "line_no"=this_line_no, "type"="lossrate", "rate"=sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)))
    }
    all_rate_stats = rbind(all_rate_stats, c("context"="GBM segments", "shift"=this_shift*step_size-step_size/2, "line_no"=this_line_no, "type"="gainrate", "rate"=sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)))
  }
}

# add in rows to represent the region inside the segments

#Rates for segments only:
target_ranges = GBM_segments.gr

#Rates for whole genes:
#target_ranges = GBM_gene_loci.gr
#levels(target_ranges@seqnames) = toupper(as.character(levels(target_ranges@seqnames)))
#target_ranges@seqinfo@seqnames = levels(target_ranges@seqnames)


#for (this_line_no in 1:length(gen3_30)) {
for (this_line_no in 1:no_lines) {
  this_line = gen3_30[this_line_no]
  #cat(paste0(this_line," "))
  target_stats = merge(table(subjectHits(findOverlaps(U_parent_call_ranges, target_ranges))), table(subjectHits(findOverlaps(M_parent_call_ranges, target_ranges))), by="Var1", all=TRUE)
  colnames(target_stats) = c("Var1", "U_count", "M_count")
  
  # find the subset of m_to_u and u_to_m which are in the target zone - should speed up the overlap searches later
  #line_u_to_m_target_ranges = line_u_to_m_call_ranges[[this_line]][queryHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
  #line_m_to_u_target_ranges = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
  
  target_stats = merge(target_stats, table(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_ranges))), by="Var1", all=TRUE)
  colnames(target_stats)[4] = "m_to_u_count"
  
  target_stats = merge(target_stats,  table(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_ranges))), by="Var1", all=TRUE)
  colnames(target_stats)[5] = "u_to_m_count"
  
  target_stats[is.na(target_stats$M_count),"M_count"] = 0
  target_stats[is.na(target_stats$U_count),"U_count"] = 0
  target_stats[is.na(target_stats$m_to_u_count),"m_to_u_count"] = 0
  target_stats[is.na(target_stats$u_to_m_count),"u_to_m_count"] = 0
  
  #library(ggplot2)
  #pdf(paste0(project_id,"_",this_line,"_density_of_loss_and_gain_rate_across_gbM_segments.pdf"))
  
  #print(ggplot(target_stats) + geom_histogram(aes(x=m_to_u_count/M_count), bins=100))
  
  #print(ggplot(target_stats) + geom_histogram(aes(x=u_to_m_count/U_count), bins=100))
  
  #dev.off() 
  
  loss_rates = target_stats$m_to_u_count/target_stats$M_count
  gain_rates = target_stats$u_to_m_count/target_stats$U_count
  
  cat(paste0("line: ", this_line," m_to_u_count:", sum(target_stats$m_to_u_count), " m_count:", sum(target_stats$M_count), " u_to_m_count:", sum(target_stats$u_to_m_count), " u_count:", sum(target_stats$U_count), "\n"))
  # These are mean rates across any/all lines
  #cat(paste0("Mean loss rate: ",mean(loss_rates, na.rm=TRUE), "\n"))
  #[1] 0.2318511
  
  #cat(paste0("Mean gain rate: ",mean(gain_rates, na.rm=TRUE), "\n"))
  #[1] 0.1749956
  
  no_gen=30
  cc_per_gen=30
  #no_lines = 4  #
  #cat(paste0("Mean loss rate: ",mean(loss_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)))
  #[1] 6.440308e-05
  
  #cat(paste0("Mean gain rate: ",mean(gain_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)))
  #[1] 4.86099e-05
  
  # Overall average loss and gain rates
  #sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))
  #[1] 0.1487853
  
  #sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))
  #[1] 0.1412831
  
  
  cat(paste0(this_shift*step_size," "))
  cat(paste0(this_line_no," "))
  # translated to per cell cycle rates:
  #cat(paste0("Loss rate: ",sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen*no_lines)))
  #[1] 4.132925e-05
  # we don't divide by no_lines any more, as we are doign line by line
  cat(paste0("Loss rate: ",sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)))
  
  #cat(paste0("Gain rate: ",sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen*no_lines)))
  #[1] 3.924532e-05
  cat(paste0("Gain rate: ",sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)))
  cat(paste0("\n"))
  loss_rate_stats[this_shift+1-first_step, this_line_no] = sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)
  gain_rate_stats[this_shift+1-first_step, this_line_no] = sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)
  if (nrow(all_rate_stats) == 0) {
    all_rate_stats = c("context"="GBM segments", "shift"=0, "line_no"=this_line_no, "type"="lossrate", "rate"=sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen))
  } else {
    all_rate_stats = rbind(all_rate_stats, c("context"="GBM segments", "shift"=0, "line_no"=this_line_no, "type"="lossrate", "rate"=sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)))
  }
  all_rate_stats = rbind(all_rate_stats, c("context"="GBM segments", "shift"=0, "line_no"=this_line_no, "type"="gainrate", "rate"=sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)))
}

#Rates for segments only:
#line: Line29Gen30 m_to_u_count:7550 m_count:186688 u_to_m_count:4354 u_count:107854
#line: Line49Gen30 m_to_u_count:7301 m_count:186688 u_to_m_count:4166 u_count:107854
#line: Line59Gen30 m_to_u_count:7627 m_count:186688 u_to_m_count:4130 u_count:107854
#line: Line119Gen30 m_to_u_count:5129 m_count:186688 u_to_m_count:3425 u_count:107854

#Rates for whole genes:
#line: Line29Gen30 m_to_u_count:7613 m_count:195642 u_to_m_count:6384_ u_count:704832
#line: Line49Gen30 m_to_u_count:7363 m_count:195642 u_to_m_count:6161 u_count:704832
#line: Line59Gen30 m_to_u_count:7659 m_count:195642 u_to_m_count:6037 u_count:704832
#line: Line119Gen30 m_to_u_count:5149 m_count:195642 u_to_m_count:4983 u_count:704832

# add in rows to represent various kinds of UM segments

target_ranges = UMR_segments.gr

# find overlaps between UMR segments and UMR genes
UMR_gene_loci=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")]
colnames(UMR_gene_loci)=c("Chromosome", "start_site", "end_site")
#UMR_gene_loci$Chromosome=mixedCaseChr(UMR_gene_loci$Chromosome)
UMR_gene_loci.gr=reduce(makeGRangesFromDataFrame(df=UMR_gene_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome"))

# find UMR segments which overlap UM genes
target_ranges = UMR_segments.gr[subjectHits(findOverlaps(UMR_gene_loci.gr, UMR_segments.gr))]

#target_ranges = UMR_gene_loci.gr


# Find genes often methylated in 1001 methylomes

# Find genes rarely methylated in 1001 methylomes

# find overlaps between UMR segments and rarely Unmethylated genes

# find overlaps between UMR segments and often Unmethylated genes

this_shift = this_shift + 1

#for (this_line_no in 1:length(gen3_30)) {
for (this_line_no in 1:no_lines) {
  this_line = gen3_30[this_line_no]
  #cat(paste0(this_line," "))
  target_stats = merge(table(subjectHits(findOverlaps(U_parent_call_ranges, target_ranges))), table(subjectHits(findOverlaps(M_parent_call_ranges, target_ranges))), by="Var1", all=TRUE)
  colnames(target_stats) = c("Var1", "U_count", "M_count")
  
  # find the subset of m_to_u and u_to_m which are in the target zone - should speed up the overlap searches later
  #line_u_to_m_target_ranges = line_u_to_m_call_ranges[[this_line]][queryHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
  #line_m_to_u_target_ranges = line_m_to_u_call_ranges[[this_line]][queryHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_segment, ignore.strand=TRUE))]
  
  target_stats = merge(target_stats, table(subjectHits(findOverlaps(line_m_to_u_call_ranges[[this_line]], target_ranges))), by="Var1", all=TRUE)
  colnames(target_stats)[4] = "m_to_u_count"
  
  target_stats = merge(target_stats,  table(subjectHits(findOverlaps(line_u_to_m_call_ranges[[this_line]], target_ranges))), by="Var1", all=TRUE)
  colnames(target_stats)[5] = "u_to_m_count"
  
  target_stats[is.na(target_stats$M_count),"M_count"] = 0
  target_stats[is.na(target_stats$U_count),"U_count"] = 0
  target_stats[is.na(target_stats$m_to_u_count),"m_to_u_count"] = 0
  target_stats[is.na(target_stats$u_to_m_count),"u_to_m_count"] = 0
  
  #library(ggplot2)
  #pdf(paste0(project_id,"_",this_line,"_density_of_loss_and_gain_rate_across_gbM_segments.pdf"))
  
  #print(ggplot(target_stats) + geom_histogram(aes(x=m_to_u_count/M_count), bins=100))
  
  #print(ggplot(target_stats) + geom_histogram(aes(x=u_to_m_count/U_count), bins=100))
  
  #dev.off() 
  
  loss_rates = target_stats$m_to_u_count/target_stats$M_count
  gain_rates = target_stats$u_to_m_count/target_stats$U_count
  
  # These are mean rates across any/all lines
  #cat(paste0("Mean loss rate: ",mean(loss_rates, na.rm=TRUE), "\n"))
  #[1] 0.2318511
  
  #cat(paste0("Mean gain rate: ",mean(gain_rates, na.rm=TRUE), "\n"))
  #[1] 0.1749956
  
  no_gen=30
  cc_per_gen=30
  #no_lines = 4  #
  #cat(paste0("Mean loss rate: ",mean(loss_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)))
  #[1] 6.440308e-05
  
  #cat(paste0("Mean gain rate: ",mean(gain_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)))
  #[1] 4.86099e-05
  
  # Overall average loss and gain rates
  #sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))
  #[1] 0.1487853
  
  #sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))
  #[1] 0.1412831
  
  
  cat(paste0(this_shift*step_size," "))
  cat(paste0(this_line_no," "))
  # translated to per cell cycle rates:
  #cat(paste0("Loss rate: ",sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen*no_lines)))
  #[1] 4.132925e-05
  # we don't divide by no_lines any more, as we are doign line by line
  cat(paste0("Loss rate: ",sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)))
  
  #cat(paste0("Gain rate: ",sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen*no_lines)))
  #[1] 3.924532e-05
  cat(paste0("Gain rate: ",sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)))
  cat(paste0("\n"))
  loss_rate_stats[this_shift+1-first_step, this_line_no] = sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)
  gain_rate_stats[this_shift+1-first_step, this_line_no] = sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)
  if (nrow(all_rate_stats) == 0) {
    all_rate_stats = c("context"="UMR genes", "shift"=this_shift*step_size-step_size/2, "line_no"=this_line_no, "type"="lossrate", "rate"=sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen))
  } else {
    all_rate_stats = rbind(all_rate_stats, c("context"="UMR genes", "shift"=this_shift*step_size-step_size/2, "line_no"=this_line_no, "type"="lossrate", "rate"=sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen)))
  }
  all_rate_stats = rbind(all_rate_stats, c("context"="UMR genes", "shift"=this_shift*step_size-step_size/2, "line_no"=this_line_no, "type"="gainrate", "rate"=sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen)))
}

all_rate_stats = as.data.frame(all_rate_stats)
all_rate_stats$shift = as.numeric(all_rate_stats$shift)
all_rate_stats$rate = as.numeric(all_rate_stats$rate)


# Plot densities of loss and gain rates across all target ranges
# We are plotting rates of change inside gbM segments at x=0, rates in 50bp bins outside segments in x>0, x<1000 and rates in UM segments of UM genes at x>1000
pdf("rates_of_loss_and_gain.pdf")
print(ggplot(as.data.frame(all_rate_stats[all_rate_stats$rate>0,])) + geom_line(aes(x=shift, y=rate, colour=line_no)) + facet_wrap(~type))
# replot gain rates zoomed in 
print(ggplot(as.data.frame(all_rate_stats[(all_rate_stats$rate>0) & (all_rate_stats$type=="gainrate"),])) + geom_line(aes(x=shift, y=rate, colour=line_no)) + facet_wrap(~type) + ylim(0,1e-5))
print(ggplot(as.data.frame(all_rate_stats[(all_rate_stats$rate>0) & (all_rate_stats$type=="lossrate"),])) + geom_line(aes(x=shift, y=rate, colour=line_no)) + facet_wrap(~type) + xlim(0,1000))
dev.off()

write.table(all_rate_stats, file="loss_gain_rates_gbM_boundaries.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)



# We want to plot rates by distance to nearest M on same scale as rates by distance to nearest gbM segment. for comparison
zM=read.table(paste0(project_id,"_dist_from_M_to_nearest_X.tsv"),sep="\t",header=TRUE)
zU=read.table(paste0(project_id,"_dist_from_U_to_nearest_X.tsv"),sep="\t",header=TRUE)

gene_stats = read.table(file="gbm_tem_gene_calls.txt", sep="\t", header=TRUE)


# gain rates
changed_site_counts = hist(zU[zU$u_to_m_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
changed_site_counts = hist(zU[zU$u_to_m_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=12000, by=10))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=12000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 

# loss rates
changed_site_counts = hist(zM[zM$m_to_u_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
changed_site_counts = hist(zM[zM$m_to_u_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=12000, by=10))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=12000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 



# Analysis of rate variation by distance to nearest U site

# gain rates
changed_site_counts = hist(zU[zU$u_to_m_count==1,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts
changed_site_counts = hist(zU[zU$u_to_m_count==1,"dist_to_nearest_U"], breaks=seq(from=0, to 12000, by=10))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"dist_to_nearest_U"], breaks=seq(from=0, to=12000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 

# loss rates
changed_site_counts = hist(zM[zM$m_to_u_count==1,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts
changed_site_counts = hist(zM[zM$m_to_u_count==1,"dist_to_nearest_U"], breaks=seq(from=0, to=12000, by=10))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"dist_to_nearest_U"], breaks=seq(from=0, to=12000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 


#Split zM and zU up by line, and plot histograms for each line of changes and unchanges in each direction


zM_backup = zM
zU_backup = zU

zM=zM_backup
zU=zU_backup

# make granges for zM and zU changing sites
zM.gr = makeGRangesFromDataFrame(zM, seqnames.field="seqnames", start.field="start", end.field="end")
zU.gr = makeGRangesFromDataFrame(zU, seqnames.field="seqnames", start.field="start", end.field="end")

# Take subset of zM and zU that represents gbM genes
gene_stats = read.table(file="gbm_tem_gene_calls.txt", sep="\t", header=TRUE)
zM = unique(zM[queryHits(findOverlaps(zM.gr, gene_ranges[gff.genes$gene_ID %in% gene_stats[gene_stats$gbm_calls==1,"gene_ID"]])),])
zU = unique(zU[queryHits(findOverlaps(zU.gr, gene_ranges[gff.genes$gene_ID %in% gene_stats[gene_stats$gbm_calls==1,"gene_ID"]])),])

# make granges for new versions of zM and zU changing sites
zM.gr = makeGRangesFromDataFrame(zM, seqnames.field="seqnames", start.field="start", end.field="end")
zU.gr = makeGRangesFromDataFrame(zU, seqnames.field="seqnames", start.field="start", end.field="end")



changing_sites_stats = NULL
bin_size=1
bin_extent=12000
for (this_line in gen3_30) {
  ##  line_m_to_u_call_ranges[[this_line]] = setdiff(line_m_to_u_call_ranges[[this_line]], DMRs.gr, ignore.strand=TRUE)
  ##  line_u_to_m_call_ranges[[this_line]] = setdiff(line_u_to_m_call_ranges[[this_line]], DMRs.gr, ignore.strand=TRUE)

  zM_line_changes = zM[queryHits(findOverlaps(zM.gr, line_m_to_u_call_ranges[[this_line]])),]
  cat(paste0("line: ", this_line, " no changes: ", length(line_m_to_u_call_ranges[[this_line]]), " overlaps: ", nrow(zM_line_changes), "\n"))
  zU_line_changes = zU[queryHits(findOverlaps(zU.gr, line_u_to_m_call_ranges[[this_line]])),]
  cat(paste0("line: ", this_line, " no changes: ", length(line_u_to_m_call_ranges[[this_line]]), " overlaps: ", nrow(zU_line_changes), "\n"))
  
  zM=merge(zM, zM_line_changes[,c(1,2,3,10)], by=c("seqnames","start","end"), all=TRUE)
  colnames(zM)[ncol(zM)]=this_line
  zU=merge(zU, zU_line_changes[,c(1,2,3,10)], by=c("seqnames","start","end"), all=TRUE)
  colnames(zU)[ncol(zU)]=this_line
  
  # pick out the zM sites corresponding to the changing sites for the line
  # pick out the zU sites corresponding to the changing sites for the line
  # add individual columns to zM and zU for each line
  
  #m_changed_site_counts = rbind.data.frame(m_changed_site_counts, cbind("line"=rep(this_line,800), "changed_site_counts"=hist(zM[zM[,ncol(zM)]==1,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts, "unchanged_site_counts"=hist(zM[is.na(zM[,ncol(zM)]),"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts))
  #u_changed_site_counts = rbind.data.frame(u_changed_site_counts, cbind("line"=rep(this_line,800), "changed_site_counts"=hist(zU[zU[,ncol(zU)]==1,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts, "unchanged_site_counts"=hist(zU[is.na(zU[,ncol(zU)]),"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts))
  changing_sites_stats = rbind.data.frame(changing_sites_stats, cbind("mid"=seq(bin_size/2,bin_extent,bin_size), "line"=rep(this_line,bin_extent/bin_size), "parental"=rep("M",bin_extent/bin_size), "nearest"=rep("M",bin_extent/bin_size), "changed_site_counts"=hist(zM[zM[,ncol(zM)]==1,"dist_to_nearest_M"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts, "unchanged_site_counts"=hist(zM[is.na(zM[,ncol(zM)]),"dist_to_nearest_M"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts))
  changing_sites_stats = rbind.data.frame(changing_sites_stats, cbind("mid"=seq(bin_size/2,bin_extent,bin_size), "line"=rep(this_line,bin_extent/bin_size), "parental"=rep("M",bin_extent/bin_size), "nearest"=rep("U",bin_extent/bin_size), "changed_site_counts"=hist(zM[zM[,ncol(zM)]==1,"dist_to_nearest_U"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts, "unchanged_site_counts"=hist(zM[is.na(zM[,ncol(zM)]),"dist_to_nearest_U"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts))
  changing_sites_stats = rbind.data.frame(changing_sites_stats, cbind("mid"=seq(bin_size/2,bin_extent,bin_size), "line"=rep(this_line,bin_extent/bin_size), "parental"=rep("U",bin_extent/bin_size), "nearest"=rep("M",bin_extent/bin_size), "changed_site_counts"=hist(zU[zU[,ncol(zU)]==1,"dist_to_nearest_M"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts, "unchanged_site_counts"=hist(zU[is.na(zU[,ncol(zU)]),"dist_to_nearest_M"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts))
  changing_sites_stats = rbind.data.frame(changing_sites_stats, cbind("mid"=seq(bin_size/2,bin_extent,bin_size), "line"=rep(this_line,bin_extent/bin_size), "parental"=rep("U",bin_extent/bin_size), "nearest"=rep("U",bin_extent/bin_size), "changed_site_counts"=hist(zU[zU[,ncol(zU)]==1,"dist_to_nearest_U"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts, "unchanged_site_counts"=hist(zU[is.na(zU[,ncol(zU)]),"dist_to_nearest_U"], breaks=seq(from=0, to=bin_extent, by=bin_size))$counts))
}
changing_sites_stats$mid = as.numeric(changing_sites_stats$mid)
changing_sites_stats$changed_site_counts = as.numeric(changing_sites_stats$changed_site_counts)
changing_sites_stats$unchanged_site_counts = as.numeric(changing_sites_stats$unchanged_site_counts)

changing_sites_stats$rate = changing_sites_stats$changed_site_counts/(changing_sites_stats$changed_site_counts+changing_sites_stats$unchanged_site_counts)

ggplot(changing_sites_stats, aes(x=mid, y=rate)) + geom_line(aes(colour=line)) + facet_wrap(nearest ~ parental) +xlim(0,500) +ylim(0,0.25)
#attempting to do fancy ribbon, but need to precompute min max and mean first to make this work:
#ggplot(changing_sites_stats, aes(x=mid, y=mean_rate, ymin=min_rate, ymax=max_rate)) + geom_line() +geom_ribbon() + facet_wrap(nearest ~ parental) +xlim(0,500) +ylim(0,0.25)

write.table(changing_sites_stats, file="changed_and_unchanged_by_line.txt", sep="\t", quote=FALSE, col.names =TRUE)





# save the new zM and zU instances
# export new histograms by line
# plot histograms with each line plotted separately

# gain rates
changed_site_counts = hist(zU[zU$u_to_m_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 

# loss rates
changed_site_counts = hist(zM[zM$m_to_u_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 



# Save things for Lizzie to use in comparisons
write.table(parental_consensus, file="parental_consensus.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
saveRDS(line_m_to_u_call_ranges,file="line_m_to_u_call_ranges.rds")
saveRDS(line_m_to_u_call_ranges,file="line_m_to_u_call_ranges.rds")

write.table(all_samples_meth_status[,c("Chromosome", "Locus", all_samples_meth_status$Locus+1)], file="CG_sites.bed", sep="\t", col.names=FALSE, row.names=FALSE)

levels(exon_ranges@seqnames@values) = mixedCaseChr(as.character(levels(exon_ranges@seqnames@values)))
exon_ranges@seqinfo@seqnames = levels(exon_ranges@seqnames@values)
levels(intron_ranges@seqnames@values) = mixedCaseChr(as.character(levels(intron_ranges@seqnames@values)))
intron_ranges@seqinfo@seqnames = levels(intron_ranges@seqnames@values)

#exon and intron length distributions:
relevant_exons = exon_ranges[queryHits(findOverlaps(exon_ranges, GBM_gene_loci.gr))]
ggplot(as.data.frame(relevant_exons@ranges@width)) + geom_histogram(aes(x=relevant_exons@ranges@width), binwidth = 1) +xlim(0,5000)
write.table(as.data.frame(relevant_exons@ranges@width), file="gbM_gene_exon_lengths.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
relevant_introns = intron_ranges[queryHits(findOverlaps(intron_ranges, GBM_gene_loci.gr))]
ggplot(as.data.frame(relevant_introns@ranges@width)) + geom_histogram(aes(x=relevant_introns@ranges@width), binwidth = 1) +xlim(0,5000)
write.table(as.data.frame(relevant_introns@ranges@width), file="gbM_gene_intron_lengths.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

#exon and intron CG site densities:
levels(CG_site_ranges@seqnames@values) = mixedCaseChr(as.character(levels(CG_site_ranges@seqnames@values)))
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
exon_site_counts = table(queryHits(findOverlaps(relevant_exons, CG_site_ranges)))
exon_CG_sites = rep(0, length(relevant_exons))
exon_CG_sites[as.data.frame(exon_site_counts)$Var1]=as.data.frame(exon_site_counts)$Freq
ggplot(as.data.frame(cbind(relevant_exons@ranges@width, exon_CG_sites))) + geom_point(aes(x=V1, y=exon_CG_sites))
ggplot(as.data.frame(cbind(relevant_exons@ranges@width, exon_CG_sites))) + geom_density_2d(aes(x=V1, y=exon_CG_sites))
smoothScatter(relevant_exons@ranges@width, exon_CG_sites)

cor(relevant_exons@ranges@width, exon_CG_sites)

intron_site_counts = table(queryHits(findOverlaps(relevant_introns, CG_site_ranges)))
intron_CG_sites = rep(0, length(relevant_introns))
intron_CG_sites[as.data.frame(intron_site_counts)$Var1]=as.data.frame(intron_site_counts)$Freq
ggplot(as.data.frame(cbind(relevant_introns@ranges@width, intron_CG_sites))) + geom_point(aes(x=V1, y=intron_CG_sites))
ggplot(as.data.frame(cbind(relevant_introns@ranges@width, intron_CG_sites))) + geom_density_2d(aes(x=V1, y=intron_CG_sites))
smoothScatter(relevant_introns@ranges@width, intron_CG_sites)

cor(relevant_introns@ranges@width, intron_CG_sites)

sum(relevant_exons@ranges@width)/sum(exon_CG_sites)
#[1] 34.74018
sum(relevant_introns@ranges@width)/sum(intron_CG_sites)
#[1] 59.30603

  

# distribution of nearest sites

pdf("CG_site_densities_GBM_TEM_UMR_all.pdf")

target_segments.gr = GBM_segments.gr
target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_segments.gr))]

nearest_sites = nearest(target_sites)
dist_to_nearest = distance(target_sites, target_sites[nearest_sites])+1

print(ggplot(as.data.frame(dist_to_nearest)) + geom_histogram(aes(x=dist_to_nearest), binwidth = 1) +xlim(0,100))

target_segments.gr = TEM_segments.gr
target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_segments.gr))]

nearest_sites = nearest(target_sites)
dist_to_nearest = distance(target_sites, target_sites[nearest_sites])+1

print(ggplot(as.data.frame(dist_to_nearest)) + geom_histogram(aes(x=dist_to_nearest), binwidth = 1) +xlim(0,100))

target_segments.gr = UMR_segments.gr
target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_segments.gr))]

nearest_sites = nearest(target_sites)
dist_to_nearest = distance(target_sites, target_sites[nearest_sites])+1

print(ggplot(as.data.frame(dist_to_nearest)) + geom_histogram(aes(x=dist_to_nearest), binwidth = 1) +xlim(0,100))

target_sites = CG_site_ranges

nearest_sites = nearest(target_sites)
dist_to_nearest = distance(target_sites, target_sites[nearest_sites])+1

print(ggplot(as.data.frame(dist_to_nearest)) + geom_histogram(aes(x=dist_to_nearest), binwidth = 1) +xlim(0,100))

dev.off()


# Find general density of sites, nearest or not, in GBM segments and others
# Identify GenomicRegions covering GBM range of such GBM genes (from first to last methylated base in gen 3)
target_ranges = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])))]

target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
#z1=cbind(as.data.frame(target_sites), dist_to_nearest = rep(NA,length(target_sites)))
z2 = c()

site_count = 0
# for each range:
for(i in 1:length(target_ranges)) {
  target_range = target_ranges[i]
  # find overlapping sites
  near_sites = target_sites[queryHits(findOverlaps(target_sites, target_range))]
  if (length(near_sites) > 1) {
    # for each site in range
    for(j in 1:length(near_sites)-1) {
      site_count = site_count + 1
      z2 = c(z2, abs(as.data.frame(near_sites[j])$start-as.data.frame(near_sites[j+1])$start)) 
    }
  }
}

pdf(paste0(project_id,"_dist_between_CG_sites_in_GBM_segments.pdf"))
# plot distribution of gap sizes between adjacent sites
print(ggplot(as.data.frame(z2)) + geom_histogram(aes(x=z2), binwidth = 1) + xlim(c(0,500)))
dev.off()

# Find general density of sites, nearest or not, in TEM segments
target_ranges = TEM_segments.gr

target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
#z1=cbind(as.data.frame(target_sites), dist_to_nearest = rep(NA,length(target_sites)))
z2 = c()

site_count = 0
# for each range:
for(i in 1:length(target_ranges)) {
  target_range = target_ranges[i]
  # find overlapping sites
  near_sites = target_sites[queryHits(findOverlaps(target_sites, target_range))]
  if (length(near_sites) > 1) {
    # for each site in range
    for(j in 1:length(near_sites)-1) {
      site_count = site_count + 1
      z2 = c(z2, abs(as.data.frame(near_sites[j])$start-as.data.frame(near_sites[j+1])$start)) 
    }
  }
}

pdf(paste0(project_id,"_dist_between_CG_sites_in_TEM_segments.pdf"))
# plot distribution of gap sizes between adjacent sites
print(ggplot(as.data.frame(z2)) + geom_histogram(aes(x=z2), binwidth = 1) + xlim(c(0,500)))
dev.off()

# Find general density of sites, nearest or not, in UMR segments
target_ranges = UMR_segments.gr

target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
#z1=cbind(as.data.frame(target_sites), dist_to_nearest = rep(NA,length(target_sites)))
z2 = c()

site_count = 0
# for each range:
for(i in 1:length(target_ranges)) {
  target_range = target_ranges[i]
  # find overlapping sites
  near_sites = target_sites[queryHits(findOverlaps(target_sites, target_range))]
  if (length(near_sites) > 1) {
    # for each site in range
    for(j in 1:length(near_sites)-1) {
      site_count = site_count + 1
      z2 = c(z2, abs(as.data.frame(near_sites[j])$start-as.data.frame(near_sites[j+1])$start)) 
    }
  }
}

pdf(paste0(project_id,"_dist_between_CG_sites_in_UMR_segments.pdf"))
# plot distribution of gap sizes between adjacent sites
print(ggplot(as.data.frame(z2)) + geom_histogram(aes(x=z2), binwidth = 1) + xlim(c(0,500)))
dev.off()



# Get list of gbM segments, and which gene they are in for Amy
gbM_segment_gene_loci = cbind(gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5","gene_ID")][queryHits(findOverlaps(GBM_gene_loci.gr, GBM_segments.gr)),], as.data.frame(GBM_segments.gr)[subjectHits(findOverlaps(GBM_gene_loci.gr, GBM_segments.gr)),c("start", "end")])
colnames(gbM_segment_gene_loci) = c("Chromosome", "gene_start", "gene_end", "gene_ID", "segment_start", "segment_end")
write.table(gbM_segment_gene_loci, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="gbM_segment_gene_loci.txt")
     

GBM_gene_loci=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")]
colnames(GBM_gene_loci)=c("Chromosome", "start_site", "end_site")
#GBM_gene_loci$Chromosome=mixedCaseChr(GBM_gene_loci$Chromosome)
#GBM_gene_loci.gr=reduce(makeGRangesFromDataFrame(df=GBM_gene_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome"))
GBM_gene_loci.gr=makeGRangesFromDataFrame(df=GBM_gene_loci, start.field="start_site", end.field="end_site", seqnames.field="Chromosome")




parental_state = parental_consensus[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
nn_parental_state = parental_consensus[queryHits(findOverlaps(CG_site_ranges, target_ranges))][nearest_sites]
type_of_pair = ifelse(parental_state=="M" & nn_parental_state=="M", "MM", ifelse(parental_state=="U" & nn_parental_state=="U", "UU", ifelse(parental_state=="M" & nn_parental_state=="U", "MU", ifelse(parental_state=="U" & nn_parental_state=="M", "MU", "XX"))))
z1=cbind(as.data.frame(target_sites), dist_to_nearest, type_of_pair)
