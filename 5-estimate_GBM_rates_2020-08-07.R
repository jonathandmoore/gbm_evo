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
if(!require(mixtools)){
  install.packages("mixtools")
  library(mixtools)
}
if(!require(rootSolve)){
  install.packages("rootSolve")
  library(rootSolve)
}
install.packages("ggpubr")


#### This part sets up libraries and contextual information and metadata about the project

library(data.table)
library(stringr)
library(ggplot2)
library(plyr)
library(ggpubr)


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
  gen3_30=c("Line4Gen3","Line8Gen3","Line29Gen31","Line39Gen31","Line49Gen31","Line59Gen31","Line79Gen31","Line89Gen31","Line99Gen31","Line109Gen31","Line119Gen31")
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
num_lines_varying_from_parentals = NA
num_lines_m_to_u_from_parentals = NA
num_lines_u_to_m_from_parentals = NA
num_lines_m_or_u_excl_parentals = NA
for(line_gen in gen3_30) {
  line_gen_no = line_gen_no + 1
  if(opt$verbose) {cat(paste0("Checking variation in methylation calls from parental consensus at Gen 30: ",line_gen,"\n"))}
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
    #first_line_gen=FALSE
  }
}



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
install.packages("sqldf")
library(sqldf)
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
m_to_u_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_m_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
# Make a GRanges for U->M sites
u_to_m_call_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status[site_varies_from_parentals & (no_u_3>=2),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
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
gff.genes$average_methylation=table(gff.genes[subjectHits(olaps_CG_sites),]$gene_ID)[gff.genes$gene_ID]

# Find out which transposons are really transposons, and which are 'euchromatic transposons' as per my own method.
gff.transposons$average_methylation=table(gff.transposons[subjectHits(tolaps_CG_sites),]$gene_ID)[gff.transposons$gene_ID]


# Create a data frame with gene ID and average methylation for each overlap between a site and a gene model
olaps_average_methylation=data.frame(cbind("gene_ID"=gff.genes[subjectHits(olaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(olaps_CG_sites)]))
tolaps_average_methylation=data.frame(cbind("gene_ID"=gff.transposons[subjectHits(tolaps_CG_sites),"gene_ID"],"site_average_methylation"=(average_methylation[,1]/average_methylation[,2])[queryHits(tolaps_CG_sites)]))

# R 3 version:
# Convert average values to numerics
#olaps_average_methylation$site_average_methylation=as.numeric(levels(olaps_average_methylation$site_average_methylation)[olaps_average_methylation$site_average_methylation])
#tolaps_average_methylation$site_average_methylation=as.numeric(levels(tolaps_average_methylation$site_average_methylation)[tolaps_average_methylation$site_average_methylation])

# R 4 version:
# Convert average values to numerics
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



# Plot densities of loss and gain rates across all gbM genes

library(ggplot2)
pdf(paste0(project_id,"_density_of_loss_and_gain_rate_across_genes.pdf"))

print(ggplot(gff.genes[gff.genes$m_class=="Gene-body Methylated",]) + geom_histogram(aes(x=m_to_u_count/(all_M_count+m_to_u_count)), bins=100))

print(ggplot(gff.genes[gff.genes$m_class=="Gene-body Methylated",]) + geom_histogram(aes(x=u_to_m_count/(all_U_count+u_to_m_count)), bins=100))

dev.off() 
 
loss_rates = gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count/(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count+gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_M_count)
gain_rates = gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count/(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count+gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_U_count)

# These are mean rates across any/all lines
mean(loss_rates, na.rm=TRUE)
#Schmitz [1] 0.3317303   Becker [1] 0.5405472

mean(gain_rates, na.rm=TRUE)
#Schmitz [1] 0.06574726  Becker [1] 0.08396383

no_gen=30
cc_per_gen=34
no_lines = 4  # 9
mean(loss_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)
#Schmitz [1] 9.214732e-05   Becker [1] [1] 5.888314e-05

mean(gain_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)
#Schmitz [1] 1.826313e-05  Becker [1] [1] 9.146387e-06

# Overall average loss and gain rates
sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_M_count))
sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count, na.rm=TRUE)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count, na.rm=TRUE)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_M_count, na.rm=TRUE))
#Schmitz [1] 0.2187601  Becker [1] 0.4183502

sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_U_count))
sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count, na.rm=TRUE)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count, na.rm=TRUE)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_U_count, na.rm=TRUE))
#Schmitz [1] 0.05959142  Becker [1] 0.06611501

# translated to per cell cycle rates:
sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_M_count))/(no_gen*cc_per_gen*no_lines)
sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count, na.rm=TRUE)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$m_to_u_count, na.rm=TRUE)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_M_count, na.rm=TRUE))/(no_gen*cc_per_gen*no_lines)
#Schmitz [1] 6.07667e-05  Becker [1] [1] 4.557192e-05

sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_U_count))/(no_gen*cc_per_gen*no_lines)
sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count, na.rm=TRUE)/(sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$u_to_m_count, na.rm=TRUE)+sum(gff.genes[gff.genes$m_class=="Gene-body Methylated",]$all_U_count, na.rm=TRUE))/(no_gen*cc_per_gen*no_lines)
#Schmitz [1] 1.655317e-05  Becker [1] [1] 7.202071e-06



gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE))
gene_body_loci.gr=makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)

unmethylated_gene_loci.gr=makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)


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

  this_sample="parental_consensus"

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
  this_meth = all_samples_meth_status[parental_consensus=="M",c("Chromosome", "Locus")]
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
  
    gene_gbm_overlap = as.data.frame(gbm_genes)
    gene_tem_overlap = as.data.frame(tem_genes)
    gene_umr_overlap = as.data.frame(umr_genes)
	colnames(gene_gbm_overlap) = c("gene_ID", this_sample)
	colnames(gene_tem_overlap) = c("gene_ID", this_sample)
	colnames(gene_umr_overlap) = c("gene_ID", this_sample)
    gene_gbm_mCG_sites = this_gene_GBM_mCG_sites
    gene_tem_mCG_sites = this_gene_TEM_mCG_sites

gene_gbm_overlap[is.na(gene_gbm_overlap)] = 0
gene_tem_overlap[is.na(gene_tem_overlap)] = 0
gene_umr_overlap[is.na(gene_umr_overlap)] = 0
gene_gbm_overlap[,2] = as.numeric(gene_gbm_overlap[,2])
gene_tem_overlap[,2] = as.numeric(gene_tem_overlap[,2])
gene_umr_overlap[,2] = as.numeric(gene_umr_overlap[,2])
gene_gbm_mCG_sites[is.na(gene_gbm_mCG_sites)] = 0
gene_tem_mCG_sites[is.na(gene_tem_mCG_sites)] = 0

gene_gbm_mCG_sites$gene_ID = as.character(gene_gbm_mCG_sites$gene_ID)
gene_tem_mCG_sites$gene_ID = as.character(gene_tem_mCG_sites$gene_ID)

# work out how many genes are both gbM and TEM

gene_both_overlap = NULL
gene_both_overlap_filtered = NULL
gene_both_overlap_resolved = NULL
gene_both_overlap_resolved_gbm = NULL
gene_both_overlap_resolved_tem = NULL
gene_gbm_calls = NULL
gene_tem_calls = NULL
no_samples = 0

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
  
    gene_both_overlap = as.data.frame(gene_stats[,c(1,6)])
    gene_both_overlap_filtered = as.data.frame(gene_stats[,c(1,7)])
    gene_both_overlap_resolved = as.data.frame(gene_stats[,c(1,8)])
    gene_both_overlap_resolved_gbm = as.data.frame(gene_stats[,c(1,9)])
    gene_both_overlap_resolved_tem = as.data.frame(gene_stats[,c(1,10)])
    gene_gbm_calls = as.data.frame(gene_stats[,c(1,11)])
    gene_tem_calls = as.data.frame(gene_stats[,c(1,12)])

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
  

# plot distributions of numbers of gbM genes, TEM genes and both per accession
gbm_gene_counts = colSums(gene_gbm_overlap[2], na.rm=TRUE)
tem_gene_counts = colSums(gene_tem_overlap[2], na.rm=TRUE)
both_counts = colSums(gene_both_overlap[2], na.rm=TRUE)
both_filtered_counts = colSums(gene_both_overlap_filtered[2], na.rm=TRUE)
both_resolved_counts = colSums(gene_both_overlap_resolved[2], na.rm=TRUE)
both_resolved_gbm_counts = colSums(gene_both_overlap_resolved_gbm[2], na.rm=TRUE)
both_resolved_tem_counts = colSums(gene_both_overlap_resolved_tem[2], na.rm=TRUE)

write.table(gene_stats, file="gbm_tem_gene_calls.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# gene_stats = read.table(file="gbm_tem_gene_calls.txt", sep="\t", header=TRUE)

gbm_gene_counts
#parental_consensus 
#             15241 
tem_gene_counts
#parental_consensus 
#              3501 
both_counts
#parental_consensus 
#               702 
both_filtered_counts
#parental_consensus 
#               187 
both_resolved_counts
#parental_consensus 
#               100 
both_resolved_gbm_counts
#parental_consensus 
#                42 
both_resolved_tem_counts
#parental_consensus 
#                58 

colSums(gene_stats[,2:12])
#Schmitz
#  gene_gbm_overlap   gene_tem_overlap gene_gbm_mCG_sites gene_tem_mCG_sites       both_overlap           filtered           resolved       resolved_gbm       resolved_tem          gbm_calls          tem_calls 
#             15241               3501             218741              51840                702                187                100                 42                 58              14581               2857 
#Becker
#  gene_gbm_overlap   gene_tem_overlap gene_gbm_mCG_sites gene_tem_mCG_sites       both_overlap           filtered           resolved       resolved_gbm       resolved_tem          gbm_calls          tem_calls 
#             15241               3501             115513              31052                702                 40                 27                 10                 17              14549               2816 




# New definition of gene_body_loci based on overlaps with gbM segments, but no significant TEM segments

gene_body_loci.gr=makeGRangesFromDataFrame(df=gff.genes[gff.genes$gene_ID %in% gene_gbm_calls$gene_ID, c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)

# The above definition is wrong (but was used to create Schmitz zM and zU so kept here in code.
# It should have excluded entries in gene_gbm_calls where parental_consensus=0
gene_body_loci.gr=makeGRangesFromDataFrame(df=gff.genes[gff.genes$gene_ID %in% gene_gbm_calls[gene_gbm_calls$parental_consensus==1,]$gene_ID, c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1", ignore.strand=TRUE)


# Three different versions of target_ranges:

# Identify GenomicRegions covering GBM genes which don't have a TE segment overlapping
###  WE DON'T DO THIS ANY MORE AS WE ARE ALLOWING A SMALL TEM SEGMENT OVERLAP AS LONG AS GBM SEGMENT IS MORE THAN 4x AS LONG
#target_ranges = setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])
#target_ranges = gene_body_loci.gr
target_ranges = reduce(gene_body_loci.gr)

### should we exclude gene_body_loci.gr if they only have 1 or 2 GBM sites?
# no because then we can't look at cooperative maintenance which occurs in sparsely methylated zones


# Identify GenomicRegions covering 'body' of such GBM genes (excluding 2kb 5' and 1kb 3')
#target_ranges = 

# Identify GenomicRegions covering GBM range of such GBM genes (from first to last methylated base in gen 3)
#target_ranges = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])))]

# Identify gbM segments
#target_ranges = GBM_segments.gr



target_stats = merge(table(subjectHits(findOverlaps(U_parent_call_ranges, target_ranges))), table(subjectHits(findOverlaps(M_parent_call_ranges, target_ranges))), by="Var1", all=TRUE)
colnames(target_stats) = c("Var1", "U_count", "M_count")
target_stats = merge(target_stats, table(subjectHits(findOverlaps(m_to_u_call_ranges, target_ranges))), by="Var1", all=TRUE)
colnames(target_stats)[4] = "m_to_u_count"

target_stats = merge(target_stats,  table(subjectHits(findOverlaps(u_to_m_call_ranges, target_ranges))), by="Var1", all=TRUE)
colnames(target_stats)[5] = "u_to_m_count"

target_stats[is.na(target_stats$M_count),"M_count"] = 0
target_stats[is.na(target_stats$U_count),"U_count"] = 0
target_stats[is.na(target_stats$m_to_u_count),"m_to_u_count"] = 0
target_stats[is.na(target_stats$u_to_m_count),"u_to_m_count"] = 0

# Plot densities of loss and gain rates across all target ranges

#library(ggplot2)
pdf(paste0(project_id,"_density_of_loss_and_gain_rate_across_gbM_segments.pdf"))

print(ggplot(target_stats) + geom_histogram(aes(x=m_to_u_count/M_count), bins=100))

print(ggplot(target_stats) + geom_histogram(aes(x=u_to_m_count/U_count), bins=100))

dev.off() 
 
loss_rates = target_stats$m_to_u_count/target_stats$M_count
gain_rates = target_stats$u_to_m_count/target_stats$U_count

# These are mean rates across any/all lines
mean(loss_rates, na.rm=TRUE)
# old gbm genes definition: [1] 0.2318511
# new gbm genes definition: [1] 0.180847
# loss rates have nearly halved

mean(gain_rates, na.rm=TRUE)
# old gbm genes definition: [1] 0.1749956
# new gbm genes definition: [1] 0.04021518
# gain rates have gone down by factor of 4


no_gen=30
cc_per_gen=34
no_lines = 4 # 9
mean(loss_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)
#old[1] 6.440308e-05
#new 4.432525e-05
#Becker [1] 3.415943e-05

mean(gain_rates, na.rm=TRUE)/(no_gen*cc_per_gen*no_lines)
#old[1] 4.86099e-05
#new 9.856662e-06
#Becker [1] 6.460099e-06

# Overall average loss and gain rates
sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))
#[1] 0.1487853
#[1] 0.12191
#Becker [1] 0.231335

sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))
#[1] 0.1412831
#[1] 0.034488
#Becker [1] 0.0520796

# translated to per cell cycle rates:
sum(target_stats$m_to_u_count)/(sum(target_stats$M_count))/(no_gen*cc_per_gen*no_lines)
#[1] 4.132925e-05
#[1] 2.987989e-05
#Becker [1] 2.519989e-05

sum(target_stats$u_to_m_count)/(sum(target_stats$U_count))/(no_gen*cc_per_gen*no_lines)
#[1] 3.924532e-05
#[1] 8.452942e-06
#Becker [1] 5.673159e-06
  
# Find U sites in target_ranges 
target_U_sites = U_parent_call_ranges[queryHits(findOverlaps(U_parent_call_ranges, target_ranges))]

# Find M sites in target_ranges
target_M_sites = M_parent_call_ranges[queryHits(findOverlaps(M_parent_call_ranges, target_ranges))]

# Find m_to_u_sites in target_ranges
#target_m_to_u_sites = m_to_u_call_ranges[queryHits(findOverlaps(m_to_u_call_ranges, target_ranges))]

# Find u_to_m_sites in target_ranges
#target_u_to_m_sites = u_to_m_call_ranges[queryHits(findOverlaps(u_to_m_call_ranges, target_ranges))]


#z=cbind(as.data.frame(target_U_sites), dist_to_nearest_M = rep(NA,length(target_U_sites)), u_to_m_count = rep(NA,length(target_U_sites)))
zU=cbind(as.data.frame(target_U_sites), dist_to_nearest_M = rep(NA,length(target_U_sites)), sites_to_nearest_M = rep(NA,length(target_U_sites)), dist_to_nearest_U = rep(NA,length(target_U_sites)), sites_to_nearest_U = rep(NA,length(target_U_sites)), u_to_m_count = rep(NA,length(target_U_sites)))


###These loops take ages (a day or two)
# for each U site in target_ranges:
i = 0
for (k in 1:length(target_ranges)) {
  cat(k)
  target_range = target_ranges[k]
  this_target_U_sites = target_U_sites[queryHits(findOverlaps(target_U_sites, target_range))]
  if (length(this_target_U_sites) > 0) {
    local_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_range))]
    for (j in 1:length(this_target_U_sites)) {
      target_site = this_target_U_sites[j]
      i = i + 1
      
      #for(i in 1:length(target_U_sites)) {
      #for(target_site in target_U_sites) {
      #   Find the target_range it is part of
      #target_range = target_ranges[subjectHits(findOverlaps(target_site, target_ranges))]
      #  target_range = target_ranges[subjectHits(findOverlaps(target_U_sites[i], target_ranges))]
      #   Find nearest M site in target_range, capture distance in bp and No. sites
      near_Ms = M_parent_call_ranges[queryHits(findOverlaps(M_parent_call_ranges, target_range))]
      if (length(near_Ms) > 0) {
        # Find the nearest in bp
        dist_to_Ms = abs(as.data.frame(near_Ms)$start - as.data.frame(target_site)$start)
        #dist_to_Ms = abs(as.data.frame(near_Ms)$start-as.data.frame(target_sites[i])$start)
        min_dist = dist_to_Ms[dist_to_Ms == min(dist_to_Ms)][1] # [1] takes care of case if there are two equally near
        # Find the nearest in sites
        #sites_to_Ms = abs(queryHits(findOverlaps(CG_site_ranges,near_Ms))-queryHits(findOverlaps(CG_site_ranges,target_U_sites[i])))
        sites_to_Ms = abs(queryHits(findOverlaps(local_sites, near_Ms)) - queryHits(findOverlaps(local_sites, target_site)))
        min_sites = sites_to_Ms[sites_to_Ms == min(sites_to_Ms)][1] # [1] takes care of case if there are two equally near
      } else {
        # There are no M sites nearby
        min_dist = NA
		min_sites = NA
      }
      zU[i, "dist_to_nearest_M"] = min_dist
      zU[i, "sites_to_nearest_M"] = min_sites
      #   Find number of times target_site is found M in gen 30 lines
      zU[i, "u_to_m_count"] = length(findOverlaps(target_site, u_to_m_call_ranges))
      
      #   Find nearest U site in target_range, capture distance in bp and No. sites
      # exclude self using setdiff
      #non_self_U_parent_ranges = setdiff(U_parent_call_ranges, target_U_sites[i])
      #near_Us = non_self_U_parent_ranges[queryHits(findOverlaps(non_self_U_parent_ranges, target_range))]
      near_Us = setdiff(U_parent_call_ranges[queryHits(findOverlaps(U_parent_call_ranges, target_range))], target_site)
      if (length(near_Us) > 0) {
        # Find the nearest in bp
        #dist_to_Ms = abs(as.data.frame(near_Ms)$start-as.data.frame(target_site)$start)
        dist_to_Us = abs(as.data.frame(near_Us)$start - as.data.frame(target_site)$start)
        min_dist = dist_to_Us[dist_to_Us == min(dist_to_Us)][1] # [1] takes care of case if there are two equally near
        # Find the nearest in sites
        sites_to_Us = abs(queryHits(findOverlaps(local_sites, near_Us)) - queryHits(findOverlaps(local_sites, target_U_sites[i])))
        min_sites = sites_to_Us[sites_to_Us == min(sites_to_Us)][1] # [1] takes care of case if there are two equally near
      } else {
        # There are no U sites nearby
        min_dist = NA
		min_sites = NA
      }
      zU[i, "dist_to_nearest_U"] = min_dist
      zU[i, "sites_to_nearest_U"] = min_sites
      #   Find number of times target_site is found U in gen 30 lines
      #zU[i, "m_to_u_count"] = length(findOverlaps(target_U_sites[i], m_to_u_call_ranges))
    }
  }
}

# Write out zU to save for later
write.table(zU, file=paste0(project_id,"_dist_from_U_to_nearest_X.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
#zU=read.table(paste0(project_id,"_dist_from_U_to_nearest_X.tsv"),sep="\t",header=TRUE)

zM = cbind(
  as.data.frame(target_M_sites),
  dist_to_nearest_M = rep(NA, length(target_M_sites)),
  sites_to_nearest_M = rep(NA, length(target_M_sites)),
  dist_to_nearest_U = rep(NA, length(target_M_sites)),
  sites_to_nearest_U = rep(NA, length(target_M_sites)),
  m_to_u_count = rep(NA, length(target_M_sites))
)

i = 0
# for each M site in target_ranges:
#for(i in 1:length(target_M_sites)) {
for (k in 1:length(target_ranges)) {
  target_range = target_ranges[k]
  this_target_M_sites = target_M_sites[queryHits(findOverlaps(target_M_sites, target_range))]
  if (length(this_target_M_sites) > 0) {
    local_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_range))]
    for (j in 1:length(this_target_M_sites)) {
      target_site = this_target_M_sites[j]
      i = i + 1
      
      #for(target_site in target_U_sites) {
      #   Find the target_range it is part of
      #target_range = target_ranges[subjectHits(findOverlaps(target_site, target_ranges))]
      target_range = target_ranges[subjectHits(findOverlaps(target_site, target_ranges))]
      #   Find nearest M site in target_range, capture distance in bp and No. sites
      # exclude self using setdiff
      #non_self_M_parent_ranges = setdiff(M_parent_call_ranges, target_M_sites[i])
      #near_Ms = non_self_M_parent_ranges[queryHits(findOverlaps(non_self_M_parent_ranges, target_range))]
      near_Ms = setdiff(M_parent_call_ranges[queryHits(findOverlaps(M_parent_call_ranges, target_range))], target_site)
      if (length(near_Ms) > 0) {
        # Find the nearest in bp
        #dist_to_Ms = abs(as.data.frame(near_Ms)$start-as.data.frame(target_site)$start)
        dist_to_Ms = abs(as.data.frame(near_Ms)$start - as.data.frame(target_site)$start)
        min_dist = dist_to_Ms[dist_to_Ms == min(dist_to_Ms)][1] # [1] takes care of case if there are two equally near
        # Find the nearest in sites
        sites_to_Ms = abs(queryHits(findOverlaps(local_sites, near_Ms)) - queryHits(findOverlaps(local_sites, target_site)))
        min_sites = sites_to_Ms[sites_to_Ms == min(sites_to_Ms)][1] # [1] takes care of case if there are two equally near
        
      } else {
        # There are no M sites nearby
        min_dist = NA
		min_sites = NA
      }
      zM[i, "dist_to_nearest_M"] = min_dist
      zM[i, "sites_to_nearest_M"] = min_sites
      #   Find number of times target_site is found M in gen 30 lines
      #zM[i, "u_to_m_count"] = length(findOverlaps(target_M_sites[i], u_to_m_call_ranges))
      
      #   Find nearest U site in target_range, capture distance in bp and No. sites
      near_Us = U_parent_call_ranges[queryHits(findOverlaps(U_parent_call_ranges, target_range))]
      if (length(near_Us) > 0) {
        # Find the nearest in bp
        #dist_to_Ms = abs(as.data.frame(near_Ms)$start-as.data.frame(target_site)$start)
        dist_to_Us = abs(as.data.frame(near_Us)$start - as.data.frame(target_site)$start)
        min_dist = dist_to_Us[dist_to_Us == min(dist_to_Us)][1] # [1] takes care of case if there are two equally near
        # Find the nearest in sites
        sites_to_Us = abs(queryHits(findOverlaps(local_sites, near_Us)) - queryHits(findOverlaps(local_sites, target_site)))
        min_sites = sites_to_Us[sites_to_Us == min(sites_to_Us)][1] # [1] takes care of case if there are two equally near
      } else {
        # There are no U sites nearby
        min_dist = NA
		min_sites = NA
      }
      zM[i, "dist_to_nearest_U"] = min_dist
      zM[i, "sites_to_nearest_U"] = min_sites
      #   Find number of times target_site is found U in gen 30 lines
      zM[i, "m_to_u_count"] = length(findOverlaps(target_site, m_to_u_call_ranges))
    }
  }
}

# Write out zM to save for later
write.table(zM, file=paste0(project_id,"_dist_from_M_to_nearest_X.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
#zM=read.table(paste0(project_id,"_dist_from_M_to_nearest_X.tsv"),sep="\t",header=TRUE)


library(ggplot2)
pdf(paste0(project_id,"_dist_to_nearest_X_all_X_sites_in_GBM_genes_bp_new.pdf"))

# plot changing and unchanging sites separately
print(ggplot(zU) + geom_histogram(aes(x=dist_to_nearest_M), binwidth = 10) + facet_wrap(~u_to_m_count) +xlim(c(0,2000)))

# plot changing and unchanging sites together
print(ggplot(zU) + geom_histogram(aes(x=dist_to_nearest_M), binwidth = 10) +xlim(c(0,2000)))

changed_site_counts = hist(zU[zU$u_to_m_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 

# plot changing and unchanging sites separately
print(ggplot(zU) + geom_histogram(aes(x=dist_to_nearest_U), binwidth = 10) + facet_wrap(~u_to_m_count) +xlim(c(0,2000)))

# plot changing and unchanging sites together
print(ggplot(zU) + geom_histogram(aes(x=dist_to_nearest_U), binwidth = 10) +xlim(c(0,2000)))

changed_site_counts = hist(zU[zU$u_to_m_count==1,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 


# plot changing and unchanging sites separately
print(ggplot(zM) + geom_histogram(aes(x=dist_to_nearest_M), binwidth = 10) + facet_wrap(~m_to_u_count) +xlim(c(0,2000)))

# plot changing and unchanging sites together
print(ggplot(zM) + geom_histogram(aes(x=dist_to_nearest_M), binwidth = 10) +xlim(c(0,2000)))

changed_site_counts = hist(zM[zM$m_to_u_count==1,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"dist_to_nearest_M"], breaks=seq(from=0, to=8000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 

# plot changing and unchanging sites separately
print(ggplot(zM) + geom_histogram(aes(x=dist_to_nearest_U), binwidth = 10) + facet_wrap(~m_to_u_count) +xlim(c(0,2000)))

# plot changing and unchanging sites together
print(ggplot(zM) + geom_histogram(aes(x=dist_to_nearest_U), binwidth = 10) +xlim(c(0,2000)))

changed_site_counts = hist(zM[zM$m_to_u_count==1,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"dist_to_nearest_U"], breaks=seq(from=0, to=8000, by=10))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=5, to=7995, by=10), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,1000))) 


dev.off()


# Plot gaps in sites rather than bp

pdf(paste0(project_id,"_sites_to_nearest_X_all_X_sites_in_GBM_genes_sites_new.pdf"))

# plot changing and unchanging sites separately
print(ggplot(zU) + geom_histogram(aes(x=sites_to_nearest_M), binwidth = 1) + facet_wrap(~u_to_m_count) +xlim(c(0,200)))

# plot changing and unchanging sites together
print(ggplot(zU) + geom_histogram(aes(x=sites_to_nearest_M), binwidth = 1) +xlim(c(0,200)))

changed_site_counts = hist(zU[zU$u_to_m_count==1,"sites_to_nearest_M"], breaks=seq(from=-0.5, to=200.5, by=1))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"sites_to_nearest_M"], breaks=seq(from=-0.5, to=200.5, by=1))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=-0.5, to=200.5, by=1), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,100))) 

# plot changing and unchanging sites separately
print(ggplot(zU) + geom_histogram(aes(x=sites_to_nearest_U), binwidth = 1) + facet_wrap(~u_to_m_count) +xlim(c(0,200)))

# plot changing and unchanging sites together
print(ggplot(zU) + geom_histogram(aes(x=sites_to_nearest_U), binwidth = 1) +xlim(c(0,200)))

changed_site_counts = hist(zU[zU$u_to_m_count==1,"sites_to_nearest_U"], breaks=seq(from=-0.5, to=200.5, by=1))$counts
unchanged_site_counts = hist(zU[zU$u_to_m_count==0,"sites_to_nearest_U"], breaks=seq(from=-0.5, to=200.5, by=1))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=-0.5, to=200.5, by=1), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,100))) 


# plot changing and unchanging sites separately
print(ggplot(zM) + geom_histogram(aes(x=sites_to_nearest_M), binwidth = 1) + facet_wrap(~m_to_u_count) +xlim(c(0,200)))

# plot changing and unchanging sites together
print(ggplot(zM) + geom_histogram(aes(x=sites_to_nearest_M), binwidth = 1) +xlim(c(0,200)))

changed_site_counts = hist(zM[zM$m_to_u_count==1,"sites_to_nearest_M"], breaks=seq(from=-0.5, to=200.5, by=1))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"sites_to_nearest_M"], breaks=seq(from=-0.5, to=200.5, by=1))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=-0.5, to=200.5, by=1), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,100))) 

# plot changing and unchanging sites separately
print(ggplot(zM) + geom_histogram(aes(x=sites_to_nearest_U), binwidth = 1) + facet_wrap(~m_to_u_count) +xlim(c(0,200)))

# plot changing and unchanging sites together
print(ggplot(zM) + geom_histogram(aes(x=sites_to_nearest_U), binwidth = 1) +xlim(c(0,200)))

changed_site_counts = hist(zM[zM$m_to_u_count==1,"sites_to_nearest_U"], breaks=seq(from=-0.5, to=200.5, by=1))$counts
unchanged_site_counts = hist(zM[zM$m_to_u_count==0,"sites_to_nearest_U"], breaks=seq(from=-0.5, to=200.5, by=1))$counts

# plot rate of change
print(ggplot(as.data.frame(cbind("x"=seq(from=-0.5, to=200.5, by=1), "y"=changed_site_counts/(changed_site_counts+unchanged_site_counts))))  +geom_point(aes(x=x, y=y)) + xlim(c(0,100))) 


dev.off()





# Find general density of nearest sites in GBM segments
# Identify GenomicRegions covering GBM range of such GBM genes (from first to last methylated base in gen 3)
target_ranges = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])))]
target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_ranges))]

# slow method with loop:
z1=cbind(as.data.frame(target_sites), dist_to_nearest = rep(NA,length(target_sites)), type_of_pair = rep(NA,length(target_sites)))

# for each site:
for(i in 1:length(target_sites)) {
  #for(target_site in target_U_sites) {
  #   Find the target_range it is part of
  #target_range = target_ranges[subjectHits(findOverlaps(target_site, target_ranges))]
  target_range = target_ranges[subjectHits(findOverlaps(target_sites[i], target_ranges))]
  #   Find nearest M site in target_range, capture distance in bp and No. sites
  near_sites = setdiff(target_sites[queryHits(findOverlaps(target_sites, target_range))], target_sites[i])
  if (length(near_sites) > 0) {
    # Find the nearest
    #dist_to_Ms = abs(as.data.frame(near_Ms)$start-as.data.frame(target_site)$start)
    dist_to_nearest = abs(as.data.frame(near_sites)$start-as.data.frame(target_sites[i])$start)
    min_dist = dist_to_nearest[dist_to_nearest == min(dist_to_nearest)][1] # [1] takes care of case if there are two equally near
  } else {
    # There are no M sites nearby
    min_dist = NA
  }
  z1[i,"dist_to_nearest"] = min_dist
  #   Find number of times target_site is found M in gen 30 lines
  #z[i, "u_to_m_count"] = length(findOverlaps(target_U_sites[i], u_to_m_call_ranges))

  # What type of pair is this?  MM, MU, UU
    
}

# alternative fast method:
nearest_sites = nearest(target_sites)
dist_to_nearest = distance(target_sites, target_sites[nearest_sites])+1
parental_state = parental_consensus[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
nn_parental_state = parental_consensus[queryHits(findOverlaps(CG_site_ranges, target_ranges))][nearest_sites]
type_of_pair = ifelse(parental_state=="M" & nn_parental_state=="M", "MM", ifelse(parental_state=="U" & nn_parental_state=="U", "UU", ifelse(parental_state=="M" & nn_parental_state=="U", "MU", ifelse(parental_state=="U" & nn_parental_state=="M", "MU", "XX"))))
z1=cbind(as.data.frame(target_sites), dist_to_nearest, type_of_pair)

table(type_of_pair)
#type_of_pair
#    MM     MU     UU     XX 
#165091  63171  66400  76401 


pdf(paste0(project_id,"_dist_to_nearest_CG_site_in_GBM_segments.pdf"))

# plot distribution of gap sizes between adjacent sites
print(ggplot(z1) + geom_histogram(aes(x=dist_to_nearest), binwidth = 1) +facet_wrap(~type_of_pair) + xlim(c(0,200)))
dev.off()

write.table(z1, file="gaps_between_adjacent_GBM_CG_sites.txt", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

# Do the same for UMR segments
target_ranges = UMR_segments.gr
target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_ranges))]

nearest_sites = nearest(target_sites)
dist_to_nearest = distance(target_sites, target_sites[nearest_sites])+1
parental_state = parental_consensus[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
nn_parental_state = parental_consensus[queryHits(findOverlaps(CG_site_ranges, target_ranges))][nearest_sites]
type_of_pair = ifelse(parental_state=="M" & nn_parental_state=="M", "MM", ifelse(parental_state=="U" & nn_parental_state=="U", "UU", ifelse(parental_state=="M" & nn_parental_state=="U", "MU", ifelse(parental_state=="U" & nn_parental_state=="M", "MU", "XX"))))
z1=cbind(as.data.frame(target_sites), dist_to_nearest, type_of_pair)

table(type_of_pair)
#type_of_pair
#     MM      MU      UU      XX 
#     40      64 1498914  240024 

pdf(paste0(project_id,"_dist_to_nearest_CG_site_in_UMR_segments.pdf"))

# plot distribution of gap sizes between adjacent sites
print(ggplot(z1) + geom_histogram(aes(x=dist_to_nearest), binwidth = 1) +facet_wrap(~type_of_pair) + xlim(c(0,200)))
dev.off()

write.table(z1, file="gaps_between_adjacent_UMR_CG_sites.txt", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)


# Find general density of sites, nearest or not, in GBM segments
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
print(ggplot(z2) + geom_histogram(aes(x=dist_to_nearest), binwidth = 1) + xlim(c(0,500)))
dev.off()



# Split GBM-gene average methylation values into quartiles, and replot

meth_quartiles = c(0,quantile(gff.genes[gff.genes$m_class=="Gene-body Methylated","average_methylation"],c(0.25, 0.5, 0.75, type=1)))
meth_quartiles
#                  25%        50%        75%       100% 
#0.00000000 0.08759728 0.18859130 0.32126244 0.77777778 

iteration = 0
for (change_type in c("U_to_M", "M_to_U")) {
  for (nearest in c("M", "U")) {
    for (meth_quartile in 1:4) {
      iteration = iteration + 1
      
      # Make granges for GBM genes with average methylation in the relevant quartile, and find which 'z' sites overlap them
      if (change_type == "U_to_M") {
        z3 = zU[queryHits(findOverlaps(
          makeGRangesFromDataFrame(
            df = zU,
            start.field = "start",
            end.field = "end",
            seqnames.field = "seqnames"
          ),
          makeGRangesFromDataFrame(
            df = gff.genes[(gff.genes$m_class == "Gene-body Methylated") &
                             (gff.genes$average_methylation > meth_quartiles[meth_quartile]) &
                             (gff.genes$average_methylation <= meth_quartiles[meth_quartile + 1]),],
            start.field = "V4",
            end.field = "V5",
            seqnames.field = "V1"
          )
        )),]
        
      } else {
        z3 = zM[queryHits(findOverlaps(
          makeGRangesFromDataFrame(
            df = zM,
            start.field = "start",
            end.field = "end",
            seqnames.field = "seqnames"
          ),
          makeGRangesFromDataFrame(
            df = gff.genes[(gff.genes$m_class == "Gene-body Methylated") &
                             (gff.genes$average_methylation > meth_quartiles[meth_quartile]) &
                             (gff.genes$average_methylation <= meth_quartiles[meth_quartile + 1]),],
            start.field = "V4",
            end.field = "V5",
            seqnames.field = "V1"
          )
        )),]
      }
      
      # Work out histograms of changing and unchanging sites 
      if (change_type == "U_to_M") {
        changed_site_counts = hist(z3[z3$u_to_m_count == 1, paste0("dist_to_nearest_", nearest)], breaks =
                                     seq(
                                       from = 0,
                                       to = 8000,
                                       by = 10
                                     ))$counts
        unchanged_site_counts = hist(z3[z3$u_to_m_count == 0, paste0("dist_to_nearest_", nearest)], breaks =
                                       seq(
                                         from = 0,
                                         to = 8000,
                                         by = 10
                                       ))$counts
        changed_site_counts_detail = hist(z3[z3$u_to_m_count == 1, paste0("dist_to_nearest_", nearest)], breaks =
                                     seq(
                                       from = 0,
                                       to = 8000,
                                       by = 1
                                     ))$counts
        unchanged_site_counts_detail = hist(z3[z3$u_to_m_count == 0, paste0("dist_to_nearest_", nearest)], breaks =
                                       seq(
                                         from = 0,
                                         to = 8000,
                                         by = 1
                                       ))$counts
      } else {
        changed_site_counts = hist(z3[z3$m_to_u_count == 1, paste0("dist_to_nearest_", nearest)], breaks =
                                     seq(
                                       from = 0,
                                       to = 8000,
                                       by = 10
                                     ))$counts
        unchanged_site_counts = hist(z3[z3$m_to_u_count == 0, paste0("dist_to_nearest_", nearest)], breaks =
                                       seq(
                                         from = 0,
                                         to = 8000,
                                         by = 10
                                       ))$counts
        changed_site_counts_detail = hist(z3[z3$m_to_u_count == 1, paste0("dist_to_nearest_", nearest)], breaks =
                                     seq(
                                       from = 0,
                                       to = 8000,
                                       by = 1
                                     ))$counts
        unchanged_site_counts_detail = hist(z3[z3$m_to_u_count == 0, paste0("dist_to_nearest_", nearest)], breaks =
                                       seq(
                                         from = 0,
                                         to = 8000,
                                         by = 1
                                       ))$counts
      }
      
      # Store rates of change by bp bins in accumulation table
      if (iteration == 1) {
        changed_site_props = as.data.frame(
          cbind(
            "dist_metric" = rep("bp", 800),
            "meth_quartile" = rep(meth_quartile, 800),
            "change_type" = rep(change_type, 800),
            "nearest" = rep(nearest, 800),
            "x" = seq(
              from = 5,
              to = 7995,
              by = 10
            ),
            "y" = changed_site_counts / (changed_site_counts + unchanged_site_counts)
          )
        )
        site_detail = as.data.frame(
          cbind(
            "dist_metric" = rep("bp", 8000),
            "meth_quartile" = rep(meth_quartile, 8000),
            "change_type" = rep(change_type, 8000),
            "nearest" = rep(nearest, 8000),
            "x" = seq(
              from = 1,
              to = 8000,
              by = 1
            ),
            "changed" = changed_site_counts_detail,
            "unchanged" = unchanged_site_counts_detail
          )
        )
      } else {
        changed_site_props = rbind(
          changed_site_props,
          cbind(
            "dist_metric" = rep("bp", 800),
            "meth_quartile" = rep(meth_quartile, 800),
            "change_type" = rep(change_type, 800),
            "nearest" = rep(nearest, 800),
            "x" = seq(
              from = 5,
              to = 7995,
              by = 10
            ),
            "y" = changed_site_counts / (changed_site_counts + unchanged_site_counts)
          )
        )
        site_detail = rbind(
          site_detail,
          cbind(
            "dist_metric" = rep("bp", 8000),
            "meth_quartile" = rep(meth_quartile, 8000),
            "change_type" = rep(change_type, 8000),
            "nearest" = rep(nearest, 8000),
            "x" = seq(
              from = 1,
              to = 8000,
              by = 1
            ),
            "changed" = changed_site_counts_detail,
            "unchanged" = unchanged_site_counts_detail
          )
        )
      }
      
      
      # Recalculate changed_site_counts in no. sites bins rather than bp
      if (change_type == "U_to_M") {
        changed_site_counts = hist(z3[z3$u_to_m_count == 1, paste0("sites_to_nearest_", nearest)], breaks =
                                     seq(
                                       from = -0.5,
                                       to = 200.5,
                                       by = 1
                                     ))$counts
        unchanged_site_counts = hist(z3[z3$u_to_m_count == 0, paste0("sites_to_nearest_", nearest)], breaks =
                                       seq(
                                         from = -0.5,
                                         to = 200.5,
                                         by = 1
                                       ))$counts
      } else {
        changed_site_counts = hist(z3[z3$m_to_u_count == 1, paste0("sites_to_nearest_", nearest)], breaks =
                                     seq(
                                       from = -0.5,
                                       to = 200.5,
                                       by = 1
                                     ))$counts
        unchanged_site_counts = hist(z3[z3$m_to_u_count == 0, paste0("sites_to_nearest_", nearest)], breaks =
                                       seq(
                                         from = -0.5,
                                         to = 200.5,
                                         by = 1
                                       ))$counts
      }
      
      # Store rates of change by site in accumulation table
      changed_site_props = rbind(
        changed_site_props,
        cbind(
          "dist_metric" = rep("site", 201),
          "meth_quartile" = rep(meth_quartile, 201),
          "change_type" = rep(change_type, 201),
          "nearest" = rep(nearest, 201),
          "x" = seq(
            from = 0,
            to = 200,
            by = 1
          ),
          "y" = changed_site_counts / (changed_site_counts + unchanged_site_counts)
        )
      )
      
    }
  }
}
changed_site_props$x = as.numeric(levels(changed_site_props$x)[changed_site_props$x])
changed_site_props$y = as.numeric(levels(changed_site_props$y)[changed_site_props$y])
site_detail$x = as.numeric(levels(site_detail$x)[site_detail$x])
site_detail$changed = as.numeric(levels(site_detail$changed)[site_detail$changed])
site_detail$unchanged = as.numeric(levels(site_detail$unchanged)[site_detail$unchanged])

write.table(site_detail, file="changed_and_unchanged_by_meth_quartile.txt", sep="\t")
# Do the plots in no. sites
pdf(paste0(project_id,"_dist_to_nearest_X_all_X_sites_in_GBM_genes_quartiles.pdf"))

# plot rate of change by site
#print(ggplot(changed_site_props[changed_site_props$dist_metric=="site",])  + geom_line(aes(x = x, y = y, colour = meth_quartile))  + xlim(c(0, 100))  + ylim(c(0,0.2)) + xlab("No. CG sites between uCG and nearest mCG at gen3") + ylab("proportion of uCG sites becoming mCG at gen30") + labs(colour = "Gene mean mCG quartile"))
print(ggplot(changed_site_props[changed_site_props$dist_metric=="site",])  + geom_line(aes(x = x, y = y, colour = meth_quartile))  + xlim(c(0, 25))  + ylim(c(0,0.2)) + xlab("No. CG sites between site and nearest X at gen3") + ylab("proportion of sites changing at gen30") + labs(colour = "Gene mean mCG quartile") + facet_grid(rows=vars(nearest), cols=vars(change_type)))

# plot rate of change by bp
#print(ggplot(changed_site_props[changed_site_props$dist_metric=="bp",])  + geom_line(aes(x = x, y = y, colour = meth_quartile))  + xlim(c(0, 2000))  + ylim(c(0,0.75)) + xlab("No. bp between uCG and nearest mCG  at gen3 (10bp bins)") + ylab("proportion of uCG sites becoming mCG at gen30") + labs(colour = "Gene mean mCG quartile"))
print(ggplot(changed_site_props[changed_site_props$dist_metric=="bp",])  + geom_line(aes(x = x, y = y, colour = meth_quartile))  + xlim(c(0, 500))  + ylim(c(0,0.75)) + xlab("No. bp between site and nearest X at gen3 (10bp bins)") + ylab("proportion of sites changing at gen30") + labs(colour = "Gene mean mCG quartile") + facet_grid(rows=vars(nearest), cols=vars(change_type)))

dev.off()


# Iterate through various sizes of masks, plotting the resulting rates of change as it varies with distance to nearest sites
pdf(
  paste0(
    project_id,
    "_dist_to_nearest_X_all_X_sites_in_GBM_genes_quartiles_masked.pdf"
  )
)

# Redo some of the plots, but masking out regions at 5' and 3' of genes
min_segment_length = 100
plot_histograms=TRUE
#for (mask_5_prime in c(100, 200, 500, 1000, 2000)) {
for (mask_5_prime in c(100, 2000)) {
  #for (mask_3_prime in c(100, 200, 500, 1000, 2000)) {
  for (mask_3_prime in c(100, 1000)) {
  #for (mask_3_prime in c(1000)) {
    # Calculate the masked start and end points for all GBM genes
      start_sites = gff.genes$V4 + ifelse(gff.genes$V7 == "+", mask_5_prime, mask_3_prime)
    end_sites = gff.genes$V5 - ifelse(gff.genes$V7 == "+", mask_3_prime, mask_5_prime)
    cat(mask_5_prime, " ", mask_3_prime, " ", nrow(
      cbind(
        "chromosome" = gff.genes$V1,
        "start" = start_sites,
        "end" = end_sites
      )[(end_sites > (start_sites + min_segment_length - 1)) &
          (gff.genes$m_class == "Gene-body Methylated"),]
    ), "\n")
    
    masked_ranges = makeGRangesFromDataFrame(df = as.data.frame(
      cbind(
        "chromosome" = gff.genes$V1,
        "start" = start_sites,
        "end" = end_sites
      )[(end_sites > (start_sites + min_segment_length - 1)) &
          (gff.genes$m_class == "Gene-body Methylated"), ]
    )) #,
    #      start_field = "start",
    #      end_field = "end"
    #    )
    
    
    
    
    # Do the plotting
    iteration = 0
    for (change_type in c("U_to_M", "M_to_U")) {
      for (nearest in c("M", "U")) {
        for (meth_quartile in 1:4) {
          iteration = iteration + 1
          
          # Make granges for masked ranges of GBM genes with average methylation in the relevant quartile, and find which 'z' sites overlap them
          if (change_type == "U_to_M") {
            z3 = zU[queryHits(findOverlaps(makeGRangesFromDataFrame(df = zU), masked_ranges[queryHits(findOverlaps(
              masked_ranges,
              makeGRangesFromDataFrame(
                df = gff.genes[(gff.genes$m_class == "Gene-body Methylated") &
                                 (gff.genes$average_methylation > meth_quartiles[meth_quartile]) &
                                 (gff.genes$average_methylation <= meth_quartiles[meth_quartile + 1]), ],
                start.field = "V4",
                end.field = "V5",
                seqnames.field = "V1"
              )
            )), ])),]
          } else {
            z3 = zM[queryHits(findOverlaps(makeGRangesFromDataFrame(df = zM), masked_ranges[queryHits(findOverlaps(
              masked_ranges,
              makeGRangesFromDataFrame(
                df = gff.genes[(gff.genes$m_class == "Gene-body Methylated") &
                                 (gff.genes$average_methylation > meth_quartiles[meth_quartile]) &
                                 (gff.genes$average_methylation <= meth_quartiles[meth_quartile + 1]), ],
                start.field = "V4",
                end.field = "V5",
                seqnames.field = "V1"
              )
            )), ])),]
          }
          
          # Work out histograms of changing and unchanging sites
          if (change_type == "U_to_M") {
            changed_site_counts = hist(z3[z3$u_to_m_count == 1, paste0("dist_to_nearest_", nearest)], breaks =
                                         seq(
                                           from = 0,
                                           to = 8000,
                                           by = 10
                                         ), plot=plot_histograms)$counts
            unchanged_site_counts = hist(z3[z3$u_to_m_count == 0, paste0("dist_to_nearest_", nearest)], breaks =
                                           seq(
                                             from = 0,
                                             to = 8000,
                                             by = 10
                                           ), plot=plot_histograms)$counts
          } else {
            changed_site_counts = hist(z3[z3$m_to_u_count == 1, paste0("dist_to_nearest_", nearest)], breaks =
                                         seq(
                                           from = 0,
                                           to = 8000,
                                           by = 10
                                         ), plot=plot_histograms)$counts
            unchanged_site_counts = hist(z3[z3$m_to_u_count == 0, paste0("dist_to_nearest_", nearest)], breaks =
                                           seq(
                                             from = 0,
                                             to = 8000,
                                             by = 10
                                           ), plot=plot_histograms)$counts
          }
          
          # Store rates of change by bp bins in accumulation table
          if (iteration == 1) {
            changed_site_props = as.data.frame(
              cbind(
                "dist_metric" = rep("bp", 800),
                "meth_quartile" = rep(meth_quartile, 800),
                "change_type" = rep(change_type, 800),
                "nearest" = rep(nearest, 800),
                "x" = seq(
                  from = 5,
                  to = 7995,
                  by = 10
                ),
                "y" = changed_site_counts / (changed_site_counts + unchanged_site_counts)
              )
            )
          } else {
            changed_site_props = rbind(
              changed_site_props,
              cbind(
                "dist_metric" = rep("bp", 800),
                "meth_quartile" = rep(meth_quartile, 800),
                "change_type" = rep(change_type, 800),
                "nearest" = rep(nearest, 800),
                "x" = seq(
                  from = 5,
                  to = 7995,
                  by = 10
                ),
                "y" = changed_site_counts / (changed_site_counts + unchanged_site_counts)
              )
            )
          }
          
          
          # Recalculate changed_site_counts in no. sites bins rather than bp
          if (change_type == "U_to_M") {
            changed_site_counts = hist(z3[z3$u_to_m_count == 1, paste0("sites_to_nearest_", nearest)], breaks =
                                         seq(
                                           from = -0.5,
                                           to = 200.5,
                                           by = 1
                                         ), plot=plot_histograms)$counts
            unchanged_site_counts = hist(z3[z3$u_to_m_count == 0, paste0("sites_to_nearest_", nearest)], breaks =
                                           seq(
                                             from = -0.5,
                                             to = 200.5,
                                             by = 1
                                           ), plot=plot_histograms)$counts
          } else {
            changed_site_counts = hist(z3[z3$m_to_u_count == 1, paste0("sites_to_nearest_", nearest)], breaks =
                                         seq(
                                           from = -0.5,
                                           to = 200.5,
                                           by = 1
                                         ), plot=plot_histograms)$counts
            unchanged_site_counts = hist(z3[z3$m_to_u_count == 0, paste0("sites_to_nearest_", nearest)], breaks =
                                           seq(
                                             from = -0.5,
                                             to = 200.5,
                                             by = 1
                                           ), plot=plot_histograms)$counts
          }
          
          # Store rates of change by site in accumulation table
          changed_site_props = rbind(
            changed_site_props,
            cbind(
              "dist_metric" = rep("site", 201),
              "meth_quartile" = rep(meth_quartile, 201),
              "change_type" = rep(change_type, 201),
              "nearest" = rep(nearest, 201),
              "x" = seq(
                from = 0,
                to = 200,
                by = 1
              ),
              "y" = changed_site_counts / (changed_site_counts + unchanged_site_counts)
            )
          )
          
        }
      }
    }
    changed_site_props$x = as.numeric(levels(changed_site_props$x)[changed_site_props$x])
    changed_site_props$y = as.numeric(levels(changed_site_props$y)[changed_site_props$y])
    
    # plot rate of change by site
    #print(ggplot(changed_site_props[changed_site_props$dist_metric=="site",])  + geom_line(aes(x = x, y = y, colour = meth_quartile))  + xlim(c(0, 100))  + ylim(c(0,0.2)) + xlab("No. CG sites between uCG and nearest mCG at gen3") + ylab("proportion of uCG sites becoming mCG at gen30") + labs(colour = "Gene mean mCG quartile"))
    print(
      ggplot(changed_site_props[changed_site_props$dist_metric == "site",])  + geom_line(aes(
        x = x, y = y, colour = meth_quartile
      ))  + xlim(c(0, 25))  + ylim(c(0, 0.2)) + xlab("No. CG sites between site and nearest X at gen3") + ylab("proportion of sites changing at gen30") + labs(colour = "Gene mean mCG quartile") + facet_grid(rows =
                                                                                                                                                                                                                  vars(nearest), cols = vars(change_type)) + ggtitle(paste0(
                                                                                                                                                                                                                    "Masks:", mask_5_prime, "bp 5', ", mask_3_prime, "bp 3', N=", nrow(cbind(
                                                                                                                                                                                                                      "chromosome" = gff.genes$V1,
                                                                                                                                                                                                                      "start" = start_sites,
                                                                                                                                                                                                                      "end" = end_sites
                                                                                                                                                                                                                    )[(end_sites > (start_sites + min_segment_length - 1)) &
                                                                                                                                                                                                                        (gff.genes$m_class == "Gene-body Methylated"),]) 
                                                                                                                                                                                                                  ))
    )
    
    # plot rate of change by bp
    #print(ggplot(changed_site_props[changed_site_props$dist_metric=="bp",])  + geom_line(aes(x = x, y = y, colour = meth_quartile))  + xlim(c(0, 2000))  + ylim(c(0,0.75)) + xlab("No. bp between uCG and nearest mCG  at gen3 (10bp bins)") + ylab("proportion of uCG sites becoming mCG at gen30") + labs(colour = "Gene mean mCG quartile"))
    print(
      ggplot(changed_site_props[changed_site_props$dist_metric == "bp",])  + geom_line(aes(
        x = x, y = y, colour = meth_quartile
      ))  + xlim(c(0, 500))  + ylim(c(0, 0.75)) + xlab("No. bp between site and nearest X at gen3 (10bp bins)") + ylab("proportion of sites changing at gen30") + labs(colour = "Gene mean mCG quartile") + facet_grid(rows =
                                                                                                                                                                                                                          vars(nearest), cols = vars(change_type)) + ggtitle(
                                                                                                                                                                                                                            paste0("Masks:", mask_5_prime, "bp 5', ", mask_3_prime, "bp 3', N=", nrow(cbind(
                                                                                                                                                                                                                              "chromosome" = gff.genes$V1,
                                                                                                                                                                                                                              "start" = start_sites,
                                                                                                                                                                                                                              "end" = end_sites
                                                                                                                                                                                                                            )[(end_sites > (start_sites + min_segment_length - 1)) &
                                                                                                                                                                                                                                (gff.genes$m_class == "Gene-body Methylated"),]))
                                                                                                                                                                                                                          )
    )
  }
}
dev.off()

# There is a 3bp periodicity in distance to nearest M data from both U->U sites and U->M sites.
# Rates of change of U->M peak at 2, 7, 11, 18, and 167 bp

# This raises the hypothesis that sites which change U->M which are distant from the nearest M site by 160-170 bp are enriched in sites occurring in linker DNA, in comparison with, say, sites at other distances.

# Can we plot distance from nearest dyad as a function of distance to nearest M site, among sites that change?

### Comparison with nucleosome positions
### Read in the nucleosome positions as a dump of genomicRanges which was previously prepared in a VM (Windows cannot import BigWig files)
nucleosome_positions.gr = makeGRangesFromDataFrame(df = read.table(reference_nucleosomes), start.field = "V2", end.field = "V3", seqnames.field = "V1")

levels(nucleosome_positions.gr@seqnames)=toupper(levels(nucleosome_positions.gr@seqnames))
nucleosome_positions.gr@seqinfo@seqnames=levels(nucleosome_positions.gr@seqnames)

# Identify central site in each nucleosome (location of dyad)
dyad_positions.gr = makeGRangesFromDataFrame(df = cbind("start"=nucleosome_positions.gr@ranges@start+(nucleosome_positions.gr@ranges@width-1)/2, "chrom"=as.data.frame(nucleosome_positions.gr)$seqnames), start.field = "start", end.field = "start", seqnames.field = "chrom")


#target_ranges = GBM_segments.gr[queryHits(findOverlaps(GBM_segments.gr, setdiff(gene_body_loci.gr,gene_body_loci.gr[queryHits(findOverlaps(gene_body_loci.gr, TEM_segments.gr))])))]
#target_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, target_ranges))]
#target_sites = makeGRangesFromDataFrame(df = zU[zU$u_to_m_count==1,], start.field = "start", end.field = "end", seqnames.field = "seqnames")
target_sites = makeGRangesFromDataFrame(df = zU, start.field = "start", end.field = "end", seqnames.field = "seqnames")
z1=cbind(as.data.frame(target_sites), dist_to_nearest = rep(NA,length(target_sites)))

	
# for each site:
window_width = 500
for(i in 1:length(target_sites)) {
  cat(paste(i,":\n"))
  #target_site = makeGRangesFromDataFrame(df = zU[i], start.field = "start", end.field = "end", seqnames.field = "seqnames")
  target_range = makeGRangesFromDataFrame(df = cbind("seqnames" = zU[i,c("seqnames")],"start" = zU[i,c("start")]-window_width/2, "end" = zU[i,c("start")]+window_width/2), start.field = "start", end.field = "end", seqnames.field = "seqnames")

  near_sites = dyad_positions[queryHits(findOverlaps(dyad_positions, target_range))]
  if (length(near_sites) > 0) {
    dist_to_nearest = abs(as.data.frame(near_sites)$start-as.data.frame(target_sites[i])$start)
    min_dist = dist_to_nearest[dist_to_nearest == min(dist_to_nearest)][1] # [1] takes care of case if there are two equally near
  } else {
    # There are no dyads nearby
    min_dist = NA
  }
  z1[i,"dist_to_nearest"] = min_dist
  cat(paste(length(near_sites),min_dist,":\n"))
}


	
	
	
	
	# Plot nuclosome coverage +- 1kb relative to variable sites
	
	# Remake the boundary plots by calculating the number of times a relevant window overlaps with a nucleosome
	# This time, when we calculate rates we will use No. nucleosome overlaps / No. hemi-segments in boundary

	window_size = 1 # window length in nt
	step_size=window_size # No. nt to slide window per step
	window_range = 250 # No steps to take each side of central site
	
	window_site_summary=data.frame()
	first_window = TRUE

	no_width_quantiles = 3
	for (width_quantile in 1:no_width_quantiles) {
		for (focus_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
			min_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile-1)/no_width_quantiles)
			max_width = quantile(as.data.frame(segmentation_model.gr)[as.data.frame(segmentation_model.gr)$segment.mStatus==focus_mStatus,]$width,(width_quantile)/no_width_quantiles)
			
			for (adjacent_mStatus in c("GBM", "GBM-like", "TEM", "UMR")) {
	
				if (focus_mStatus != adjacent_mStatus) {
					#focus_mStatus = "UMR"
					#adjacent_mStatus = "TEM"
	
					# Create set of central sites to be the segments on the right hand side of the requested boundaries
					central_sites = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@second)
				
					# Create corresponding set of left hand segments for each boundary
					left_hand_segments = as.data.frame(subset(findOverlapPairs(segmentation_model.gr[segmentation_model.gr$segment.mStatus==focus_mStatus], segmentation_model.gr[(segmentation_model.gr$segment.mStatus==adjacent_mStatus) & (segmentation_model.gr@ranges@width>min_width) & (segmentation_model.gr@ranges@width<=max_width)], maxgap=1L), end(first) == start(second) - 1L)@first)

					# We want to avoid accumulating data from outside of the segments under consideration.  Accordingly, create genomicranges corresponding to the adjoining halfs of each of the two adjoining segments.  We only use adjoining halfs of the segments to avoid edge effects from the other end of the segments
					relevant_sites_end = central_sites$end-round((central_sites$end-central_sites$start)/2)
					relevant_sites_start = left_hand_segments$end-round((left_hand_segments$end-left_hand_segments$start)/2)

					cat(paste0(nrow(central_sites)," boundaries between ",focus_mStatus," and ",adjacent_mStatus," segments.\n"))
					#	levels(central_sites@seqnames)=toupper(levels(central_sites@seqnames))
					#	central_sites@seqinfo@seqnames=levels(central_sites@seqnames)


					for (window_no in -window_range:window_range) {
						# for upstream and downstream:
						# make a genomic ranges for the window
						# calculate the No. times the window is present in a relevant segment boundary, and No. times the window overlaps a nucleosome
						# store the resulting values 
		
						window_sites=data.frame(cbind("seqnames"=as.character(central_sites$seqnames), "start"=central_sites$start + (window_no*step_size) - round((window_size/2) - 1.25), "end"=central_sites$start + (window_no*step_size) +round((window_size/2) + 0.25), "width"=window_size, "strand"=as.character(central_sites$strand)))

						window_sites$start=as.numeric(levels(window_sites$start)[window_sites$start])
						window_sites$end=as.numeric(levels(window_sites$end)[window_sites$end])
		
						window_sites$start = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$start<relevant_sites_start, relevant_sites_start, window_sites$start)))
						window_sites$end = ifelse(window_sites$end<=relevant_sites_start,NA,ifelse(window_sites$start>=relevant_sites_end,NA,ifelse(window_sites$end>relevant_sites_end, relevant_sites_end, window_sites$end)))

						# Remove any windows which did not overlap relvant sites at all
						window_sites = window_sites[!is.na(window_sites$start),]
					
						# Remove any windows which overlap ends of chromosomes
						if(window_no<0) {
							window_sites = window_sites[window_sites$end>=1,]
							window_sites = window_sites[window_sites$start<=sLengths[window_sites$seqnames],]
						} else {
							window_sites = window_sites[window_sites$start>=1,]
							window_sites = window_sites[window_sites$end<=sLengths[window_sites$seqnames],]
						}
		
						window_sites.gr = makeGRangesFromDataFrame(df = window_sites, start.field = "start", end.field = "end", seqnames.field = "seqnames")

						window_sites.gr$ID=row.names(as.data.frame(window_sites.gr))
						#window_sites.gr$CG_site_count=table(window_sites.gr[subjectHits(findOverlaps(CG_site_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$CG_site_count[is.na(window_sites.gr$CG_site_count)]=0
						#window_sites.gr$variant_count=table(window_sites.gr[subjectHits(findOverlaps(variant_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$variant_count[is.na(window_sites.gr$variant_count)]=0
						window_sites.gr$nucleosome_count=table(window_sites.gr[subjectHits(findOverlaps(nucleosome_positions.gr, window_sites.gr)),]$ID)[window_sites.gr$ID]
						window_sites.gr$nucleosome_count[is.na(window_sites.gr$nucleosome_count)]=0
						#window_sites.gr$all_M_count=table(window_sites.gr[subjectHits(findOverlaps(all_M_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$all_M_count[is.na(window_sites.gr$all_M_count)]=0
						#window_sites.gr$all_U_count=table(window_sites.gr[subjectHits(findOverlaps(all_U_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$all_U_count[is.na(window_sites.gr$all_U_count)]=0
						#window_sites.gr$m_to_u_count=table(window_sites.gr[subjectHits(findOverlaps(m_to_u_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$m_to_u_count[is.na(window_sites.gr$m_to_u_count)]=0
						#window_sites.gr$u_to_m_count=table(window_sites.gr[subjectHits(findOverlaps(u_to_m_call_ranges, window_sites.gr)),]$ID)[window_sites.gr$ID]
						#window_sites.gr$u_to_m_count[is.na(window_sites.gr$u_to_m_count)]=0
		
						# ensure all window sites overlap with one or more CHG sites
						#window_sites.gr = window_sites.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, window_sites.gr)))]

						#window_meths=meth_by_segment(meth_parents_CHG.gr, segment_model=window_sites.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(0.2), mMeth.classes=c("UMR","MR"))
		
						cat(paste0("Window: ",window_no," nucleosome proportion: ",mean(window_sites.gr$nucleosome_count, na.rm=TRUE), "variance: ",var(window_sites.gr$nucleosome_count, na.rm=TRUE),"\n"))
			
						if (first_window) {
							#window_site_summary=rbind.data.frame(window_site_summary, c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$nucleosome_count, na.rm=TRUE), var(window_sites.gr$nucleosome_count, na.rm=TRUE), mean(window_sites.gr$m_to_u_count, na.rm=TRUE), var(window_sites.gr$m_to_u_count, na.rm=TRUE), mean(window_sites.gr$u_to_m_count, na.rm=TRUE), var(window_sites.gr$u_to_m_count, na.rm=TRUE), mean(window_sites.gr$CG_site_count, na.rm=TRUE), var(window_sites.gr$CG_site_count, na.rm=TRUE), mean(window_sites.gr$all_M_count, na.rm=TRUE), var(window_sites.gr$all_M_count, na.rm=TRUE), mean(window_sites.gr$all_U_count, na.rm=TRUE), var(window_sites.gr$all_U_count, na.rm=TRUE), nrow(window_sites), width_quantile))
							#names(window_site_summary) = c("category","position","mean_vmCG","var_vmCG","mean_vMU","var_vMU","mean_vUM","var_vUM","mean_nCG","var_nCG","mean_all_M","var_all_M","mean_all_U","var_all_U","seg_count","width_quantile")
							window_site_summary=rbind.data.frame(window_site_summary, c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$nucleosome_count, na.rm=TRUE), var(window_sites.gr$nucleosome_count, na.rm=TRUE), sum(window_sites.gr$nucleosome_count, na.rm=TRUE), nrow(window_sites), width_quantile))
							names(window_site_summary) = c("category","position","mean_nucl","var_nucl","nucl_count","seg_count","width_quantile")
							window_site_summary$category = as.character(window_site_summary$category)
							window_site_summary$position = as.integer(levels(window_site_summary$position))[window_site_summary$position]
							window_site_summary$mean_nucl = as.numeric(levels(window_site_summary$mean_nucl))[window_site_summary$mean_nucl]
							window_site_summary$var_nucl = as.numeric(levels(window_site_summary$var_nucl))[window_site_summary$var_nucl]
							#window_site_summary$mean_vmCG = as.numeric(levels(window_site_summary$mean_vmCG))[window_site_summary$mean_vmCG]
							#window_site_summary$var_vmCG = as.numeric(levels(window_site_summary$var_vmCG))[window_site_summary$var_vmCG]
							#window_site_summary$mean_vMU = as.numeric(levels(window_site_summary$mean_vMU))[window_site_summary$mean_vMU]
							#window_site_summary$var_vMU = as.numeric(levels(window_site_summary$var_vMU))[window_site_summary$var_vMU]
							#window_site_summary$mean_vUM = as.numeric(levels(window_site_summary$mean_vUM))[window_site_summary$mean_vUM]
							#window_site_summary$var_vUM = as.numeric(levels(window_site_summary$var_vUM))[window_site_summary$var_vUM]
							#window_site_summary$mean_nCG = as.numeric(levels(window_site_summary$mean_nCG))[window_site_summary$mean_nCG]
							#window_site_summary$var_nCG = as.numeric(levels(window_site_summary$var_nCG))[window_site_summary$var_nCG]
							#window_site_summary$mean_all_M = as.numeric(levels(window_site_summary$mean_all_M))[window_site_summary$mean_all_M]
							#window_site_summary$var_all_M = as.numeric(levels(window_site_summary$var_all_M))[window_site_summary$var_all_M]
							#window_site_summary$mean_all_U = as.numeric(levels(window_site_summary$mean_all_U))[window_site_summary$mean_all_U]
							#window_site_summary$var_all_U = as.numeric(levels(window_site_summary$var_all_U))[window_site_summary$var_all_U]
							window_site_summary$nucl_count = as.integer(levels(window_site_summary$nucl_count))[window_site_summary$nucl_count]
							window_site_summary$seg_count = as.integer(levels(window_site_summary$seg_count))[window_site_summary$seg_count]
							window_site_summary$width_quantile = as.integer(levels(window_site_summary$width_quantile))[window_site_summary$width_quantile]
							first_window = FALSE
						} else {
							window_site_summary=rbind.data.frame(window_site_summary, setNames(c(paste0("Variable ",focus_mStatus,"-",adjacent_mStatus," CG sites"), window_no*step_size, mean(window_sites.gr$nucleosome_count, na.rm=TRUE), var(window_sites.gr$nucleosome_count, na.rm=TRUE), sum(window_sites.gr$nucleosome_count, na.rm=TRUE), nrow(window_sites), width_quantile), names(window_site_summary)))
						}
				#		window_summary=rbind.data.frame(window_summary, setNames(c("Random UMR CG sites", window_no*step_size, mean(window_sites.gr$variant_count, na.rm=TRUE), var(window_sites.gr$variant_count, na.rm=TRUE)), names(window_summary)))
				#		colnames(window_summary) = c("category","position","mean_vmCG","var_vmCG")

					}
	
					#window_summary$category = as.character(window_summary$category)
					window_site_summary$position = as.integer(window_site_summary$position)
					window_site_summary$mean_nucl = as.numeric(window_site_summary$mean_nucl)
					window_site_summary$var_nucl = as.numeric(window_site_summary$var_nucl)
					#window_site_summary$mean_vmCG = as.numeric(window_site_summary$mean_vmCG)
					#window_site_summary$var_vmCG = as.numeric(window_site_summary$var_vmCG)
					#window_site_summary$mean_vMU = as.numeric(window_site_summary$mean_vMU)
					#window_site_summary$var_vMU = as.numeric(window_site_summary$var_vMU)
					#window_site_summary$mean_vUM = as.numeric(window_site_summary$mean_vUM)
					#window_site_summary$var_vUM = as.numeric(window_site_summary$var_vUM)
					#window_site_summary$mean_nCG = as.numeric(window_site_summary$mean_nCG)
					#window_site_summary$var_nCG = as.numeric(window_site_summary$var_nCG)
					#window_site_summary$mean_all_M = as.numeric(window_site_summary$mean_all_M)
					#window_site_summary$var_all_M = as.numeric(window_site_summary$var_all_M)
					#window_site_summary$mean_all_U = as.numeric(window_site_summary$mean_all_U)
					#window_site_summary$var_all_U = as.numeric(window_site_summary$var_all_U)
					window_site_summary$nucl_count = as.integer(window_site_summary$nucl_count)
					window_site_summary$seg_count = as.integer(window_site_summary$seg_count)
					window_site_summary$width_quantile = as.integer(window_site_summary$width_quantile)

					# Plot mean and variance of vmCG in windows centred at segment boundaries
					#pdf(paste0(project_id,"_",meth_context,"_",focus_mStatus,"_",adjacent_mStatus,"_variable_CG_sites_windows_",CHG_window_size,".pdf"))
					#print(ggplot() +geom_pointrange(data=window_site_summary, aes(x=position, y=mean_vmCG, ymin=mean_vmCG-var_vmCG/2, ymax=mean_vmCG+var_vmCG/2, colour=category)))
					#dev.off()
				}
			}
		}	
	} # for each quantile
	
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−TEM CG sites", "Variable TEM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
	#ggplot() +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vmCG, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vMU, colour=category)) +geom_line(data=window_site_summary[window_site_summary$category %in% c("Variable UMR−GBM CG sites", "Variable GBM-UMR CG sites"),], aes(x=position, y=mean_vUM, colour=category))
		
		

	# Reshape data into one row per quantity plotted per window (count data)
	paneltext=array()
	paneltext[1]=paste0("Proportion of segments overlapping nucleosomes per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Feature"=rep("Nuclosomes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$nucl_count, "site_count"=window_site_summary$seg_count, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Feature"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$feature_count=as.numeric(levels(window_site_detail$feature_count)[window_site_detail$feature_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
		
	# Plot counts of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_count_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}	
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_count_by_segment_boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()
	
	
	
	# Reshape data into one row per quantity plotted per window (rate data - nucleosomes per segment)
	window_site_detail = rbind.data.frame(cbind("Feature"=rep("Nuclosomes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$nucl_count/window_site_summary$seg_count, "site_count"=window_site_summary$seg_count, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile), cbind("Feature"=rep("Segments", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$seg_count, "site_count"=rep(1, nrow(window_site_summary)), "panel"=rep(paneltext[2],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$feature_count=as.numeric(levels(window_site_detail$feature_count)[window_site_detail$feature_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])
		
	# Plot counts of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop_by_segment_boundary_by_width_windows_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}	
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()

	# Plot counts of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop_by_segment_boundary_by_width_windows_mirror_",window_size,".pdf"))
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	for (width_quantile in 1:no_width_quantiles) {
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")) & (window_site_detail$width_quantile==width_quantile),], aes(x=position, y=feature_count, colour=Feature)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nSegment length quantile ",width_quantile)) + facet_grid(panel~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	}
	dev.off()
	
	

	# Reshape data into one row per quantity plotted per window (rate data - nucleosomes per segment), but without seg_count data
	paneltext=array()
	paneltext[1]=paste0("Proportion of segments overlapping nucleosomes per ",window_size,"nt sliding window")
	paneltext[2]="Number of segments containing window"

	window_site_detail = rbind.data.frame(cbind("Feature"=rep("Nuclosomes", nrow(window_site_summary)), "position"=window_site_summary$position, "category"=window_site_summary$category, "feature_count"=window_site_summary$nucl_count, "site_count"=window_site_summary$seg_count, "panel"=rep(paneltext[1],nrow(window_site_summary)), "width_quantile"=window_site_summary$width_quantile))

	window_site_detail$position=as.integer(window_site_summary$position)
	window_site_detail$feature_count=as.numeric(levels(window_site_detail$feature_count)[window_site_detail$feature_count])
	window_site_detail$site_count=as.numeric(levels(window_site_detail$site_count)[window_site_detail$site_count])

	# Plot rates of transitions in each direction
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop4_by_segment_boundary_by_width_windows_",window_size,".pdf"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM Segments\nby GBM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−TEM CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | TEM Segments\nby TEM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable UMR−GBM-like CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("UMR Segments | GBM-like Segments\nby GBM-like segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()

	# Plot rates of transitions in each direction at opposite segment boundaries
	# We anticipate these will be mirror images of the above plots
	pdf(paste0(project_id,"_",meth_context,"_nucleosome_prop4_by_segment boundary_by_width_windows_mirror_",window_size,".pdf"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("GBM Segments | UMR Segments\nby GBM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable TEM-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("TEM Segments | UMR Segments\nby TEM segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
		print(ggplot() +geom_line(data=window_site_detail[(window_site_detail$category %in% c("Variable GBM-like-UMR CG sites")),], aes(x=position, y=ifelse(site_count==0,0,feature_count/site_count), colour=Feature)) + theme_minimal() + xlab(paste0("GBM-like Segments | UMR Segments\nby GBM-like segment length quantile")) + facet_grid(width_quantile~., scale="free", switch="y") + theme( axis.title.y = element_blank(), strip.placement="outside"))
	dev.off()
