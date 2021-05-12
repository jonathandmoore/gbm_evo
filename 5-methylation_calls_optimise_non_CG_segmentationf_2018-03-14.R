
	
	
# Known bugs:

# Genes with 'Unknown' methylation status (too high or low coverage) are masked out.  Some of these are microRNAs, and probably should not be masked out.  In fact, should we really be masking out any of these genes?  (Fixed by removing 'Unknown' genes from the mask.)
# Could consider, instead, masking out high- and low-coverage zones based on a segmentation of coverage.

# We should use both CHH and CHG thresholds (or in combination) to make a determination between GBM and TEM segments.  Currently, there are a number of short segments with low mCHG (and few CHG sites) and high mCHH.

# There are some cases where GBM and TEM segments adjoin, and where, really, both should be TEM - get a list and go through them.

# Where gaps occur due to short segments with either no CG, or no CHG sites, first, check if the segment overlaps any TEs and exons.  If no and yes, respectively, then classify it as GBM 

	
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
	#meth_context="CG"
	#meth_context="CHG"
	meth_context="CHH"
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

install.packages("ggpubr")

	
#### This part sets up libraries and contextual information and metadata about the project

library(data.table)
library(stringr)
library(ggplot2)
library(plyr)
library(ggpubr)

# Platform-specific stuff:
#  Base of path to project
#  Where to find bioconductor (unresolved on cluster at the moment)
pathroot = ""
if (.Platform$OS.type=="windows") {
	pathroot="D:"
	source("http://bioconductor.org/biocLite.R")
} else if (.Platform$OS.type=="unix") {
	if(substr(getwd(),1,20)=="/nbi/Research-Groups") {
		# assume we are running on cluster
		pathroot="/nbi/Research-Groups/JIC/Daniel-Zilberman"
	} else {
		# assume we are running in Virtualbox with shared folder /media/sf_D_DRIVE
		pathroot="/media/sf_D_DRIVE"
		source("http://bioconductor.org/biocLite.R")
	}
}

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

# For nucleosome positioning data need to first read in from BigWig file, but Windows can't do this so did it using VM, then exported resulting genomicranges object.
# This is the code to run in R in the VM:
	#source("http://www.bioconductor.org/biocLite.R")
	#biocLite("rtracklayer")
	#library(rtracklayer)
	#saveRDS(import("/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.bw", format="BigWig"), "/media/sf_D_DRIVE/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.RDS")
reference_nucleosomes = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Nucleosomes/Ath_leaf_NPS.RDS")  # Nucleosomes position data from http://plantdhs.org/Download (Jiming Jiang lab) 


# get a list of the samples from the subdirectories of the source directory
samples = list.dirs(path = source_dir, full.names = FALSE, recursive = FALSE)
no_samples=length(samples)

# Read the sample metadata in from a text file in the source directory
sample_metadata=read.delim(paste0(source_dir,project_id,"_metadata.tsv"), header=TRUE, sep="\t")
rownames(sample_metadata)=sample_metadata$Identifier

sample_info <- readRDS(paste0(project_id,"_",meth_context,"_sample_info.rds"))

	sample_info=cbind(sample_info,sample_metadata)

	# Load up genomic ranges and get what we need from the annotation
	biocLite("GenomicRanges")
	library(GenomicRanges)
	
	gff.genes=readRDS(paste0(project_id,"_CG_gff.genes.rds"))
	gff.transposons=readRDS(paste0(project_id,"_CG_gff.transposons.rds"))
	
	# Load the BSGenome version of the reference 
	biocLite("BSgenome")
	library(BSgenome.Athaliana.TAIR.TAIR9)

	# Find the sequence lengths for each chromosome
	sLengths=seqlengths(Athaliana)


	# Use the MethylSeekR package to segment the methylome
	library(MethylSeekR)

	# Define some loci to mask out from the genome (mitochondria, chloroplast, genes found previously to have coverage anomalies)
	mask_loci=data.frame(Chromosome=character(), start_site=integer(), end_site=integer())
	mask_loci=rbind(mask_loci,list(Chromosome="ChrM",start_site=1,end_site=sLengths[6]),list(Chromosome="ChrC",start_site=1,end_site=sLengths[7]))
	#colnames(mask_loci)=c("Chromosome","start_site","end_site")

	# Commented out the part which masks out genes annotated as mCG 'Unknown'
	#unknown_gene_loci=gff.genes[gff.genes$m_class=="Unknown",c("V1","V4","V5")]
	#colnames(unknown_gene_loci)=colnames(mask_loci)
	#unknown_gene_loci$Chromosome=paste0(substr(unknown_gene_loci$Chromosome,1,1),tolower(substr(unknown_gene_loci$Chromosome,2,3)),substr(unknown_gene_loci$Chromosome,4,4))
	#mask_loci=rbind.data.frame(mask_loci,unknown_gene_loci)
	#rm(unknown_gene_loci)

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

	# Set chromosome names to mixed case in genome segments object
	levels(genome_loci.gr@seqnames@values) = mixedCaseChr(as.character(levels(genome_loci.gr@seqnames@values)))
	genome_loci.gr@seqinfo@seqnames = levels(genome_loci.gr@seqnames@values)
	levels(mask_loci.gr@seqnames@values) = mixedCaseChr(as.character(levels(mask_loci.gr@seqnames@values)))
	mask_loci.gr@seqinfo@seqnames = levels(mask_loci.gr@seqnames@values)

	methylated_loci.gr=reduce(makeGRangesFromDataFrame(df=rbind(gff.transposons[,c("V1","V4","V5")],gff.genes[gff.genes$m_class=="Heterochromatic",c("V1","V4","V5")]), start.field="V4", end.field="V5", seqnames.field="V1"))
	gene_body_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Gene-body Methylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))
	unmethylated_loci.gr=reduce(makeGRangesFromDataFrame(df=gff.genes[gff.genes$m_class=="Unmethylated",c("V1","V4","V5")], start.field="V4", end.field="V5", seqnames.field="V1"))

	# Concatenate the CHG and CHH methylomes, and sort them, by position and chromosome
	meth_parents_CHG_CHH.gr <- sort(sortSeqlevels(c(readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_parents.tsv"), seqLengths=sLengths), readMethylome(FileName=paste0(project_id,"_CHH_TM_read_counts_across_parents.tsv"), seqLengths=sLengths))))


	m_n_non_CG_optimisation=data.frame(m=numeric(), n=numeric(), UMR=integer(), LMR=integer(), fp=numeric(), tp=numeric())

	# Loop through values for m and n parameters 
	#m_CHG.seq=seq(0.002,0.1, by=0.002)
	m_non_CG.seq=seq(0.05,0.5, by=0.05)
	#n_CHG.seq=seq(1,20, by=1)
	n_non_CG.seq=seq(50,500, by=50)
	for (m_non_CG.sel in m_non_CG.seq) {
		for (n_non_CG.sel in n_non_CG.seq) {
			#n.sel=4
			#UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHG.gr, meth.cutoff=m_CHG.sel, nCpG.cutoff=n_CHG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
			#UMRLMRsegments_CHG.gr <- segmentUMRsLMRs(m=meth_parents_CHH.gr, meth.cutoff=m_CHG.sel, nCpG.cutoff=n_CHG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
			UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_parents_CHG_CHH.gr, meth.cutoff=m_non_CG.sel, nCpG.cutoff=n_non_CG.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

			# Capitalise chromosome names in segmentation object
			levels(UMRLMRsegments_non_CG.gr@seqnames)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
			UMRLMRsegments_non_CG.gr@seqinfo@seqnames=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")

			# Find overlaps between segments and positive and negative 'control' loci (annotated genes and transposons)
			fp_hits=findOverlaps(methylated_loci.gr, UMRLMRsegments_non_CG.gr)
			fp_overlaps <- pintersect(methylated_loci.gr[queryHits(fp_hits)], UMRLMRsegments_non_CG.gr[subjectHits(fp_hits)])
			fp_overlap_prop <- sum(width(fp_overlaps)) / sum(width(methylated_loci.gr))
			tp_hits=findOverlaps(unmethylated_loci.gr, UMRLMRsegments_non_CG.gr)
			tp_overlaps <- pintersect(unmethylated_loci.gr[queryHits(tp_hits)], UMRLMRsegments_non_CG.gr[subjectHits(tp_hits)])
			tp_overlap_prop <- sum(width(tp_overlaps)) / sum(width(unmethylated_loci.gr))
			
			# Find overlaps between segments and variable loci
#			variable_sites_UMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments_CHG.gr[UMRLMRsegments_CHG.gr@elementMetadata$type=="UMR"])
#			variable_sites_LMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments_CHG.gr[UMRLMRsegments_CHG.gr@elementMetadata$type=="LMR"])
			m_n_non_CG_optimisation=rbind(m_n_non_CG_optimisation, list(m=m_non_CG.sel, n=n_non_CG.sel, UMR=0, LMR=0, fp=fp_overlap_prop, tp=tp_overlap_prop))
			# How many variable sites were captured?
#			cat(paste0("m.sel=",m_CHG.sel," n.sel=",n_CHG.sel," Variable sites in UMRs: ",length(variable_sites_UMR_olaps@from),"  LMRs: ",length(variable_sites_LMR_olaps@from),"\n"))
		}
	}
	m_n_non_CG_optimisation$variable_sites_captured=m_n_non_CG_optimisation$UMR+m_n_non_CG_optimisation$LMR

	# Plot ROC curves for each value of n, and estimate AUC
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_m_n_non_CG_optimisation_parents_ROC.pdf"))
	print(ggplot(m_n_non_CG_optimisation, aes(x=fp, y=tp)) + geom_point() +geom_text(aes(label=paste0(m)),hjust=0, vjust=0) +geom_line() +facet_wrap(~ n))
	dev.off()
	
	# Estimate AUC from trapezium approximation
	cat(paste0("Estimated AUC values from ROC curves generated by varying m parameter for each value of n in MethylSeekR:\n"))
	best_n_non_CG = 0
	prev_best_auc_non_CG = 0
	for (n_non_CG.sel in n_non_CG.seq) {
		m_n_non_CG_optimisation_auc=0
		prev_tp=0
		prev_fp=0
		for (m_non_CG.sel in m_non_CG.seq) {
			m_n_non_CG_optimisation_auc = m_n_non_CG_optimisation_auc + (1-m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$fp)*(m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$tp-prev_tp) + ((m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$tp-prev_tp)*(m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$fp-prev_fp))/2
			prev_tp=m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$tp
			prev_fp=m_n_non_CG_optimisation[(m_n_non_CG_optimisation$m==m_non_CG.sel) & (m_n_non_CG_optimisation$n==n_non_CG.sel),]$fp
		}
		if (m_n_non_CG_optimisation_auc > prev_best_auc_non_CG) {
			best_n_non_CG = n_non_CG.sel
			prev_best_auc_non_CG = m_n_non_CG_optimisation_auc
		}
		cat(paste0(" n=",n_non_CG.sel," AUC=",m_n_non_CG_optimisation_auc,"\n"))
	}
	cat(paste0("n=",best_n_non_CG," maximises AUC (",prev_best_auc_non_CG,").\n"))
	
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
 
	# Free up some memory by removing largest object(s) and garbage collection
	### This is not always effective (R may not defrag memory effectively).  If that happens, and you still run out of memory when trying to segment after gc(), try:
	### Quit R; rerun everything before the above segmentation m/n optimisation loop; set best_n_non_CG by hand to the value found previously; then continue here:
	rm(UMRLMRsegments_non_CG.gr)
	gc()
 
	# Segment using best_n and m=0.9 to find appropriate cutoff for m.sel
	UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_parents_CHG_CHH.gr, meth.cutoff=0.9, nCpG.cutoff=best_n_non_CG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_non_CG_segmentation_landscape_",best_n_non_CG,"_0_9.pdf"))

	saveUMRLMRSegments(segs=UMRLMRsegments_non_CG.gr, GRangesFilename=paste0("UMRsLMRs_parents_non_CG_0.9_",best_n_non_CG,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_non_CG_0.9_",best_n_non_CG,".tsv"))
	# UMRLMRsegments_non_CG.gr <- readRDS(paste0("UMRsLMRs_parents_non_CG_0.9_",best_n_non_CG,".gr.rds"))
	
	
	# More formally:
	# Plot density distribution of segment mCH*
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_non_CG_segmentation_landscape_",best_n_non_CG,"_0_9_density.pdf"))
	print(ggplot(as.data.frame(UMRLMRsegments_non_CG.gr), aes(x=pmeth)) + geom_histogram(binwidth=0.01))
	dev.off()

	# Visual inspection of SRA035939_CHH_MethylSeekR_non_CG_segmentation_landscape_500_0_9.pdf shows three clear populations of segments.  Short methylated segments (n<2^5, m>0.15), and longer unmethylated and low-methylated segments (n>2^10, m<0.1; m>0.1).

	# Fit mixture of three Gaussians to mCHG density to identify cutoff
	
	# Initially tried fitting 3 Gaussians to whole data set, but third Gaussian was used to cover the scattering of high-pMeth sgments (m>0.3).  Accordingly, fit to segments with m<0.3 only
	
	#install.packages("MASS")  # for fitting negative binomial models to data
	#library(MASS)
	install.packages("mixtools")
	library(mixtools)
	m_non_CG_mixmdl = normalmixEM(as.data.frame(UMRLMRsegments_non_CG.gr)[(!is.na(as.data.frame(UMRLMRsegments_non_CG.gr)$pmeth)) & (as.data.frame(UMRLMRsegments_non_CG.gr)$pmeth<0.3),]$pmeth, k=3)
	# Had to exclude sites with NA CHG methylation or the fitting algorithm exploded	
	
	# Print the characteristics of the fit curves
	cat(paste0("Fit lambda values (mixture proportions):\n"))
	m_non_CG_mixmdl$lambda
	# Schmitz data:   
	# Becker data: 
	
	cat(paste0("Fit mu values (means):\n"))
	m_non_CG_mixmdl$mu
	# Schmitz data:
	# Becker data:  
	
	cat(paste0("Fit sigma values (st.devs):\n"))
	m_non_CG_mixmdl$sigma
	# Schmitz data: 
	# Becker data:  
	
	# Find the intersection of the two larger components in the mixture - here we cut off

	#packages were installed and mixfn defined earlier
	#install.packages("rootSolve")
	library(rootSolve)
	mixfn <- function(x, m1, sd1, m2, sd2, p1, p2){
	dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 }

	m_non_CG_low_model=1
	m_non_CG_high_model=m_non_CG_low_model+1
	segment_m_non_CG_cutoff=uniroot.all(mixfn, lower=0, upper=1, m1=m_non_CG_mixmdl$mu[m_non_CG_low_model], sd1=m_non_CG_mixmdl$sigma[m_non_CG_low_model], m2=m_non_CG_mixmdl$mu[m_non_CG_high_model], sd2=m_non_CG_mixmdl$sigma[m_non_CG_high_model], p1=m_non_CG_mixmdl$lambda[m_non_CG_low_model], p2=m_non_CG_mixmdl$lambda[m_non_CG_high_model])
	segment_m_non_CG_cutoff
	# Schmitz data: 0.05920002 but 3 Gaussians doesn't look like a great fit
	# Becker data: 

	pdf(paste0(project_id,"_",meth_context,"_m_non_CG_segmentation_m_non_CG_cutoff_density_fitted.pdf"))
	print(plot(m_non_CG_mixmdl,which=2, breaks=seq(0,1,0.001)))
	dev.off()

	# Segment based on optimum parameters (m~0.067)
	# Clean up memory space
	rm(UMRLMRsegments_non_CG.gr)
	gc()

	m_non_CG.sel = segment_m_non_CG_cutoff
	# Run optimised segmentation
	UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_parents_CHG_CHH.gr, meth.cutoff=m_non_CG.sel, nCpG.cutoff=best_n_non_CG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

	# Save segmentation
	saveUMRLMRSegments(segs=UMRLMRsegments_non_CG.gr, GRangesFilename=paste0("UMRsLMRs_parents_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".tsv"))

	

	# Fitting Gaussians to separate the components of the pmeth distribution was generally unsuccessful.  It was able to identify the cutoff around pmeth=0.067, but this led to highly fragmented unmethylated regions, with very many very short methylated segments that were unconvincing.  Manually setting the cutoff to pmeth=0.15 in the segmentation, on the basis of visual inspection of the methylation landscape, produced a much better segmentation result.

	m_non_CG.sel=0.15

	# Clean up memory space
	rm(UMRLMRsegments_non_CG.gr)
	gc()
	
	UMRLMRsegments_non_CG.gr <- segmentUMRsLMRs(m=meth_parents_CHG_CHH.gr, meth.cutoff=m_non_CG.sel, nCpG.cutoff=best_n_non_CG, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_non_CG_segmentation_landscape_",best_n_non_CG,"_",m_non_CG.sel,".pdf"))
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_CG_segmentation_landscape_",best_n_non_CG,"_",m_non_CG.sel,"_density.pdf"))
	print(ggplot(as.data.frame(UMRLMRsegments_non_CG.gr), aes(x=pmeth)) + geom_histogram(binwidth=0.01))
	dev.off()

	# Save segmentation
	saveUMRLMRSegments(segs=UMRLMRsegments_non_CG.gr, GRangesFilename=paste0("UMRsLMRs_parents_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_non_CG_",m_non_CG.sel,"_",best_n_non_CG,".tsv"))
	
	
	#MG_segments.gr = setdiff(genome_loci.gr, mask_loci.gr)
	m_non_CG_segments.gr = setdiff(setdiff(genome_loci.gr, UMRLMRsegments_non_CG.gr), mask_loci.gr)
	UMR_non_CG_segments.gr = setdiff(UMRLMRsegments_non_CG.gr, mask_loci.gr)
	
	# Switch to CG context to complete the segmentation
	meth_context="CG"
	
	# Free up some space
	rm(meth_parents_CHG_CHH.gr)
	gc()
	
	# Read the CG methylome data in as a Granges object
	meth_parents_CG.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)

	# Read the CHG methylome data in as a Granges object
	#meth_parents_CHG.gr <- readMethylome(FileName=paste0(project_id,"_",meth_context,"_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	
	# Identify the parts of the parental mCGome which overlap with unmethylated parental mCHGome
	#potential_CG_segments.gr = meth_parents_CG.gr[unique(subjectHits(findOverlaps(UMRLMRsegments_CHG.gr, meth_parents_CG.gr)))]
	potential_CG_segments.gr = meth_parents_CG.gr[unique(subjectHits(findOverlaps(UMR_non_CG_segments.gr, meth_parents_CG.gr)))]
	
	# Try m.sel=0.95 to check out whole methylome landscape
	m.sel=0.95
	n.sel=1  # We segment using minimum segment size of 1 so that we can extend our mCHG segments by a single CG site only, if necessary (minimise likelihood of extending a mCHG segment into an adjacent mCG but non-mCHG segment).
	UMRLMRsegments_CG.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_CHG_segmentation_landscape_",n.sel,"_",m.sel,".pdf"))
	
	# Save a copy for visualisation
	saveUMRLMRSegments(segs=UMRLMRsegments_CG.gr, GRangesFilename=paste0("UMRsLMRs_parents_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_",m.sel,"_",n.sel,".tsv"))

	# Visual inspection of landscape shows clear cutoff at 0.4 between methylated and unmethylated segments among non-mCHG loci
	m.sel=0.4
	n.sel=1  # We segment using minimum segment size of 1 so that we can extend our m_non_CG segments by a single CG site only, if necessary (minimise likelihood of extending a m_non_CG segment into an adjacent mCG but non-m_non_CG segment).
	UMRLMRsegments_CG.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=UMRLMRsegments_CG.gr, GRangesFilename=paste0("UMRsLMRs_parents_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_",m.sel,"_",n.sel,".tsv"))


	# Find the mCG segments, rather than the unmethylated, and filter by masked genome:
	# First numbers = with 'Unknown' genes masked out, second numbers with only CHRM and CHRC masked
	mCG_segments_1.gr = setdiff(setdiff(genome_loci.gr, mask_loci.gr),UMRLMRsegments_CG.gr)
	length(mCG_segments_1.gr)
	# 41863, 41808 mCG segments

	# Adjust segments to remove any segments which don't overlap with any CG sites
	mCG_segments_1.gr = mCG_segments_1.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, mCG_segments_1.gr)))]
	length(mCG_segments_1.gr)
	# 41844, 41805 mCG segments
	
	# Add back the mCG values for the mCG segments, so it will save properly
	values(mCG_segments_1.gr) = cbind(values(mCG_segments_1.gr),meth_by_segment(meth_parents_CG.gr, segment_model=mCG_segments_1.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=mCG_segments_1.gr, GRangesFilename=paste0("UMRsLMRs_parents_mCG_masked_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_mCG_masked_",m.sel,"_",n.sel,".tsv"))
	
	# Find the mCG segments which overlap with m_non_CG segments:
	#mCHG_mCG_segments.gr=mCG_segments_1.gr[unique(subjectHits(findOverlaps(setdiff(setdiff(genome_loci.gr, mask_loci.gr),UMRLMRsegments_CHG.gr), mCG_segments_1.gr)))]
	m_non_CG_mCG_segments.gr=mCG_segments_1.gr[unique(subjectHits(findOverlaps(m_non_CG_segments.gr, mCG_segments_1.gr)))]
	length(m_non_CG_mCG_segments.gr)
	# 2683, 2584 mCG segments have an overlap with m_non_CG segments

	# Add back the mCG values for the mCG overlap segments, so it will save properly
	values(m_non_CG_mCG_segments.gr) = cbind(values(m_non_CG_mCG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=m_non_CG_mCG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=m_non_CG_mCG_segments.gr, GRangesFilename=paste0("UMRsLMRs_parents_non_CG_CG_overlap_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_non_CG_CG_overlap_",m.sel,"_",n.sel,".tsv"))

	# Take the Union of the m_non_CG best fit segments (subtracted from masked genome) and those mCG 1 0.4 segments which overlap them, to extend the m_non_CG best fit segments to the next closest CG site where appropriate:
	m_non_CG_segments.gr = reduce(union(setdiff(setdiff(genome_loci.gr, UMR_non_CG_segments.gr), mask_loci.gr), m_non_CG_mCG_segments.gr))
	length(m_non_CG_segments.gr)
	# 8492, 8423 segments
	
	# Adjust segments to remove any segments which don't overlap with any CG sites
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, m_non_CG_segments.gr)))]
	length(m_non_CG_segments.gr)
	# 8286, 8196 m_non_CG_segments

	# Add back the mCG values for the mCHG overlap segments, so it will save properly
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=m_non_CG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Save a copy for visualisation
	saveUMRLMRSegments(segs=m_non_CG_segments.gr, GRangesFilename=paste0("UMRsLMRs_parents_non_CG_segments.gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_non_CG_.tsv"))
	
	# Rebuild the non-m_non_CG fraction of the masked genome for subsequent segmentation in the CG context
	potential_CG_segments.gr = meth_parents_CG.gr[unique(subjectHits(findOverlaps(setdiff(setdiff(genome_loci.gr, m_non_CG_segments.gr), mask_loci.gr), meth_parents_CG.gr)))]

	### This is where we may have optimised the minimum segment length for segmenting the mCGome.  This code moved to separate script, and the numbers set here:
	
	# Schmitz numbers:
	best_n = 9
	prev_best_auc = 0.936370415538498
			
			
	# Execute the preferred segmentation model, visualise and save the results
	#UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=0.012, nCpG.cutoff=9, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape.pdf"))

	
	# Here are a couple of other ideas for setting cutoffs between U and M, especially for short segments:
	#m.sel=0.025
	# Set m.sel to 99ile of C/(C+T) for "U" sites
	m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.99)
	cat(paste0("mCG proportion which captures 99% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 99% of sites 'U' across all samples is 0.0769230769230769.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 99% of sites 'U' across all samples is 0.0833333333333333.

	# Using pmeth, rather than median.meth with .99 threshold for unmethylated segments produces some bleed into UMR from SMR particularly.  .95 is a better cutoff:
	m.sel=quantile(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C/(cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_C+cbind(coverage_data,meth_data)[(!is.na(coverage_data$Cov_C)) & (meth_data$meth_status=="U"),]$Cov_T),.95)
	cat(paste0("mCG proportion which captures 95% of sites 'U' across all samples is ",m.sel,".\n"))
	# Schmitz: mCG proportion which captures 95% of sites 'U' across all samples is 0.0377358490566038.
	# Fishers p<0.05, Binomial p<0.005: mCG proportion which captures 95% of sites 'U' across all samples is 0.0384615384615385.

	# Finally, we chose m.sel iteratively.  By building a finished model of long, medium and short segments in CG context, and assessing density distribution of mCG per segment, we identified a lot of short, lowly methylated segments, that are not credible.  We fit a mixture of Gaussians to this distribution and use their crossing point to set the mCG threshold here:
	
	m.sel=0.09309094
	
	# Segment mCG using the minimum segment which maximises AUC for ROC curve
	n.sel=best_n
	
	UMRLMRsegments.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_parents_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename=paste0("UMRsLMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))

	# We observe that n.sel=9 generates reasonable segments, not breaking up the majority of GBM loci.  However, it does miss short isolated regions of mCG which look like real GBM.
	# As a strategy to address this we carry out a new segmentation with the same m.sel, but with n.sel=3 then again with n.sel=1
	# We will use these finer segmentations to split Unmethylated segments, where necessary, but will not use resulting short unmethylated segments to further divide mCG segments.

	# Segment using minimum segment length of 3
	n.sel=3
	UMRLMRsegments_medium.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_parents_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments_medium.gr, GRangesFilename=paste0("UMRsLMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))

	# Segment using minimum segment length of 1
	n.sel=1
	UMRLMRsegments_fine.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, nCpG.smoothing = 1, myGenomeSeq=Athaliana, seqLengths=sLengths, pdfFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segmentation_landscape_parents_CG_",m.sel,"_",n.sel,".pdf"))
	#head(UMRLMRsegments.gr)
	#plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel)

	saveUMRLMRSegments(segs=UMRLMRsegments_fine.gr, GRangesFilename=paste0("UMRsLMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))

	# We have generated Unmethylated Regions with minimum lengths of 9, 3 and 1 CG sites respectively (UMR9, UMR3 and UMR1)
	# Masked Genome (MG) = Genome - Mask
	# Long Methylated regions (LMR) = MG - (UMR9 + mCHG)
	# Medium Methylated regions (MMR) = (MG - (UMR3 + mCHG)) - LMR
	# Short methylated regions (SMR) = (MG - (UMR1 + mCHG)) - (LMR + MMR)
	# Unmethylated regions (UMR) = UMR9 - (SMR + MMR +mCHG)
	# Methylated regions (MR) = MG - UMR

	# Load parental CHG and CHH methylomes separately, ready to use for annotating segments
	meth_parents_CHG.gr <- readMethylome(FileName=paste0(project_id,"_CHG_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)
	# Now the SLOW bit
	meth_parents_CHH.gr <- readMethylome(FileName=paste0(project_id,"_CHH_TM_read_counts_across_parents.tsv"), seqLengths=sLengths)

	
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
	set_aside_no_CG.gr = c(setdiff(LMR_segments.gr, LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, LMR_segments.gr)))]), setdiff(MMR_segments.gr, MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MMR_segments.gr)))]), setdiff(SMR_segments.gr, SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, SMR_segments.gr)))]))
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 45
	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 1362, 1065, 1065, 1029, 5
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 16500, 16909, 13659, 12401, 12419
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 5294, 5315, 5898, 5053, 5069
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 5635, 8385, 9949, 7955, 7648
	#MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, MR_segments.gr)))]
	#cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 27429, 30609, 29506, 25106

	
	# Annotate segments with CG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=LMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","LMR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MMR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=SMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","SMR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=MR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#replace values previously added: 
	values(m_non_CG_segments.gr) = cbind(meth_by_segment(meth_parents_CG.gr, segment_model=m_non_CG_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m_non_CG.sel), mMeth.classes=c("UMR","TEM")))

	# Adjust SMRs to remove any with too low mCG
	# Visual observation of the distribution of mCG density, and good sense (we are averaging over six samples) suggests SMRs with mCG<0.25 are not credible
	SMR_segments.gr = SMR_segments.gr[ as.vector((!is.na(SMR_segments.gr$pmeth)) & (SMR_segments.gr$pmeth>0.25)) ]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after removing segments with mCG<0.25.\n"))
	# 4619, 5255, 4570, 4562, 4564

	# Build set of UMRs corresponding to new set of LMRs, MMRs, SMRs, TEMs, set_asides:
	UMR_segments.gr = setdiff(setdiff(setdiff(setdiff(setdiff(setdiff(genome_loci.gr, SMR_segments.gr), MMR_segments.gr), LMR_segments.gr), m_non_CG_segments.gr), set_aside_no_CG.gr), mask_loci.gr)
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments.\n"))
	# 26495, 26255, 26178

	# Remove any UMR segments containing no CG sites, after setting them aside
	set_aside_no_CG.gr = c(set_aside_no_CG.gr, setdiff(UMR_segments.gr, UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, UMR_segments.gr)))]))
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no ",meth_context," sites.\n"))
	# 27188, 30383, 29026, 26389, 26137
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 86
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CG.gr, segment_model=UMR_segments.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	
	# Adjust segments to remove any segments which don't overlap with any CHG sites, after setting them aside
	set_aside_no_CHG.gr = c(setdiff(LMR_segments.gr, LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, LMR_segments.gr)))]), setdiff(MMR_segments.gr, MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MMR_segments.gr)))]), setdiff(SMR_segments.gr, SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, SMR_segments.gr)))]))
	cat(paste0("set_aside_no_CHG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 1099

	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 1344, 1047, 1047, 1012, 5
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 16483, 16892, 13621, 12304, 12329
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 5205, 5221, 5750, 4873, 4893
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 4677, 6672, 7964, 3797, 3731
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 27185, 30374, 29000, 26268, 26087
	#MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, MR_segments.gr)))]
	#cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	# 26365, 28785, 27335, 23169, 
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, m_non_CG_segments.gr)))]
	cat(paste0("m_non_CG_segments.gr contains ",length(m_non_CG_segments.gr)," segments after trimming segments with no CHG sites.\n"))
	#  ,  , 4437, 8193, 8117
	
	# Annotate segments with CHG methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=LMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=SMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=UMR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=MR_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=m_non_CG_segments.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))


	# Adjust segments to remove any segments which don't overlap with any CHH sites, after setting them aside
	set_aside_no_CHH.gr = c(setdiff(LMR_segments.gr, LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, LMR_segments.gr)))]), setdiff(MMR_segments.gr, MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, MMR_segments.gr)))]), setdiff(SMR_segments.gr, SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, SMR_segments.gr)))]))
	cat(paste0("set_aside_no_CHG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 1099

	MG_segments.gr = MG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, MG_segments.gr)))]
	cat(paste0("MG_segments.gr contains ",length(MG_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 1012
	LMR_segments.gr = LMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, LMR_segments.gr)))]
	cat(paste0("LMR_segments.gr contains ",length(LMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 12328
	MMR_segments.gr = MMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, MMR_segments.gr)))]
	cat(paste0("MMR_segments.gr contains ",length(MMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 4890
	SMR_segments.gr = SMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, SMR_segments.gr)))]
	cat(paste0("SMR_segments.gr contains ",length(SMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 3720
	UMR_segments.gr = UMR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, UMR_segments.gr)))]
	cat(paste0("UMR_segments.gr contains ",length(UMR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 26087
	#MR_segments.gr = MR_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, MR_segments.gr)))]
	#cat(paste0("MR_segments.gr contains ",length(MR_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 23116
	m_non_CG_segments.gr = m_non_CG_segments.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, m_non_CG_segments.gr)))]
	cat(paste0("m_non_CG_segments.gr contains ",length(m_non_CG_segments.gr)," segments after trimming segments with no CHH sites.\n"))
	# 8117
	
	# Annotate segments with CHH methylation values
	#values(MG_segments.gr) = cbind(values(MG_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=MG_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(LMR_segments.gr) = cbind(values(LMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=LMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(MMR_segments.gr) = cbind(values(MMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=MMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(SMR_segments.gr) = cbind(values(SMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=SMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(UMR_segments.gr) = cbind(values(UMR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=UMR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	#values(MR_segments.gr) = cbind(values(MR_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=MR_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(m_non_CG_segments.gr) = cbind(values(m_non_CG_segments.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=m_non_CG_segments.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

	# Merge the various distinct classes of segment
	segmentation_model.gr = c(LMR_segments.gr, MMR_segments.gr, SMR_segments.gr, UMR_segments.gr, m_non_CG_segments.gr)
	length(segmentation_model.gr)
	# 55873, 55142

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
	
	saveSegmentationModel(segs=LMR_segments.gr, GRangesFilename=paste0("LMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_LMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=MMR_segments.gr, GRangesFilename=paste0("MMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=SMR_segments.gr, GRangesFilename=paste0("SMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_SMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=UMR_segments.gr, GRangesFilename=paste0("UMRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_UMR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))
	saveSegmentationModel(segs=MR_segments.gr, GRangesFilename=paste0("MRs_parents_CG_",m.sel,"_",n.sel,".gr.rds"), TableFilename=paste0(project_id,"_",meth_context,"_MethylSeekR_MR_segments_parents_CG_",m.sel,"_",n.sel,".tsv"))



	# Plot mCHG vs mCG for each mCG class of segment
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_vs_mCG_median_by_mCG_class.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=type)) + theme_minimal())
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCHG_vs_mCG_mean_by_mCG_class.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=type)) + theme_minimal())
	dev.off()

	# Plot density of segment mCHG 
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_median_mCHG_density.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=median.meth.1)) + geom_density(aes(colour=type)) + theme_minimal())
	dev.off()
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mean_mCHG_density.pdf"))
	print(ggplot(data.frame(values(segmentation_model.gr)@listData), aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal())
	ggplot(data.frame(values(segmentation_model.gr)@listData)[data.frame(values(segmentation_model.gr)@listData)$type=="LMR",], aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal()
	dev.off()
	
	# Show a table of counts of segments by type

	#table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),]$type)
	table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$type)
	# Schmitz data:
	#LMR   MMR   SMR   UMR 
	#16546  5201  6769 30811 
	#16508  5205     4  6668 30643

	# Using TEM segments separately derived from segmenting mCHG:
	#LMR   MMR    MR   SMR   TEM   UMR 
	#13146  5726     1  2940  4343 25645 
	# Using TEM segments separately derived from CHG/CHH combined:
	#11813  4850     9  3724  7823 27519 
	# Not masking 'Unknown' genes
	#12081  4843    28  3719  7565 26900 

	# 4,1 segments ended up with type="MR" (short segments with relatively high levels of methylation) so we set these to be "SMR":
	segmentation_model.gr@elementMetadata@listData$type[(!is.na(segmentation_model.gr@elementMetadata@listData$type)) & (segmentation_model.gr@elementMetadata@listData$type=="MR")]="SMR"

	
	# Show a table of counts of segments by type

	table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$type)
	#Schmitz data:
	#LMR   MMR   SMR   TEM   UMR 
	#11813  4850  3733  7823 27519 
	#12081  4843  3747  7565 26900 
	
	# We set mCG_cutoff_mean here to the value previously found, and reiterate mCHG cutoffs to be used later for classifying finished segments
	mCG_cutoff_mean = 0.09309094    # derived by fitting mixture of Gaussians to density distribution of segment mCG of whole parental methylome segmented in CG context
	m.sel = 0.09309094
	segment_mCHG_cutoff = 0.1262584 # derived by fitting mixture of Gaussians to density distribution of segment mCHG of whole parental methylome segmented in CHG context
	mCHG_cutoff_mean = 0.1262584
	mCHG_cutoff_median = 0.1262584
	segment_mCHH_cutoff = 0.1      # this is somewhat arbitrary, and is only used to classify the segments as UMR/MR in CHH.  May come back to this, and give it a better cutoff later, if this is thought to be useful
	mCHH_cutoff_mean = 0.1      
	

	# Plot mCG vs. segment size, with segments coloured by type 
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCG_vs_length_by_type_revised.pdf"))
	
	# Plot distribution of segment mCG and mCHG by type
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCG"))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=type), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=type), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.1)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCHG"))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHG"))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=values(segmentation_model.gr)@listData$median.meth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.2)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.1)) + geom_density(aes(colour=type)) + theme_minimal() + facet_wrap(~ type, ncol=1, scales="free_y") + xlab("Segment Mean mCHH"))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHH_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.2)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$pmeth.2>mCHH_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.2)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHH"))
	
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$pmeth.2>mCHH_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth.1, y=pmeth.2)) + geom_point(aes(colour=type), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCHG") + ylab("Segment Mean mCHH"))

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
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1,mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=type), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	dev.off()
	
	
	# Crystalise segment methylation status as a column in the segmentation object
	# mCHG cutoff is segment_mCHG_cutoff
	# mCHH cutoff is segment_mCHH_cutoff
	# mCG cutoff is m.sel
	segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse((cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff) | (cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[,]$mCHH>segment_mCHH_cutoff),"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$pmeth>m.sel,"GBM","UMR"))

	# For Schmitz data, two of these segments have NA pmeth values
	#segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$mCHG>mCHG_cutoff_mean,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$pmeth>m.sel,"GBM","UMR"))

	# Set aside any unclassified segments
	length(segmentation_model.gr)
	# 55873, 55142
	
	set_aside_no_class.gr = segmentation_model.gr[is.na(as.data.frame(segmentation_model.gr)$segment.mStatus)]
		
	segmentation_model.gr = segmentation_model.gr[!is.na(as.data.frame(segmentation_model.gr)$segment.mStatus)]
	length(segmentation_model.gr)
	# 55867, 55136

	cat(paste0("There are ",length(set_aside_no_class.gr), " segments with no methylation status classification.\n"))
	# Schmitz data: 6
	
	# It's not exactly clear why this is not zero - if these segments don't overlap CHG or CHH sites they should have been already set aside.
	# On inspection, these segments all have mCG data (in Schmitz data set)
	# Check for overlaps with TEs or exons and characterise accordingly

	# Load TE and exon annotations so we can check for overlaps to resolve anomalies
	mixedCaseChr <- function(s, strict = FALSE) {
		paste0(toupper(substr(s,1,1)),tolower(substr(s,2,3)),toupper(substr(s,4,4)))
	}

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

	# For Schmitz data the 6 segments found all had mCG but various combos of NaN and <NA> values for mCHG and mCHH
	# These most likely arose at the end of GBM segments. Check overlaps with TE or exons to decide whether TEM or GBM if mCG>m.sel

	# Set segments overlapping exons to "GBM" if mCG supports it
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))])@listData$pmeth>m.sel,"GBM","UMR")
	# Set segments overlapping transposons to "TEM" if mCG supports it (to override exons in case segment overlaps both)
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))])@listData$pmeth>m.sel,"TEM","UMR")
	# Presumption of "UMR" for segments which don't overlap exons or TEs
	set_aside_no_class.gr[is.na(as.data.frame(set_aside_no_class.gr)$segment.mStatus)]@elementMetadata@listData$segment.mStatus = "UMR"

	
	# Classify each of the other set_aside segments, where possible, and merge back into the model
	# Check first to make sure that none of the set_aside_segments.gr overlap with each other or with segments in the model (sense check)

	length(c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr))
	# 1200
	length(reduce(c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr)))
	# 1200
	length(findOverlaps(c(set_aside_no_CG.gr, set_aside_no_CHG.gr, set_aside_no_CHH.gr), segmentation_model.gr))
	#0
	
	# Split set_aside_no_*.gr into those which have data in other contexts and those that don't
	
	set_aside_no_CG_CHG.gr = setdiff(set_aside_no_CG.gr, set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, set_aside_no_CG.gr)))])
	length(set_aside_no_CG_CHG.gr)
	# 35
	
	set_aside_no_CG_CHG_CHH.gr = setdiff(set_aside_no_CG_CHG.gr, set_aside_no_CG_CHG.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, set_aside_no_CG_CHG.gr)))])
	length(set_aside_no_CG_CHG_CHH.gr)
	# 2  
	# These cannot be annotated - they both look like TEM by eye
	set_aside_no_CG_CHG.gr = set_aside_no_CG_CHG.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, set_aside_no_CG_CHG.gr)))]
	length(set_aside_no_CG_CHG.gr)
	# 33
	# These can be annotated using CHH
	
	set_aside_no_CG.gr = set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, set_aside_no_CG.gr)))]
	length(set_aside_no_CG.gr)
	# 51
	
	set_aside_no_CG_CHH.gr = setdiff(set_aside_no_CG.gr, set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, set_aside_no_CG.gr)))])
	length(set_aside_no_CG_CHH.gr)
	# 0
	set_aside_no_CG.gr = set_aside_no_CG.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, set_aside_no_CG.gr)))]
	length(set_aside_no_CG.gr)
	# 51
	# These can be annotated using CHG & CHH
	
	set_aside_no_CHG_CG.gr = setdiff(set_aside_no_CHG.gr, set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, set_aside_no_CHG.gr)))])
	length(set_aside_no_CHG_CG.gr)
	# 0
	set_aside_no_CHG.gr = set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, set_aside_no_CHG.gr)))]
	length(set_aside_no_CHG.gr)
	# 1099
	
	set_aside_no_CHG_CHH.gr = setdiff(set_aside_no_CHG.gr, set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, set_aside_no_CHG.gr)))])
	length(set_aside_no_CHG_CHH.gr)
	# 12
	# These can be annotated using CG
	
	set_aside_no_CHG.gr = set_aside_no_CHG.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, set_aside_no_CHG.gr)))]
	length(set_aside_no_CHG.gr)
	# 1087
	# These can be annotated using CG and CHH, although CHH is unreliable as a differentiator between GBM and TEM so also use annotation. Chances are though that any anomalies will be resolved by the merging of adjacent TEMs with short GBMs

	set_aside_no_CHH_CG.gr = setdiff(set_aside_no_CHH.gr, set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, set_aside_no_CHH.gr)))])
	length(set_aside_no_CHH_CG.gr)
	# 0
	set_aside_no_CHH.gr = set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, set_aside_no_CHH.gr)))]
	length(set_aside_no_CHH.gr)
	# 15
	
	set_aside_no_CHH_CHG.gr = setdiff(set_aside_no_CHH.gr, set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, set_aside_no_CHH.gr)))])
	length(set_aside_no_CHH_CHG.gr)
	# 0
	set_aside_no_CHH.gr = set_aside_no_CHH.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, set_aside_no_CHH.gr)))]
	length(set_aside_no_CHH.gr)
	# 15
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
		values(set_aside_no_CG.gr) = cbind(values(set_aside_no_CG.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=set_aside_no_CG.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, 	mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
		values(set_aside_no_CG.gr) = cbind(values(set_aside_no_CG.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=set_aside_no_CG.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
		
		# no CG sites so segments likely arise at the end of a TEM segment.  Use mCHG and mCHH to decide status
		set_aside_no_CG.gr@elementMetadata@listData$segment.mStatus = ifelse((values(set_aside_no_CG.gr)@listData$pmeth.1>segment_mCHG_cutoff) | (values(set_aside_no_CG.gr)@listData$pmeth.2>segment_mCHH_cutoff),"TEM","UMR")
	}
	
	if (length(set_aside_no_CG_CHG.gr) > 0) {
		values(set_aside_no_CG_CHG.gr) = cbind(NA_by_segment(set_aside_no_CG_CHG.gr), NA_by_segment(set_aside_no_CG_CHG.gr))
		values(set_aside_no_CG_CHG.gr) = cbind(values(set_aside_no_CG_CHG.gr), meth_by_segment(meth_parents_CHH.gr, segment_model=set_aside_no_CG_CHG.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

		# For this one no CG so it arises at the end of a TEM segment, most likely.  Use mCHH to decide (though may be safer to refer to annotation, as CHH in absence of CHG may be unreliable)
		set_aside_no_CG_CHG.gr@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_CG_CHG.gr)@listData$pmeth.2>segment_mCHH_cutoff,"TEM","UMR")
	}
	
	if (length(set_aside_no_CG_CHH.gr) > 0) {
		values(set_aside_no_CG_CHG.gr) = cbind(NA_by_segment(set_aside_no_CG_CHH.gr), meth_by_segment(meth_parents_CHG.gr, segment_model=set_aside_no_CG_CHH.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")), NA_by_segment(set_aside_no_CG_CHH.gr))

		# For this one no CG so it arises at the end of a TEM segment, most likely.  Use mCHG to decide
		set_aside_no_CG_CHH.gr@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_CG_CHH.gr)@listData$pmeth.1>segment_mCHG_cutoff,"TEM","UMR")
	}
	
	if (length(set_aside_no_CHG_CHH.gr) > 0) {
		values(set_aside_no_CHG_CHH.gr) = meth_by_segment(meth_parents_CG.gr, segment_model=set_aside_no_CHG_CHH.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))
		values(set_aside_no_CHG_CHH.gr) = cbind(values(set_aside_no_CHG_CHH.gr), NA_by_segment(set_aside_no_CHG_CHH.gr), NA_by_segment(set_aside_no_CHG_CHH.gr))

		# Check overlaps with TE or exons to decide whether TEM or GBM if mCG>m.sel
		#set_aside_exon_hits = unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG_CHH.gr)))
		#set_aside_transposon_hits = unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG_CHH.gr)))

		set_aside_no_CHG_CHH.gr@elementMetadata@listData$segment.mStatus = as.character(rep(NA, length(set_aside_no_CHG_CHH.gr)))
		# Set all those overlapping exons to "GBM"
		set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG_CHH.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG_CHH.gr)))])@listData$pmeth>m.sel,"GBM","UMR")
		# Set all those overlapping transposons to "TEM" (to override exons in case overlaps both)
		set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG_CHH.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_CHG_CHH.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG_CHH.gr)))])@listData$pmeth>m.sel,"TEM","UMR")
		# For those overlapping neither exons nor TEs, start with a presumption of "UMR" (why?)
		set_aside_no_CHG_CHH.gr[is.na(as.data.frame(set_aside_no_CHG_CHH.gr)$segment.mStatus)]@elementMetadata@listData$segment.mStatus = rep("UMR", length(set_aside_no_CHG_CHH.gr[is.na(as.data.frame(set_aside_no_CHG_CHH.gr)$segment.mStatus)]))
	}

	if (length(set_aside_no_CHG.gr) > 0) {
		values(set_aside_no_CHG.gr) = meth_by_segment(meth_parents_CG.gr, segment_model=set_aside_no_CHG.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))
		values(set_aside_no_CHG.gr) = cbind(values(set_aside_no_CHG.gr), NA_by_segment(set_aside_no_CHG.gr))
		values(set_aside_no_CHG.gr) = cbind(values(set_aside_no_CHG.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=set_aside_no_CHG.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))

		# For this one no CHG so need to check CG and CHH
		set_aside_no_CHG.gr@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_CHG.gr)@listData$pmeth.2>segment_mCHH_cutoff,"TEM",ifelse(values(set_aside_no_CHG.gr)@listData$pmeth>m.sel,"GBM","UMR"))
		
		# In 2 cases (Schmitz) the resulting classification failed with too little coverage in CHH
		# Set all those overlapping exons to "GBM"
		set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)][unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)])))]$segment.mStatus = ifelse(values(set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)][unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)])))])@listData$pmeth>m.sel,"GBM","UMR") 
		# Set all those overlapping transposons to "TEM" (to override exons in case overlaps both)
		set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)][unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)])))]$segment.mStatus = ifelse(values(set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)][unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)])))])@listData$pmeth>m.sel,"GBM","UMR") 

		# For those overlapping neither exons nor TEs, start with a presumption of "UMR" (why?)
		set_aside_no_CHG.gr[is.na(as.data.frame(set_aside_no_CHG.gr)$segment.mStatus)]@elementMetadata@listData$segment.mStatus = "UMR"
	}
	
	if (length(set_aside_no_CHH.gr) > 0) {
		values(set_aside_no_CHH.gr) = meth_by_segment(meth_parents_CG.gr, segment_model=set_aside_no_CHH.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR"))
		values(set_aside_no_CHH.gr) = cbind(values(set_aside_no_CHH.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=set_aside_no_CHH.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
		values(set_aside_no_CHH.gr) = cbind(values(set_aside_no_CHH.gr), NA_by_segment(set_aside_no_CHH.gr))

		# For this one no CHH so need to check CG and CHG and whether overlaps exon or TE
		set_aside_no_CHH.gr@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_CHH.gr)@listData$pmeth.1>segment_mCHG_cutoff,"TEM",ifelse(values(set_aside_no_CHH.gr)@listData$pmeth>m.sel,"GBM","UMR"))
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
	# 54028, 54341

	# Adjust segments to set aside any segments which don't overlap with any CG sites. We will deal with them later
	# Counts here are first masking only CHRM and CHRC
	set_aside_no_CG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 64
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 54301
	set_aside_no_CHG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CHG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 1048
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 53253
	set_aside_no_CHH.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CHH.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 15
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 53238


	# Reannotate the segmentation model (again!)
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_parents_CG.gr, segment_model=segmentation_model.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=segmentation_model.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHG_cutoff), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=segmentation_model.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHH_cutoff), mMeth.classes=c("UMR","MR")))
	segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$pmeth>m.sel,"GBM","UMR"))
	
	# Check whether any unannotated segments remain in the model
	set_aside_no_class.gr = segmentation_model.gr[is.na(as.data.frame(segmentation_model.gr)$segment.mStatus)]
	length(set_aside_no_class.gr)
	# 4
	segmentation_model.gr = segmentation_model.gr[!is.na(as.data.frame(segmentation_model.gr)$segment.mStatus)]
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
	merged_segments.gr = reduce(c(segmentation_model.gr[segmentation_model.gr$segment.mStatus == "TEM"], segmentation_model.gr[(segmentation_model.gr$segment.mStatus == "GBM") & (as.data.frame(segmentation_model.gr)$width < GBM_distinct_min_width)]))
	length(merged_segments.gr)
	# 16195, 16527

	# remove tha annotations from the segmentation model remainder, so it will merge properly with the new, unannotated, merged segments
	segmentation_model_remainder.gr = segmentation_model.gr[!((segmentation_model.gr$segment.mStatus == "TEM") | ((segmentation_model.gr$segment.mStatus == "GBM") & (as.data.frame(segmentation_model.gr)$width < GBM_distinct_min_width)))]
	length(segmentation_model_remainder.gr)
	# 36042, 35382
	segmentation_model_remainder.gr@elementMetadata@listData = list()

	# put the whole model back together again by concatenating the new merged segments with all the other segments from the original model which were not considered for merging
	# May as well sort the new model as we build it, in case we need a sorted model later
	segmentation_model.gr = sort(sortSeqlevels(c(merged_segments.gr, segmentation_model_remainder.gr)))
	length(segmentation_model.gr)
	# 52237, 51909

	# Adjust segments to set aside any segments which don't overlap with any CG sites. We will deal with them later
	# Counts here are first masking only CHRM and CHRC
	set_aside_no_CG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CG.gr)," segments which contain no CG sites.\n"))
	# 40
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 51869
	set_aside_no_CHG.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CHG.gr)," segments which contain no CHG sites.\n"))
	# 504
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHG.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 51365
	set_aside_no_CHH.gr = setdiff(segmentation_model.gr,segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, segmentation_model.gr)))])
	cat(paste0("set_aside_no_CG.gr contains ",length(set_aside_no_CHH.gr)," segments which contain no CHH sites.\n"))
	# 11
	segmentation_model.gr = segmentation_model.gr[unique(subjectHits(findOverlaps(meth_parents_CHH.gr, segmentation_model.gr)))]
	length(segmentation_model.gr)
	# 51354


	
	# Reannotate the segmentation model (again!)
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_parents_CG.gr, segment_model=segmentation_model.gr, meth_context="CG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(m.sel), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_parents_CHG.gr, segment_model=segmentation_model.gr, meth_context="CHG", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHG_cutoff), mMeth.classes=c("UMR","MR")))
	values(segmentation_model.gr) = cbind(values(segmentation_model.gr),meth_by_segment(meth_parents_CHH.gr, segment_model=segmentation_model.gr, meth_context="CHH", myGenomeSeq=Athaliana, seqLengths=sLengths, mMeth.classification=c(segment_mCHH_cutoff), mMeth.classes=c("UMR","MR")))
	#segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$mCHG>segment_mCHG_cutoff,"TEM",ifelse(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[,]$pmeth>m.sel,"GBM","UMR"))
	segmentation_model.gr@elementMetadata@listData$segment.mStatus = ifelse((values(segmentation_model.gr)@listData$pmeth.1>segment_mCHG_cutoff) | (values(segmentation_model.gr)@listData$pmeth.2>segment_mCHH_cutoff),"TEM",ifelse(values(segmentation_model.gr)@listData$pmeth>m.sel,"GBM","UMR"))

	# Add back in any of the original set aside segments which don't overlap with the new model, and sort for good measure
	segmentation_model.gr = sort(sortSeqlevels(c(segmentation_model.gr, set_aside_segments.gr[!(set_aside_segments.gr %in% set_aside_segments.gr[unique(subjectHits(findOverlaps(segmentation_model.gr, set_aside_segments.gr)))])])))
	length(segmentation_model.gr)
	# 51909

	# Assign status to the set_aside_no_class.gr segments
	# Set segments overlapping exons to "GBM" if mCG supports it
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(exon_ranges, set_aside_no_class.gr)))])@listData$pmeth>m.sel,"GBM","UMR")
	# Set segments overlapping transposons to "TEM" if mCG supports it (to override exons in case segment overlaps both)
	set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))]@elementMetadata@listData$segment.mStatus = ifelse(values(set_aside_no_class.gr[unique(subjectHits(findOverlaps(transposon_ranges, set_aside_no_class.gr)))])@listData$pmeth>m.sel,"TEM","UMR")
	# Presumption of "UMR" for segments which don't overlap exons or TEs
	set_aside_no_class.gr[is.na(as.data.frame(set_aside_no_class.gr)$segment.mStatus)]@elementMetadata@listData$segment.mStatus = "UMR"


	# Set segments overlapping exons to "GBM" if mCG supports it
	segmentation_model.gr[is.na(as.data.frame(segmentation_model.gr)$segment.mStatus),] = set_aside_no_class.gr
	
	# Output a text version of the segmentation model for visualisation
	write.table(cbind("Chromosome"=as.character(as.data.frame(segmentation_model.gr)$seqnames),"Start"=as.data.frame(segmentation_model.gr)$start,"End"=as.data.frame(segmentation_model.gr)$end,"Type"=as.data.frame(segmentation_model.gr)$segment.mStatus),file=paste0(project_id,"_",meth_context,"_segmentation_model_draft.tsv"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

	# Save a copy in case we want to come back to it later
	saveUMRLMRSegments(segs=segmentation_model.gr, GRangesFilename=paste0(project_id,"_",meth_context,"_segmentation_model_draft.rds"), TableFilename=paste0(project_id,"_",meth_context,"_segmentation_model_draft.tsv"))


	table(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),]$segment.mStatus)
	# These are the counts of segments by type, based on coverage in all 3 contexts.  There are additional segments classified on the basis of two, or one contexts, plus annotation, which are included in the second table below.
	
	#Schmitz data:
	#GBM   TEM   UMR 
	#19738  7500 26313 
	#18292  7460 25609 (segmentation_model_draft files)
	
	# After annotation-driven reannotation of GBM segments to TEM if they don't overlap a gene (segmentation_model_draft2 files):
	# GBM   TEM   UMR 
	#16629  8793 25607 

	table(as.data.frame(segmentation_model.gr)$segment.mStatus)
	#18774  7484 25651  (segmentation_model_draft files)
	#17006  8919 25648  (segmentation_model_draft2 files)
	
	
	# Plot mCG vs. segment size, with segments coloured by segment methylation status
	pdf(paste0(project_id,"_",meth_context,"_genome_segment_mCG_vs_length_by_mStatus.pdf"))

	# Plot distribution of segment mCG, mCHG, mCHH by status
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Median mCG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Mean mCG"))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=1) + coord_trans(x = "log10") + theme_minimal() +xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 20000)))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.1)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.1)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Mean mCHG"))

	# Plot segment mCHG vs mCG
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$median.meth.1>mCHG_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHG"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=ifelse(values(segmentation_model.gr)@listData$pmeth.1>mCHG_cutoff_mean,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHG"))
	
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=values(segmentation_model.gr)@listData$median.meth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=median.meth.2)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=pmeth.2)) + geom_density(aes(colour=segment.mStatus)) + theme_minimal() + facet_wrap(~ segment.mStatus, ncol=1, scales="free_y") + xlab("Segment Mean mCHH"))

	# Plot segment mCHH vs mCG
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$median.meth.2>mCHH_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth, y=median.meth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCG") + ylab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$pmeth.2>segment_mCHH_cutoff,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth, y=pmeth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCG") + ylab("Segment Mean mCHH"))
	
	# Plot segment mCHH vs mCHG
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$median.meth.2>mCHH_cutoff_median,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=median.meth.1, y=median.meth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Median mCHG") + ylab("Segment Median mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$pmeth.2>segment_mCHH_cutoff,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth.1, y=pmeth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.5) + theme_minimal() + xlab("Segment Mean mCHG") + ylab("Segment Mean mCHH"))
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHH=ifelse(values(segmentation_model.gr)@listData$pmeth.2>segment_mCHH_cutoff,"High","Low"), seg_width=segmentation_model.gr@ranges@width), aes(x=pmeth.1, y=pmeth.2)) + geom_point(aes(colour=pmeth, size=(seg_width)/1000), alpha=0.5) + theme_minimal() + xlab("Segment Mean mCHG") + ylab("Segment Mean mCHH") + scale_colour_gradient("Segment\nMean mCG",low="green", high="red"))
	
	# Plot m* vs length by mStatus
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$median.meth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$median.meth)),], aes(x=seg_width, y=median.meth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>mCHG_cutoff_median,paste0("1. 'Transposon-like' segments:\n median mCHG > ",round(mCHG_cutoff_median, digits=4)),ifelse(median.meth>m.sel,paste0("2. 'Gene-body-like' segments:\n median mCG > ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)),paste0("3. Unmethylated segments:\n median mCG <= ",round(m.sel, digits=4),"\n median mCHG <= ",round(mCHG_cutoff_median, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Median mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))

	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>segment_mCHG_cutoff,paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHG>segment_mCHG_cutoff | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	#print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1,mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)),], aes(x=seg_width, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse(mCHH>segment_mCHH_cutoff,paste0("1. 'Transposon-like' segments:\n mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))

	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1,mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)),], aes(x=seg_width, y=pmeth)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))

	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1,mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)),], aes(x=seg_width, y=pmeth.1)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHG") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))
		
	print(ggplot(cbind(data.frame(values(segmentation_model.gr)@listData),mCHG=values(segmentation_model.gr)@listData$pmeth.1,mCHH=values(segmentation_model.gr)@listData$pmeth.2, seg_width=segmentation_model.gr@ranges@width)[(!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.1)) & (!is.na(data.frame(values(segmentation_model.gr)@listData)$pmeth.2)),], aes(x=seg_width, y=pmeth.2)) + geom_point(aes(colour=segment.mStatus), size=0.1, alpha=0.25) + coord_trans(x = "log10") + theme_minimal() + facet_wrap(~ ifelse((mCHG>segment_mCHG_cutoff) | (mCHH>segment_mCHH_cutoff),paste0("1. 'Transposon-like' segments:\n mean mCHG > ",round(segment_mCHG_cutoff, digits=4),"\n or mean mCHH > ",round(segment_mCHH_cutoff, digits=4)),ifelse(pmeth>m.sel,paste0("2. 'Gene-body-like' segments:\n mean mCG > ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)),paste0("3. Unmethylated segments:\n mean mCG <= ",round(m.sel, digits=4),"\n mean mCHG <= ",round(segment_mCHG_cutoff, digits=4),"\n and mean mCHH <= ",round(segment_mCHH_cutoff, digits=4)))), ncol=3) + xlab("Segment length (nt)") +ylab("Segment Mean mCHH") + scale_x_continuous(breaks = seq(0, 250000, by = 10000)))


	
	dev.off()

	