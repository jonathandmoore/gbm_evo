# parse all segmentation models, and extract 'haplotype block'-like merged model encompassing all samples

pathroot = "X:/Daniel-Zilberman/"
#setwd("W:/Daniel-Zilberman/Jay_1001_methylomes")
setwd(paste0(pathroot,"Projects/Jay-1001_methylomes/5-analysis"))
#setwd("/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis")
#setwd("../Desktop/1001")
 
library(GenomicRanges)
library(stringr)
library(ggplot2)

all_samples_meth_status = readRDS(file="m1001_CG_all_samples_meth_status_0.005_0.05_10_25_SNPs.rds")

model_dir = "segmentation_models"
# make an empty structure to hold the resulting model
merged_model = data.frame()
no_models = 0
ignore_list=c()

# load in all the segmentation models, and create the union of all the gbM segments in each gene

for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 10,19)
  if (substr(this_sample,10,10)=="_") {this_sample=substr(this_sample,1,9)}
  if (!(this_sample %in% ignore_list)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model = this_model[substr(this_model$Type,1,4)=="GBM_",]
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    if (no_models == 1) {
      overlap_model = this_model.gr
    } else {
      overlap_model = union(overlap_model, this_model.gr, ignore.strand=TRUE)
    }
  }
}

# Add "Chr" to the chromosome names
levels(overlap_model@seqnames@values) = paste0("Chr",as.character(levels(overlap_model@seqnames@values)))
overlap_model@seqinfo@seqnames = levels(overlap_model@seqnames@values)


model_dir = "segmentation_models_v2"
# make an empty structure to hold the resulting model
merged_model = data.frame()
no_models = 0
ignore_list=c()

# load in all the alternative, trimmed segmentation models, and create the union of all the modified gbM segments in each gene
# These alternative models were created using the script assess_segment_gene_collinearity.R

for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  #if (substr(this_sample,10,10)=="_") {this_sample=substr(this_sample,1,9)}
  if (!(this_sample %in% ignore_list)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    #this_model = this_model[substr(this_model$Type,1,4)=="GBM_",]
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    if (no_models == 1) {
      overlap_model = this_model.gr
    } else {
      overlap_model = union(overlap_model, this_model.gr, ignore.strand=TRUE)
    }
  }
}

# Add "Chr" to the chromosome names
levels(overlap_model@seqnames@values) = paste0("Chr",as.character(levels(overlap_model@seqnames@values)))
overlap_model@seqinfo@seqnames = levels(overlap_model@seqnames@values)


# version 3: load in all the alternative, trimmed segmentation models, discard any segments with less than 3 sites, and create the union of all the modified gbM segments in each gene

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=3

for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  #if (substr(this_sample,10,10)=="_") {this_sample=substr(this_sample,1,9)}
  if (!(this_sample %in% ignore_list)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    #this_model = this_model[substr(this_model$Type,1,4)=="GBM_",]
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    
    if (no_models == 1) {
      overlap_model = this_model_trimmed
    } else {
      overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    }
  }
}

# Add "Chr" to the chromosome names
levels(overlap_model@seqnames@values) = paste0("Chr",as.character(levels(overlap_model@seqnames@values)))
overlap_model@seqinfo@seqnames = levels(overlap_model@seqnames@values)


# Identify accessions with gbM segments collinear with end of gene AT1G01610

for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  #if (substr(this_sample,10,10)=="_") {this_sample=substr(this_sample,1,9)}
  if (!(this_sample %in% ignore_list)) {
    #no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    #this_model = this_model[substr(this_model$Type,1,4)=="GBM_",]
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]

    head(as.data.frame(this_model_trimmed[queryHits(findOverlaps(this_model_trimmed, makeGRangesFromDataFrame(df=cbind("Chromosome"=1, "Locus"=221642), seqnames.field="Chromosome", start.field="Locus", end.field="Locus")))]))
         #if (no_models == 1) {
         #  overlap_model = this_model_trimmed
         #} else {
         #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
         #}
  }
}



# load in all the alternative, trimmed segmentation models, discard any segments with less than 3 sites, and for each CG site count how many modified gbM segments overlap it

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)
min_segment_sites=3

CG_site_gbM_counts = rep(0,length(CG_site_ranges))

for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (!(this_sample %in% ignore_list)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    CG_site_gbM_counts[subjectHits(model_site_olaps)] = CG_site_gbM_counts[subjectHits(model_site_olaps)] + 1
        
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

# Add "Chr" to the chromosome names of CG_site_ranges
levels(CG_site_ranges@seqnames@values) = paste0("Chr",as.character(levels(CG_site_ranges@seqnames@values)))
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)

reference_gff = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff")  # AraPort11 annotation


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

# go through genes identifying first and last locus covered by gbM segments 3 CG sites or longer in more than X accessions
# export the 'protein coding' ones to files for comparison

#for (min_seg_cover in 1:10) {
for (min_seg_cover in 5) {
  gene_firsts = rep(0, length(gene_ranges))
  gene_lasts = rep(0, length(gene_ranges))
  for (this_gene in 1:length(gene_ranges)) {
    cat(paste0("cover: ", min_seg_cover, "  gene: ", this_gene, "\n"))
    gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
    if (length(gene_site_olaps)>0) {
      gene_site_counts = CG_site_gbM_counts[subjectHits(gene_site_olaps)]
      if (sum(gene_site_counts>=min_seg_cover)>0) {
        gene_firsts[this_gene] = CG_site_ranges[min(subjectHits(gene_site_olaps)[gene_site_counts>=min_seg_cover])]@ranges@start
        gene_lasts[this_gene] = CG_site_ranges[max(subjectHits(gene_site_olaps)[gene_site_counts>=min_seg_cover])]@ranges@start + 1
      } else {
        gene_firsts[this_gene] = NA
        gene_lasts[this_gene] = NA
      }
    } else {
      gene_firsts[this_gene] = NA
      gene_lasts[this_gene] = NA
    }
  }
  write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])

# Now identify accessions with Col-0-like gbM levels
relevant_gbM_loci = cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",]
colnames(relevant_gbM_loci)[1]="gene_ID"
relevant_gbM_loci = relevant_gbM_loci[!is.na(relevant_gbM_loci$gene_firsts),]
relevant_gbM_loci.gr = makeGRangesFromDataFrame(relevant_gbM_loci, seqnames.field = "seqnames", start.field = "gene_firsts", end.field = "gene_lasts")
relevant_CG_sites = queryHits(findOverlaps(CG_site_ranges, relevant_gbM_loci.gr))

accession_gbm_levels = rep(0, length(list.files(path=model_dir)))
this_model_no = 0
for (this_model_file in list.files(path=model_dir)) {
  this_model_no = this_model_no + 1
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  names(accession_gbm_levels)[this_model_no] = this_sample
  accession_gbm_levels[this_model_no] = sum(all_samples_meth_status[relevant_CG_sites,this_sample]=="M", na.rm=TRUE) / sum(all_samples_meth_status[relevant_CG_sites,this_sample] %in% c("M","U"), na.rm=TRUE)
}

hist(accession_gbm_levels, breaks = seq(0,1, 0.01))

mean(accession_gbm_levels)
#[1] 0.3137513
sd(accession_gbm_levels)
#[1] 0.04878895

mean(accession_gbm_levels)-sd(accession_gbm_levels)
#[1] 0.2649624
mean(accession_gbm_levels)+sd(accession_gbm_levels)
#[1] 0.3625403

# Col-0 accession levels:
accession_gbm_levels[c("SRX446038", "SRX446039", "SRX446040", "SRX248644")]
#SRX446038 SRX446039 SRX446040 SRX248644 
#0.2388367 0.2970364 0.3153711 0.3873665


# all 4 Col-0 accessions have weirdly differing levels of mCG in gbM. Need to load in Schmitz data as a check

schmitz_all_samples_meth_status = readRDS("../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_all_samples_meth_status.rds")
schmitz_parental_consensus = read.table(file="../../Jay-SMP_variation/5-analysis/SRA035939/parental_consensus.txt")
schmitz_samples_meth_status = cbind(schmitz_all_samples_meth_status, schmitz_parental_consensus)[,c("Chromosome","Locus","V1")]
colnames(schmitz_samples_meth_status)[3]="Schmitz_parental"
schmitz_samples_meth_status$Chromosome = paste0("Chr",substr(schmitz_samples_meth_status$Chromosome,4,4))
merged_samples_meth_status = merge(all_samples_meth_status[,1:2], schmitz_samples_meth_status, by=c("Chromosome","Locus"), all.x=TRUE)
schmitz_relevant_CG_sites = queryHits(findOverlaps(makeGRangesFromDataFrame(df=schmitz_samples_meth_status, seqnames.field = "Chromosome", start.field = "Locus", end.field="Locus"), relevant_gbM_loci.gr))
sum(schmitz_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"]=="M", na.rm=TRUE) / sum(schmitz_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"] %in% c("M","U"), na.rm=TRUE)
#[1] 0.2901187

# check out the values for the Schmitz GBM segments
segmentation_model.gr = readRDS(file = "../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_segmentation_model_draft3.rds")
GBM_segments = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]
relevant_CG_sites = queryHits(findOverlaps(CG_site_ranges, GBM_segments))
schmitz_relevant_CG_sites = queryHits(findOverlaps(makeGRangesFromDataFrame(df=schmitz_samples_meth_status, seqnames.field = "Chromosome", start.field = "Locus", end.field="Locus"), GBM_segments))

# Export longer gbM Schmitz segments
schmitz_relevant_gbm_segments = table(subjectHits(findOverlaps(CG_site_ranges, GBM_segments)))
schmitz_relevant_gbm_segments = schmitz_relevant_gbm_segments[schmitz_relevant_gbm_segments>3]

# Amy wants a gene ID for each segment, so do the overlapping:
relevant_genes = gff.genes[findOverlaps(GBM_segments, gene_ranges, select="first"),"gene_ID"]
# Which gene is best for each segment? We take the first one, because we expect that this is at least overlapping at 5', so in most cases, will hopefully be the longest segment

#write.table(as.data.frame(GBM_segments)[names(schmitz_relevant_gbm_segments)[schmitz_relevant_gbm_segments>3],c("seqnames","start","end")], file="Schmitz_gbM_segments_gt_3_CG_sites.tsv", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
write.table(cbind("gene_ID"=relevant_genes, as.data.frame(GBM_segments))[names(schmitz_relevant_gbm_segments)[schmitz_relevant_gbm_segments>3],c("gene_ID","seqnames","start","end")], file="Schmitz_gbM_segments_gt_3_CG_sites.tsv", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
# NB This produces duplicate entries for overlapping genes. As a fix, the first of each of such entries was chosen, in Excel, before the file was passed to Amy for modelling


accession_gbm_levels = rep(0, length(list.files(path=model_dir)))
this_model_no = 0
for (this_model_file in list.files(path=model_dir)) {
  this_model_no = this_model_no + 1
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  names(accession_gbm_levels)[this_model_no] = this_sample
  accession_gbm_levels[this_model_no] = sum(all_samples_meth_status[relevant_CG_sites,this_sample]=="M", na.rm=TRUE) / sum(all_samples_meth_status[relevant_CG_sites,this_sample] %in% c("M","U"), na.rm=TRUE)
}

hist(accession_gbm_levels, breaks = seq(0,1, 0.01))

mean(accession_gbm_levels)
#[1] 0.5772374
sd(accession_gbm_levels)
#[1] 0.06789811

mean(accession_gbm_levels)-sd(accession_gbm_levels)
#[1] 0.5093393
mean(accession_gbm_levels)+sd(accession_gbm_levels)
#[1] 0.6451355

# Col-0 accession levels:
accession_gbm_levels[c("SRX446038", "SRX446039", "SRX446040", "SRX248644")]
#SRX446038 SRX446039 SRX446040 SRX248644 
#0.5794306 0.6366895 0.6524200 0.7174731

sum(schmitz_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"]=="M", na.rm=TRUE) / sum(merged_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"] %in% c("M","U"), na.rm=TRUE)
#[1] 0.6124375

# Let's export a spreadsheet of metadata combined with accession_gbm_levels
# read sample info in
raw_data_dir = "../0-raw_data/"
sample_list = read.table(paste0(raw_data_dir,"sample_info_synthesis.txt"), sep="\t", header=TRUE)
sample_list$Name = str_split_fixed(sample_list$Title, "\\(",2)[,1]
#sample_list1$Ecotype_id = as.numeric(str_split_fixed(str_split_fixed(sample_list1$Title, "\\(",2)[,2], "\\)", 2)[,1])
sample_list$Ecotype_id = as.numeric(str_split_fixed(sample_list$Sample, "_", 2)[,1])

accession_gbm_levels = as.data.frame(accession_gbm_levels)
accession_gbm_levels$SRA.Accession = row.names(accession_gbm_levels)
sample_list = merge(sample_list, accession_gbm_levels, by="SRA.Accession", all=TRUE)

# add in proportion of sites well-covered as a proxy for coverage
accession_coverage_ratio = colSums((all_samples_meth_status[,4:ncol(all_samples_meth_status)]=="M") | (all_samples_meth_status[,4:ncol(all_samples_meth_status)]=="U"), na.rm=TRUE)
accession_coverage_ratio = accession_coverage_ratio / nrow(all_samples_meth_status)
accession_coverage_ratio = as.data.frame(accession_coverage_ratio)
accession_coverage_ratio$SRA.Accession = row.names(accession_coverage_ratio)
sample_list = merge(sample_list, accession_coverage_ratio, by="SRA.Accession", all=TRUE)

write.table(sample_list, file="m1001_samples_mCG_density_coverage.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#How many accessions have methylation within 1SD of the mean?
#accession_gbm_levels is now a data frame, let's turn it back into a vector
accession_gbm_level_names = accession_gbm_levels$SRA.Accession
accession_gbm_levels=accession_gbm_levels$accession_gbm_levels
names(accession_gbm_levels)=accession_gbm_level_names
sum(accession_gbm_levels>=mean(accession_gbm_levels)-sd(accession_gbm_levels) & accession_gbm_levels<=mean(accession_gbm_levels)+sd(accession_gbm_levels))  
#sum(accession_gbm_levels$accession_gbm_levels>=mean(accession_gbm_levels$accession_gbm_levels)-sd(accession_gbm_levels$accession_gbm_levels) & accession_gbm_levels$accession_gbm_levels<=mean(accession_gbm_levels$accession_gbm_levels)+sd(accession_gbm_levels$accession_gbm_levels))  
#[1] 891  
  
#These 891 accessions are the 'Col-0-like' ones, or at least, Col-0 and these all lie within 1SD of the mean of the methylation level of the Col-0 gbM segments
#So we are interested in genes which contain gbM segments in at least 20% of this group

relevant_samples = names(accession_gbm_levels[accession_gbm_levels>=mean(accession_gbm_levels)-sd(accession_gbm_levels) & accession_gbm_levels<=mean(accession_gbm_levels)+sd(accession_gbm_levels)])

gene_gbm_counts = rep(0, length(gene_ranges))
no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (this_sample %in% relevant_samples) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model$Chromosome = paste0("Chr",this_model$Chromosome)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    model_gene_olaps = findOverlaps(this_model_trimmed, gene_ranges)
    gene_gbm_counts[subjectHits(model_gene_olaps)] = gene_gbm_counts[subjectHits(model_gene_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

relevant_gbM_loci = (cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts, gene_gbm_counts)[substr(gff.genes$gene_name,1,2)=="AT",])
write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts, gene_gbm_counts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,"_gbm_counts.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# limit it to genes with gbM segments in at least 20% of Col-0-like accessions
relevant_gbM_loci = relevant_gbM_loci[(!is.na(relevant_gbM_loci$gene_firsts)) & (relevant_gbM_loci$gene_gbm_counts>891/5),]
relevant_gbM_loci.gr = makeGRangesFromDataFrame(relevant_gbM_loci, seqnames.field = "seqnames", start.field = "gene_firsts", end.field = "gene_lasts")


# The above method calculates relevant gbM_loci by taking first site to last site which are covered by gbM segments of at least 3 sites in at least 5 accessions.
# We now want to only include segments which meet that criteria, rather than potentially unmethylated segments between two methylated ones

# We traverse CG_site_gbM_counts one gene at a time (which already contains counts for each site of overlapping gbM segments longer than 3 sites) and identify contiguous segments within each gene which are covered by relevant gbM segments in at least 5 accessions

#for (min_seg_cover in 1:10) {
for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  gene_gbM_segments = data.frame()
  for (this_gene in 1:length(gene_ranges)) {
    cat(paste0("cover: ", min_seg_cover, "  gene: ", this_gene, "\n"))
    doing_one = FALSE
    this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
    this_start = 0
    this_end = 0
    gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
    gene_site_counts = CG_site_gbM_counts[subjectHits(gene_site_olaps)]
    no_sites = length(gene_site_olaps)
    if (no_sites>0) {
      prev_site = 0
      last_site = FALSE
      save_this = FALSE
      for (this_site in 1:no_sites) {
        if (this_site==no_sites) {}
        if (doing_one) {
          if (gene_site_counts[this_site]>=min_seg_cover) {
            if (this_site==no_sites) {
              this_end = this_site
              save_this = TRUE
            }
          } else {
            this_end = this_site -1
            doing_one = FALSE
            save_this = TRUE
          }
        } else {
          if (gene_site_counts[this_site]>=min_seg_cover) {
            doing_one=TRUE
            this_start = this_site
            if (this_site==no_sites) {
              this_end = this_site
              save_this = TRUE
            }
          }
        }
        if (save_this) {
          gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start, "end"=CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1))
          colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
          save_this=FALSE
        }
      }
    } else {
      #gene_firsts[this_gene] = NA
      #gene_lasts[this_gene] = NA
    }
  }
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])



# partition the sample list to find some for modelling
# no need to bother with anything below accession_coverage_ratio of 0.6 - that still allows us to include the most ad least methylated ones (apart from the three lowly methylated accessions with coverage<10%)
relevant_samples = sample_list[sample_list$accession_coverage_ratio>0.6,]

# Pick the 10 most and least methylated accessions
# 10 least methylated accessions
least_methylated = head(relevant_samples[order(relevant_samples$accession_gbm_levels),], 10)
# 10 most methylated
most_methylated = head(relevant_samples[order(-relevant_samples$accession_gbm_levels),], 10)

# remainder in the middle:
middle_methylated = relevant_samples[(relevant_samples$accession_gbm_levels>max(least_methylated$accession_gbm_levels)) & (relevant_samples$accession_gbm_levels<min(most_methylated$accession_gbm_levels)),]

#0.2466-0.4159
#0.4197-0.4505

#0.7081-0.7314
#0.7358-0.7932

#deciles for the in-between accessions:
#1	0.4159	0.44745
#2	0.44745	0.479
#3	0.479	0.51055
#4	0.51055	0.5421
#5	0.5421	0.57365
#6	0.57365	0.6052
#7	0.6052	0.63675
#8	0.63675	0.6683
#9	0.6683	0.69985
#10	0.69985	0.7314


# pick the accessions for export to Amy for modelling
# we are going to choose 10 sets of 10 accessions each, each set of 10 to represent a segment from the distribution of mCG density in Col-0 gbM segments

site_segment_olaps = findOverlaps(CG_site_ranges, relevant_gbM_loci.gr)
relevant_CG_sites = queryHits(site_segment_olaps)
relevant_gene_IDs = relevant_gbM_loci[subjectHits(site_segment_olaps),1]

middle_start = min(middle_methylated$accession_gbm_levels)
middle_end = max(middle_methylated$accession_gbm_levels)

#write.table(cbind("gene_ID"=relevant_gene_IDs, all_samples_meth_status[relevant_CG_sites, c("Chromosome", "Locus", least_methylated$SRA.Accession)]), file=paste0("m1001_samples_meth_status_mCG_density_low.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(cbind("gene_ID"=relevant_gene_IDs, all_samples_meth_status[relevant_CG_sites, c("Chromosome", "Locus", most_methylated$SRA.Accession)]), file=paste0("m1001_samples_meth_status_mCG_density_high.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

modelling_accessions=list()
for (this_section in 1:10) {
  section_start = middle_start + (this_section-1)*(middle_end-middle_start)/10
  section_end = middle_start + (this_section)*(middle_end-middle_start)/10
  modelling_accessions[[this_section]] = head(middle_methylated[order(-middle_methylated$accession_coverage_ratio),][(middle_methylated$accession_gbm_levels>=section_start) & (middle_methylated$accession_gbm_levels<section_end),"SRA.Accession"],10)
  cat(paste0(modelling_accessions[[this_section]]), "\n")  
}

modelling_accessions[[11]] = least_methylated$SRA.Accession
modelling_accessions[[12]] = most_methylated$SRA.Accession

for (this_section in 1:12) {
  #Need to make this gene, chromosome, start, end, status
  for (this_accession in modelling_accessions[[this_section]]) {
    write.table(cbind("gene_ID"=relevant_gene_IDs, "Chromosome"=all_samples_meth_status[relevant_CG_sites,"Chromosome"], "Start"=all_samples_meth_status[relevant_CG_sites,"Locus"], "End"=all_samples_meth_status[relevant_CG_sites,"Locus"]+1, "Status"=all_samples_meth_status[relevant_CG_sites,this_accession]), file=paste0("m1001_samples_meth_status_mCG_density_",this_section,"_",this_accession,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}

# Initial extremities for modelling:

#Least methylated (both Ecker):
#SRX248646 Cvi_0
#SRX2190740 UKID116 (5822)-2

#Most methylated (Dubin) and very methylated (Ecker):
#SRX445897 5856_10C
#SRX1664724 Gron 12 (9386)


# CG_site_gbM_counts contains counts for all accessions. We want it just for Col0-like ones
CG_site_gbM_counts_Col0like = rep(0,length(CG_site_ranges))

# at this point, check that CG_site_ranges have chromosomes identified just as numbers rather than ChrX
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)

for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if ((this_sample %in% relevant_samples)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    CG_site_gbM_counts_Col0like[subjectHits(model_site_olaps)] = CG_site_gbM_counts_Col0like[subjectHits(model_site_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}



# identify protein_coding genes
#gff.genes$V9 contains "locus_type=protein_coding"

gene_locus_type = rep("", nrow(gff.genes))
for (this_gene in 1:nrow(gff.genes)) {
  gene_locus_type_info = strsplit(gff.genes[this_gene, "V9"], ";")
  this_locus_type = gene_locus_type_info[[1]][substr(gene_locus_type_info[[1]],1,11)=="locus_type="]
  
  if (length(this_locus_type)>0) {
    this_locus_type = substr(this_locus_type,12,str_length(this_locus_type))
    gene_locus_type[this_gene] = this_locus_type
  }  
  cat(paste0(gff.genes[this_gene,"gene_ID"], " ", this_locus_type, "\n"))
}

write.table(cbind(gff.genes, gene_locus_type), file="gene_types.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# We traverse CG_site_gbM_counts_Col0like one gene at a time (which already contains counts for each site of overlapping gbM segments longer than 3 sites) and identify contiguous segments within each gene which are covered by relevant gbM segments in at least x% of accessions
# Make sure gene_ranges has numbers rather than ChrX

for (min_quantile_cov in 1:9) {
  min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  gene_gbM_segments = data.frame()
  for (this_gene in 1:length(gene_ranges)) {
    if(gene_locus_type[this_gene]=="protein_coding") {
      cat(paste0("cover: ", min_seg_cover, "  gene: ", this_gene, "\n"))
      doing_one = FALSE
      this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
      this_start = 0
      this_end = 0
      gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
      gene_site_counts = CG_site_gbM_counts_Col0like[subjectHits(gene_site_olaps)]
      no_sites = length(gene_site_olaps)
      if (no_sites>0) {
        prev_site = 0
        last_site = FALSE
        save_this = FALSE
        for (this_site in 1:no_sites) {
          if (this_site==no_sites) {}
          if (doing_one) {
            if (gene_site_counts[this_site]>=min_seg_cover) {
              if (this_site==no_sites) {
                this_end = this_site
                save_this = TRUE
              }
            } else {
              this_end = this_site -1
              doing_one = FALSE
              save_this = TRUE
            }
          } else {
            if (gene_site_counts[this_site]>=min_seg_cover) {
              doing_one=TRUE
              this_start = this_site
              if (this_site==no_sites) {
                this_end = this_site
                save_this = TRUE
              }
            }
          }
          if (save_this) {
            gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start, "end"=CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1))
            colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
            save_this=FALSE
          }
        }
      } else {
        #gene_firsts[this_gene] = NA
      #gene_lasts[this_gene] = NA
      }
    }
  }
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])


pdf("annotation_summary_stats.pdf")
# For each annotation - plot CG density vs. mCG in Col-0, average number of accessions each CG site is annotated as gbM
for (min_quantile_cov in 1:9) {
  min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  #gene_gbM_segments = data.frame()
  
  segment_site_stats = data.frame()
  
  gene_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
  
  for (this_segment in 1:nrow(gene_gbM_segments)) {
    this_segment.gr = makeGRangesFromDataFrame(gene_gbM_segments[this_segment,], seqnames.field="chromosome", start.field="start", end.field="end")
    segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    these_sites = queryHits(segment_site_olaps)
    gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    if (no_sites>0) {
      segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      colnames(segment_site_stats) = c("gene_ID", "chromosome", "start", "end", "no_sites", "mean_coverage")
    }
  }
  segment_site_stats$start = as.numeric(segment_site_stats$start)
  segment_site_stats$end = as.numeric(segment_site_stats$end)
  segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
  segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
  
  print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/segment_site_stats$no_sites, y=segment_site_stats$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("R=", cor((1+segment_site_stats$end-segment_site_stats$start)/segment_site_stats$no_sites, segment_site_stats$mean_coverage))))
  
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}
dev.off()





# The stats plots show there are a bunch of very CG-site-dense segments which are less methylated than the 'other' 'normal' segments
# Going to try filtering for 5% minimum number of accessions generally, but then removing any segments with high density (mean CG site spacing below 20ish) with less than 20% or 30% of accessions
# Make sure gene_ranges has numbers rather than ChrX

#for (min_quantile_cov in 1:9) {
min_quantile_cov=0.5
for (min_dense_quantile_cov in 2:3) {
  for (min_site_spacing in 15:25) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    gene_gbM_segments = data.frame()
    for (this_gene in 1:length(gene_ranges)) {
      if(gene_locus_type[this_gene]=="protein_coding") {
        cat(paste0("cover: ", min_seg_cover, "spacing: ", min_site_spacing, "  gene: ", this_gene, "\n"))
        doing_one = FALSE
        this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
        this_start = 0
        this_end = 0
        gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
        gene_site_counts = CG_site_gbM_counts_Col0like[subjectHits(gene_site_olaps)]
        no_sites = length(gene_site_olaps)
        if (no_sites>0) {
          prev_site = 0
          last_site = FALSE
          save_this = FALSE
          for (this_site in 1:no_sites) {
            if (this_site==no_sites) {}
            if (doing_one) {
              if (gene_site_counts[this_site]>=min_seg_cover) {
                if (this_site==no_sites) {
                  this_end = this_site
                  save_this = TRUE
                }
              } else {
                this_end = this_site -1
                doing_one = FALSE
                save_this = TRUE
              }
            } else {
              if (gene_site_counts[this_site]>=min_seg_cover) {
                doing_one=TRUE
                this_start = this_site
                if (this_site==no_sites) {
                  this_end = this_site
                  save_this = TRUE
                }
              }
            }
            if (save_this) {
              segment_site_spacing = (1+CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1-CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start)/(no_sites-1)
              if ((segment_site_spacing>=min_site_spacing) | (mean(gene_site_counts[this_start:this_end])>=min_dense_seg_cover)) {
                gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start, "end"=CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1))
                colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
              }
              save_this=FALSE
            }
          }
        } else {
          #gene_firsts[this_gene] = NA
          #gene_lasts[this_gene] = NA
        }
      }
    }
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,"_min_dense_",min_dense_quantile_cov,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}

View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])


pdf("annotation_summary_stats_2.pdf")
# For each annotation - plot CG density vs. mCG in Col-0, average number of accessions each CG site is annotated as gbM
min_quantile_cov=0.5
for (min_dense_quantile_cov in 2:3) {
  for (min_site_spacing in 15:25) {
    #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    segment_site_stats = data.frame()
    
    gene_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,"_min_dense_",min_dense_quantile_cov,".tsv"), sep="\t", header=TRUE)
    
    for (this_segment in 1:nrow(gene_gbM_segments)) {
      this_segment.gr = makeGRangesFromDataFrame(gene_gbM_segments[this_segment,], seqnames.field="chromosome", start.field="start", end.field="end")
      segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
      these_sites = queryHits(segment_site_olaps)
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
      if (no_sites>0) {
        segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats) = c("gene_ID", "chromosome", "start", "end", "no_sites", "mean_coverage")
      }
    }
    segment_site_stats$start = as.numeric(segment_site_stats$start)
    segment_site_stats$end = as.numeric(segment_site_stats$end)
    segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
    segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
    
    print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/segment_site_stats$no_sites, y=segment_site_stats$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,"_min_dense_",min_dense_quantile_cov,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }}
dev.off()




# The stats plots show there are a bunch of very CG-site-dense segments which are less methylated than the 'other' 'normal' segments
# Going to try filtering for 5% minimum number of accessions generally, but then removing any segments with high density (mean CG site spacing below 20ish) regardless of how frequently they are methylated among accessions
# In the end 5% didn't do it, and the min_site_spacings we were trying were too soft. Went down to 0.5%, and upped spacings - 40 is good
# Make sure gene_ranges has numbers rather than ChrX

#for (min_quantile_cov in 1:9) {
for (min_quantile_cov in c(0.05, 0.1, 0.2, 0.5)) {
#for (min_dense_quantile_cov in 2:3) {
  for (min_site_spacing in c(20,30,40)) {
  #for (min_site_spacing in 15:25) {
      min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    gene_gbM_segments = data.frame()
    for (this_gene in 1:length(gene_ranges)) {
      if(gene_locus_type[this_gene]=="protein_coding") {
        cat(paste0("cover: ", min_seg_cover, "spacing: ", min_site_spacing, "  gene: ", this_gene, "\n"))
        doing_one = FALSE
        this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
        this_start = 0
        this_end = 0
        gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
        gene_site_counts = CG_site_gbM_counts_Col0like[subjectHits(gene_site_olaps)]
        no_sites = length(gene_site_olaps)
        if (no_sites>0) {
          prev_site = 0
          last_site = FALSE
          save_this = FALSE
          for (this_site in 1:no_sites) {
            if (this_site==no_sites) {}
            if (doing_one) {
              if (gene_site_counts[this_site]>=min_seg_cover) {
                if (this_site==no_sites) {
                  this_end = this_site
                  save_this = TRUE
                }
              } else {
                this_end = this_site -1
                doing_one = FALSE
                save_this = TRUE
              }
            } else {
              if (gene_site_counts[this_site]>=min_seg_cover) {
                doing_one=TRUE
                this_start = this_site
                if (this_site==no_sites) {
                  this_end = this_site
                  save_this = TRUE
                }
              }
            }
            if (save_this) {
              segment_site_spacing = (1+CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1-CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start)/(no_sites-1)
              #if ((segment_site_spacing>=min_site_spacing) | (mean(gene_site_counts[this_start:this_end])>=min_dense_seg_cover)) {
              if ((segment_site_spacing>=min_site_spacing)) {
                  gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start, "end"=CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1))
                colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
              }
              save_this=FALSE
            }
          }
        } else {
          #gene_firsts[this_gene] = NA
          #gene_lasts[this_gene] = NA
        }
      }
    }
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}


View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])


pdf("annotation_summary_stats_3.pdf")
# For each annotation - plot CG density vs. mCG in Col-0, average number of accessions each CG site is annotated as gbM
min_quantile_cov=0.5
#for (min_quantile_cov in 1:9) {
for (min_quantile_cov in c(0.05, 0.1, 0.2, 0.5)) {
  #for (min_dense_quantile_cov in 2:3) {
  for (min_site_spacing in c(20,30,40)) {
    #for (min_site_spacing in 15:25) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    segment_site_stats = data.frame()
    
    gene_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,".tsv"), sep="\t", header=TRUE)
    
    for (this_segment in 1:nrow(gene_gbM_segments)) {
      this_segment.gr = makeGRangesFromDataFrame(gene_gbM_segments[this_segment,], seqnames.field="chromosome", start.field="start", end.field="end")
      segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
      these_sites = queryHits(segment_site_olaps)
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
      if (no_sites>0) {
        segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats) = c("gene_ID", "chromosome", "start", "end", "no_sites", "mean_coverage")
      }
    }
    segment_site_stats$start = as.numeric(segment_site_stats$start)
    segment_site_stats$end = as.numeric(segment_site_stats$end)
    segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
    segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
    
    print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/segment_site_stats$no_sites, y=segment_site_stats$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}
dev.off()


# Generate the 5% threshold set of segments (also 2% and 1%). Find overlap with Col-0 segments. Plot no accessions vs. density for Col-0 segments and non-Col-0 segments

for (min_quantile_cov in c(0.05, 0.1, 0.2, 0.5)) {
for (min_quantile_cov in c(0.05)) {
#for (min_dense_quantile_cov in 2:3) {
#for (min_site_spacing in 15:25) {
  min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
  #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  gene_gbM_segments = data.frame()
  for (this_gene in 1:length(gene_ranges)) {
    if(gene_locus_type[this_gene]=="protein_coding") {
      cat(paste0("cover: ", min_seg_cover, "  gene: ", this_gene, "\n"))
      doing_one = FALSE
      this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
      this_start = 0
      this_end = 0
      gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
      gene_site_counts = CG_site_gbM_counts_Col0like[subjectHits(gene_site_olaps)]
      no_sites = length(gene_site_olaps)
      if (no_sites>0) {
        prev_site = 0
        last_site = FALSE
        save_this = FALSE
        for (this_site in 1:no_sites) {
          if (this_site==no_sites) {}
          if (doing_one) {
            if (gene_site_counts[this_site]>=min_seg_cover) {
              if (this_site==no_sites) {
                this_end = this_site
                save_this = TRUE
              }
            } else {
              this_end = this_site -1
              doing_one = FALSE
              save_this = TRUE
            }
          } else {
            if (gene_site_counts[this_site]>=min_seg_cover) {
              doing_one=TRUE
              this_start = this_site
              if (this_site==no_sites) {
                this_end = this_site
                save_this = TRUE
              }
            }
          }
          if (save_this) {
            segment_site_spacing = (1+CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1-CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start)/(no_sites-1)
            #if ((segment_site_spacing>=min_site_spacing) | (mean(gene_site_counts[this_start:this_end])>=min_dense_seg_cover)) {
            #if ((segment_site_spacing>=min_site_spacing)) {
              gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start, "end"=CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1))
              colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
            #}
            save_this=FALSE
          }
        }
      } else {
        #gene_firsts[this_gene] = NA
        #gene_lasts[this_gene] = NA
      }
    }
  }
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}
#}


View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])

# Fix chromosome names in GBM_segments
levels(GBM_segments@seqnames@values) = substr(as.character(levels(GBM_segments@seqnames@values)),4,4)
GBM_segments@seqinfo@seqnames = levels(GBM_segments@seqnames@values)


pdf("annotation_summary_stats_4.pdf")
# For each annotation - plot CG density vs. mCG in Col-0, average number of accessions each CG site is annotated as gbM
for (min_quantile_cov in c(0.5, 0.2, 0.1, 0.05)) {
  #min_quantile_cov=0.05
#for (min_dense_quantile_cov in 2:3) {
#for (min_site_spacing in 15:25) {
#for (min_quantile_cov in 1:9) {
  min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
  #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
  
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  #gene_gbM_segments = data.frame()
  
  
  these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
  these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  # do the bit not present in the Schmitz annotation
  gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)

  segment_site_stats = data.frame()
  
  for (this_segment in 1:length(gene_gbM_segments)) {
    this_segment.gr = gene_gbM_segments[this_segment]
    segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    these_sites = queryHits(segment_site_olaps)
    gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    if (no_sites>0) {
      #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      segment_site_stats = rbind(segment_site_stats, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      colnames(segment_site_stats) = c("start", "end", "no_sites", "mean_coverage")
    }
  }
  segment_site_stats$start = as.numeric(segment_site_stats$start)
  segment_site_stats$end = as.numeric(segment_site_stats$end)
  segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
  segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
  
  #print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), y=segment_site_stats$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))))
  print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), y=segment_site_stats$mean_coverage, colour=log(1+segment_site_stats$end-segment_site_stats$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))) + xlim(0,100) +ylim(0,900))
  
  # do the bit which is present in the Schmitz annotation
  gene_gbM_segments = setdiff(these_gbM_segments.gr, setdiff(these_gbM_segments.gr, GBM_segments))
  
  segment_site_stats2 = data.frame()
  
  for (this_segment in 1:length(gene_gbM_segments)) {
    this_segment.gr = gene_gbM_segments[this_segment]
    segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    these_sites = queryHits(segment_site_olaps)
    gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    if (no_sites>0) {
      #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    }
  }
  segment_site_stats2$start = as.numeric(segment_site_stats2$start)
  segment_site_stats2$end = as.numeric(segment_site_stats2$end)
  segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
  segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
  
  #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
  print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
  
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}
#}
dev.off()



# For each annotation - keep the segments overlapping with Col-0, merge with the segments not overlapping but with CG density below threshold, and plot CG density vs. mCG in Col-0, average number of accessions each CG site is annotated as gbM
for (min_quantile_cov in c(0.2, 0.5, 0.1, 0.05)) {
  #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 40, 30)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
      #for (min_quantile_cov in 1:9) {
  min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
  #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
  
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  #gene_gbM_segments = data.frame()
  
  
  these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
  these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  # do the bit not present in the Schmitz annotation
  gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)
  gene_gbM_segments = gene_gbM_segments[unique(queryHits(findOverlaps(gene_gbM_segments, CG_site_ranges)))]
  
  segment_site_stats = data.frame()

  segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
  for (this_segment in 1:length(gene_gbM_segments)) {
    this_segment.gr = gene_gbM_segments[this_segment]
    these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
    gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
    if (length(these_sites)>0) {
      #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      segment_site_stats = rbind(segment_site_stats, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      colnames(segment_site_stats) = c("start", "end", "no_sites", "mean_coverage")
    }
  }
  
  #for (this_segment in 1:length(gene_gbM_segments)) {
  #  this_segment.gr = gene_gbM_segments[this_segment]
  #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
  #  these_sites = queryHits(segment_site_olaps)
  #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
  #  if (length(these_sites)>0) {
  #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
  #    segment_site_stats = rbind(segment_site_stats, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
  #    colnames(segment_site_stats) = c("start", "end", "no_sites", "mean_coverage")
  #  }
  #}
  segment_site_stats$start = as.numeric(segment_site_stats$start)
  segment_site_stats$end = as.numeric(segment_site_stats$end)
  segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
  segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
  
  # Pick the relevant non-Col-0 segments to append to the Col-0 segments
  relevant_non_Col0_segments = gene_gbM_segments[(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1)>=min_site_spacing]
  # Write out the non-Col-0 segments
  write.table(relevant_non_Col0_segments, file=paste0("relevant_non_Col0_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  
  #print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), y=segment_site_stats$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))))
  #print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), y=segment_site_stats$mean_coverage, colour=log(1+segment_site_stats$end-segment_site_stats$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))) + xlim(0,100) +ylim(0,900))
  
  # do the bit which is present in the Schmitz annotation
  gene_gbM_segments = setdiff(these_gbM_segments.gr, setdiff(these_gbM_segments.gr, GBM_segments))
  # Write out the Col-0 segments
  write.table(gene_gbM_segments, file=paste0("relevant_Col0_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  
  # Merge in the relevant non-Col-0 segments
  gene_gbM_segments = reduce(union(gene_gbM_segments, relevant_non_Col0_segments))
  
  # Write out the merged model
  write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  
  segment_site_stats2 = data.frame()
  
  segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
  for (this_segment in 1:length(gene_gbM_segments)) {
    this_segment.gr = gene_gbM_segments[this_segment]
    these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
    gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
    if (length(these_sites)>0) {
      #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    }
  }
  #for (this_segment in 1:length(gene_gbM_segments)) {
  #  this_segment.gr = gene_gbM_segments[this_segment]
  #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
  #  these_sites = queryHits(segment_site_olaps)
  #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
  #  if (length(these_sites)>0) {
  #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
  #    segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
  #    colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
  #  }
  #}
  segment_site_stats2$start = as.numeric(segment_site_stats2$start)
  segment_site_stats2$end = as.numeric(segment_site_stats2$end)
  segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
  segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
  
  #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
  #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
  print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900) + theme(legend.title = element_blank()))
  
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(segment_site_stats2, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}
}


# For each annotation - take the ones prepared above, and filter out any segments with CG site density above threshold
for (min_quantile_cov in c(0.2, 0.5, 0.1, 0.05)) {
  #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 40, 30)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
    #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    
    these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", header=TRUE)
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="seqnames", start.field="start", end.field="end")
    # do the bit not present in the Schmitz annotation
    #gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)
    gene_gbM_segments = these_gbM_segments.gr
    gene_gbM_segments = gene_gbM_segments[unique(queryHits(findOverlaps(gene_gbM_segments, CG_site_ranges)))]
    
    segment_site_stats = data.frame()
    
    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
    for (this_segment in 1:length(gene_gbM_segments)) {
      this_segment.gr = gene_gbM_segments[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats = rbind(segment_site_stats, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats = rbind(segment_site_stats, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats$start = as.numeric(segment_site_stats$start)
    segment_site_stats$end = as.numeric(segment_site_stats$end)
    segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
    segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
    
    # Pick the relevant non-Col-0 segments to append to the Col-0 segments
    relevant_Col0_or_not_segments = gene_gbM_segments[(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1)>=min_site_spacing]
    # Write out the non-Col-0 segments
    write.table(relevant_Col0_or_not_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
    #print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), y=segment_site_stats$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))))
    #print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), y=segment_site_stats$mean_coverage, colour=log(1+segment_site_stats$end-segment_site_stats$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1), segment_site_stats$mean_coverage))) + xlim(0,100) +ylim(0,900))
    

    gene_gbM_segments = relevant_Col0_or_not_segments
    
    segment_site_stats2 = data.frame()
    
    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
    for (this_segment in 1:length(gene_gbM_segments)) {
      this_segment.gr = gene_gbM_segments[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats2$start = as.numeric(segment_site_stats2$start)
    segment_site_stats2$end = as.numeric(segment_site_stats2$end)
    segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
    segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
    
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
    print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900) + theme(legend.title = element_blank()))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats2, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}

# Fix chromosome names in gene_ranges
levels(gene_ranges@seqnames@values) = substr(as.character(levels(gene_ranges@seqnames@values)),4,4)
gene_ranges@seqinfo@seqnames = levels(gene_ranges@seqnames@values)

# Load in the segments, and add gene IDs
for (min_quantile_cov in c(0.2, 0.1, 0.5, 0.05)) {
    #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 30, 40)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
      #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()

    #these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", header=TRUE)
    these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_spacing_",min_site_spacing,".tsv"), sep="\t", header=TRUE)
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="seqnames", start.field="start", end.field="end")

    # do the bit not present in the Schmitz annotation
    gene_gbM_segments = these_gbM_segments.gr

    gene_gbM_segments = cbind(gff.genes[findOverlaps(gene_gbM_segments, gene_ranges, select="first"),c("gene_ID","V1","V4","V5")], gene_gbM_segments@ranges@start, gene_gbM_segments@ranges@start+gene_gbM_segments@ranges@width-1)
    
    colnames(gene_gbM_segments) = c("gene_ID","Chromosome","gene_start","gene_end","start","end")

    #write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,"_2.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_spacing_",min_site_spacing,"_2.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }    
}

pdf("annotation_summary_stats_5.pdf")
# For each annotation - keep the segments overlapping with Col-0, merge with the segments not overlapping but with CG density below threshold, and plot CG density vs. mCG in Col-0, average number of accessions each CG site is annotated as gbM
for (min_quantile_cov in c(0.2, 0.5, 0.1, 0.05)) {
  #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 30, 40)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
    #for (min_site_spacing in c(20, 30, 40)) {
    #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    paste0("relevant_Col0_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv")
    
    #these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", header=TRUE)
    #these_gbM_segments = read.table(file=paste0("relevant_Col0_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", header=TRUE)
    these_gbM_segments = read.table(file=paste0("relevant_non_Col0_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", header=TRUE)
    
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="seqnames", start.field="start", end.field="end")
    # do the bit not present in the Schmitz annotation
    gene_gbM_segments = these_gbM_segments.gr

    segment_site_stats2 = data.frame()

    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
    for (this_segment in 1:length(gene_gbM_segments)) {
      this_segment.gr = gene_gbM_segments[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats2$start = as.numeric(segment_site_stats2$start)
    segment_site_stats2$end = as.numeric(segment_site_stats2$end)
    segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
    segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
    
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
    print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900) + theme(legend.title = element_blank()))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(segment_site_stats2, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}
dev.off()



# For each annotation - take the ones prepared above, and filter out any genes with Col-0 segments with CG site density above threshold
for (min_quantile_cov in c(0.2, 0.5, 0.1, 0.05)) {
  #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 40, 30)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
      #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    # do the bit not present in the Schmitz annotation
    #gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)
    #gene_gbM_segments = gene_gbM_segments[unique(queryHits(findOverlaps(gene_gbM_segments, CG_site_ranges)))]
    
    # do the bit which is present in the Schmitz annotation
    gene_gbM_segments = setdiff(these_gbM_segments.gr, setdiff(these_gbM_segments.gr, GBM_segments))
    
    segment_site_stats = data.frame()
    
    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
    for (this_segment in 1:length(gene_gbM_segments)) {
      this_segment.gr = gene_gbM_segments[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats = rbind(segment_site_stats, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats = rbind(segment_site_stats, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats$start = as.numeric(segment_site_stats$start)
    segment_site_stats$end = as.numeric(segment_site_stats$end)
    segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
    segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
    
    # Pick the relevant Col-0 segments 
    relevant_Col0_segments = gene_gbM_segments[(1+segment_site_stats$end-segment_site_stats$start)/(segment_site_stats$no_sites-1)>=min_site_spacing]
    
    # Find the genes corresponding to these
    relevant_Col0_segments = cbind(gff.genes[findOverlaps(relevant_Col0_segments, gene_ranges, select="first"),c("gene_ID","V1","V4","V5")], relevant_Col0_segments@ranges@start, relevant_Col0_segments@ranges@start+relevant_Col0_segments@ranges@width-1)
    
    colnames(relevant_Col0_segments) = c("gene_ID","Chromosome","gene_start","gene_end","start","end")
    
    genes_to_keep = unique(relevant_Col0_segments$gene_ID)
    
    these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_spacing_",min_site_spacing,"_2.tsv"), sep="\t", header=TRUE)

    these_gbM_segments = these_gbM_segments[these_gbM_segments$gene_ID %in% genes_to_keep,]
    write.table(these_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_strict_spacing_",min_site_spacing,"_2.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="Chromosome", start.field="start", end.field="end")
    # do the bit not present in the Schmitz annotation
    #gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)
    gene_gbM_segments = these_gbM_segments.gr

    # Fix chromosome names in gene_gbM_segments
    levels(gene_gbM_segments@seqnames@values) = substr(as.character(levels(gene_gbM_segments@seqnames@values)),4,4)
    gene_gbM_segments@seqinfo@seqnames = levels(gene_gbM_segments@seqnames@values)
    
    gene_gbM_segments = gene_gbM_segments[unique(queryHits(findOverlaps(gene_gbM_segments, CG_site_ranges)))]

    segment_site_stats2 = data.frame()
    
    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
    for (this_segment in 1:length(gene_gbM_segments)) {
      this_segment.gr = gene_gbM_segments[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats2$start = as.numeric(segment_site_stats2$start)
    segment_site_stats2$end = as.numeric(segment_site_stats2$end)
    segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
    segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
    
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
    print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900) + theme(legend.title = element_blank()))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats2, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}


accession_gene_status = read.table(file="m1001_gene_gbm_tem_umr_calls_2.txt", sep="\t", header=TRUE)


# Triangle plot of proportion of accessions TEM, gbM, UM for each gene

tem_counts = NULL
gbm_counts = NULL
um_counts = NULL

for (this_gene in 1:nrow(accession_gene_status)) {
  status_counts = table(t(relevant_accession_gene_status[this_gene,2:ncol(relevant_accession_gene_status)]))
  if (!("T" %in% names(status_counts))) {
    if ("B" %in% names(status_counts)) {
      tem_counts = c(tem_counts, status_counts[names(status_counts)=="B"])
    } else {
      tem_counts = c(tem_counts, 0)
    }
  } else {
    if ("B" %in% names(status_counts)) {
      tem_counts = c(tem_counts, status_counts[names(status_counts)=="B"] + status_counts[names(status_counts)=="T"])
    } else {
      tem_counts = c(tem_counts, status_counts[names(status_counts)=="T"])
    }
  }
  if (!("G" %in% names(status_counts))) {
    if ("B" %in% names(status_counts)) {
      gbm_counts = c(gbm_counts, status_counts[names(status_counts)=="B"])
    } else {
      gbm_counts = c(gbm_counts, 0)
    }
  } else {
    if ("B" %in% names(status_counts)) {
      gbm_counts = c(gbm_counts, status_counts[names(status_counts)=="B"] + status_counts[names(status_counts)=="G"])
    } else {
      gbm_counts = c(gbm_counts, status_counts[names(status_counts)=="G"])
    }
  }
  if (!("U" %in% names(status_counts))) {
    um_counts = c(um_counts, 0)
  } else {
    um_counts = c(um_counts, status_counts[names(status_counts)=="U"])
  }
}

library(ggtern)

gene_status_counts = as.data.frame(cbind(tem_counts, gbm_counts, um_counts))
gene_status_counts$perc_gbm = gene_status_counts$gbm_counts*100/891
gene_status_counts$perc_tem = gene_status_counts$tem_counts*100/891
gene_status_counts$perc_um = gene_status_counts$um_counts*100/891

ggtern(data=gene_status_counts[!(gene_status_counts$tem_counts==0 & gene_status_counts$gbm_counts==0 & gene_status_counts$um_counts==0),],aes(x=perc_tem,y=perc_gbm, z=perc_um)) + geom_point(alpha=0.1) + theme_showarrows()

ggtern(data=gene_status_counts,aes(x=perc_tem,y=perc_gbm, z=perc_um)) + geom_point() + geom_density_tern(aes(color=..level..),bins=5)
ggtern(data=gene_status_counts,aes(x=tem_counts,y=gbm_counts, z=um_counts)) + geom_point() + geom_density_tern(aes(color=..level..),bins=5)
ggtern(data=gene_status_counts,aes(x=tem_counts,y=gbm_counts, z=um_counts)) + geom_point() + stat_density_tern(geom='polygon', aes(fill=..level..), bins=5, color='grey', bdl=0.001)

# plot Zaigham's density stuff
zz=read.table(file="UM_gbm_TEm_conservation.txt", sep="\t", header=TRUE)
install.packages("viridis")  # Install
library("viridis")  
ggplot(zz) + geom_point(aes(x=UM_conservation, y=Tem_conservation, colour=gbM_conservation)) +scale_color_viridis()
ggplot(zz) + geom_point(aes(x=UM_conservation, y=gbM_conservation, colour=Tem_conservation)) +scale_color_viridis()


zz=read.table(file="UM_gbm_tem_conservation_across_948_accessions.txt", sep="\t", header=TRUE)
zz=zz[zz$Total>0.2*948,c(1,8,6,7)]
zz=zz[sample(nrow(zz), 2000),]
zzz=reshape2::melt(cbind(zz, "bar_order"=-10000000*(zz$Tem_conservation/(zz$gbM_conservation))-zz$gbM_conservation), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")

gbm/(gbm+tem):
  zzz=reshape2::melt(cbind(zz, "bar_order"=zz$gbM_conservation/(zz$gbM_conservation+zz$Tem_conservation)), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")
gbm then (gbm+tem):
  zzz=reshape2::melt(cbind(zz, "bar_order"=(zz$gbM_conservation+zz$Tem_conservation) + 1000*zz$gbM_conservation), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")
gbm-tem:
zzz=reshape2::melt(cbind(zz, "bar_order"=zz$gbM_conservation-zz$Tem_conservation), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")
tem then gbm:
zzz=reshape2::melt(cbind(zz, "bar_order"=zz$gbM_conservation + 1000*zz$Tem_conservation), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")
tem-um:
  zzz=reshape2::melt(cbind(zz, "bar_order"=zz$Tem_conservation-zz$UM_conservation), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")
tem then um:
  zzz=reshape2::melt(cbind(zz, "bar_order"=zz$Tem_conservation*1000 +zz$UM_conservation), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")
um then tem:
  zzz=reshape2::melt(cbind(zz, "bar_order"=zz$UM_conservation*1000 +zz$Tem_conservation), id.vars=c("gene_ID", "bar_order"), variable.name="methylation")

ggplot(zzz, aes(x=reorder(gene_ID, bar_order), y=value, fill=methylation)) + geom_bar(position="fill", stat="identity")


# investigate coverage of At1G19410
this_gene_range = makeGRangesFromDataFrame(gff.genes[gff.genes$gene_ID=="AT1G19410",], seqnames.field = "V1", start.field = "V4", end.field = "V5")
# Fix chromosome names in gene_ranges
levels(this_gene_range@seqnames@values) = substr(as.character(levels(this_gene_range@seqnames@values)),4,4)
this_gene_range@seqinfo@seqnames = levels(this_gene_range@seqnames@values)

View(as.data.frame(CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, this_gene_range))],))
View(all_samples_meth_status[queryHits(findOverlaps(CG_site_ranges, this_gene_range)),])


# There are a number of segments with 100% methylation in the experimental data. On inspection these look like transposon methylation in Schmitz data.
# We are going to filter out any genes which are TEM in more than 10% of Col-0-like accessions

accession_gene_status = read.table(file="m1001_gene_gbm_tem_umr_calls_5.txt", sep="\t", header=TRUE)

relevant_accession_gene_status = accession_gene_status[,c("gene_ID",relevant_samples)]


max_te_count = length(relevant_samples)*0.1
non_te_genes = NULL
for (this_gene in 1:nrow(relevant_accession_gene_status)) {
  status_counts = table(t(relevant_accession_gene_status[this_gene,2:ncol(relevant_accession_gene_status)]))
  if (!("T" %in% names(status_counts))) {
    non_te_genes = c(non_te_genes, relevant_accession_gene_status$gene_ID[this_gene])
  } else {
    if (status_counts[names(status_counts)=="T"]<max_te_count) {
      non_te_genes = c(non_te_genes, relevant_accession_gene_status$gene_ID[this_gene])
    } 
  }
}




# For each annotation - take the ones prepared above, and filter out any genes with TEM in more than X% of Col-0-like accessions
for (min_quantile_cov in c(0.2, 0.5, 0.1, 0.05)) {
  #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 40, 30)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
      #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_strict_spacing_",min_site_spacing,"_2.tsv"), sep="\t", header=TRUE)

    cat(paste0(min_quantile_cov, " ", min_site_spacing, " ", nrow(these_gbM_segments))) 
    these_gbM_segments = these_gbM_segments[these_gbM_segments$gene_ID %in% non_te_genes,]
    cat(paste0(" ", nrow(these_gbM_segments), "\n"))
    write.table(these_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_strict_non_te_spacing_",min_site_spacing,"_2.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="Chromosome", start.field="start", end.field="end")
    # do the bit not present in the Schmitz annotation
    #gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)
    gene_gbM_segments = these_gbM_segments.gr
    
    # Fix chromosome names in gene_gbM_segments
    levels(gene_gbM_segments@seqnames@values) = substr(as.character(levels(gene_gbM_segments@seqnames@values)),4,4)
    gene_gbM_segments@seqinfo@seqnames = levels(gene_gbM_segments@seqnames@values)
    
    gene_gbM_segments = gene_gbM_segments[unique(queryHits(findOverlaps(gene_gbM_segments, CG_site_ranges)))]
    
    segment_site_stats2 = data.frame()
    
    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
    for (this_segment in 1:length(gene_gbM_segments)) {
      this_segment.gr = gene_gbM_segments[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats2$start = as.numeric(segment_site_stats2$start)
    segment_site_stats2$end = as.numeric(segment_site_stats2$end)
    segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
    segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
    
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
    print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900) + theme(legend.title = element_blank()))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats2, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}


# This all went very well, but we now have no segments with zero gbM in Col-0. To put them back in we now want:
# Segments less dense than threshold but present in 0.5% of accessions, no TEs in more than 10% of accessions
# Went back to the 'Col-0_or_not' set, before it went strict 

paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_spacing_",min_site_spacing,"_2.tsv")

# For each annotation - take the ones prepared above, and filter out any genes with TEM in more than X% of Col-0-like accessions
for (min_quantile_cov in c(0.2, 0.5, 0.1, 0.05)) {
  #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 40, 30)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
      #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_spacing_",min_site_spacing,"_2.tsv"), sep="\t", header=TRUE)
    
    cat(paste0(min_quantile_cov, " ", min_site_spacing, " ", nrow(these_gbM_segments))) 
    these_gbM_segments = these_gbM_segments[these_gbM_segments$gene_ID %in% non_te_genes,]
    cat(paste0(" ", nrow(these_gbM_segments), "\n"))
    write.table(these_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_non_te_spacing_",min_site_spacing,"_2.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="Chromosome", start.field="start", end.field="end")
    # do the bit not present in the Schmitz annotation
    #gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)
    gene_gbM_segments = these_gbM_segments.gr
    
    # Fix chromosome names in gene_gbM_segments
    levels(gene_gbM_segments@seqnames@values) = substr(as.character(levels(gene_gbM_segments@seqnames@values)),4,4)
    gene_gbM_segments@seqinfo@seqnames = levels(gene_gbM_segments@seqnames@values)
    
    gene_gbM_segments = gene_gbM_segments[unique(queryHits(findOverlaps(gene_gbM_segments, CG_site_ranges)))]
    
    segment_site_stats2 = data.frame()
    
    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_segments)
    for (this_segment in 1:length(gene_gbM_segments)) {
      this_segment.gr = gene_gbM_segments[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats2$start = as.numeric(segment_site_stats2$start)
    segment_site_stats2$end = as.numeric(segment_site_stats2$end)
    segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
    segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
    
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
    print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900) + theme(legend.title = element_blank()))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats2, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}


# This probably is good (awaiting modelling) . we now want to find out whether it is just as good if we allow multiple segments within a gene to be merged into a single trans segment. we now want:
# Segments less dense than threshold but present in 0.5% of accessions, no TEs in more than 10% of accessions
# Went back to the 'Col-0_or_not' set, before it went strict, and then allowed multiple segments per gene to be merged by a spanning bridge segment

paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_spacing_",min_site_spacing,"_2.tsv")

# For each annotation - take the ones prepared above, and filter out any genes with TEM in more than X% of Col-0-like accessions
for (min_quantile_cov in c(0.2, 0.5, 0.1, 0.05)) {
  #min_quantile_cov=0.05
  #for (min_dense_quantile_cov in 2:3) {
  #for (min_site_spacing in c(20, 40, 30)) {
  for (min_site_spacing in c(33, 37, 43, 47, 50, 60)) {
      #for (min_quantile_cov in 1:9) {
    min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
    #min_dense_seg_cover = round(min_dense_quantile_cov*length(relevant_samples)/10)
    
    #for (min_seg_cover in 5) {
    #gene_firsts = rep(0, length(gene_ranges))
    #gene_lasts = rep(0, length(gene_ranges))
    #gene_gbM_segments = data.frame()
    
    these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_non_te_spacing_",min_site_spacing,"_2.tsv"), sep="\t", header=TRUE)
    
    cat(paste0(min_quantile_cov, " ", min_site_spacing, " ", nrow(these_gbM_segments))) 
    #these_gbM_segments = these_gbM_segments[these_gbM_segments$gene_ID %in% non_te_genes,]

    # for each gene
    # find segments from these_gbM_segments which overlap
    # find the minimum start site and maximum end site
    # if there are any segments, make a new spanning segment and add to the new_gbM_segments
    
    cat(paste0(" ", nrow(new_gbM_segments), "\n"))
    
    these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="Chromosome", start.field="start", end.field="end")

    # do the bit not present in the Schmitz annotation
    #gene_gbM_segments = setdiff(these_gbM_segments.gr, GBM_segments)
    gene_gbM_segments = these_gbM_segments.gr
    
    # Fix chromosome names in gene_gbM_segments
    levels(gene_gbM_segments@seqnames@values) = substr(as.character(levels(gene_gbM_segments@seqnames@values)),4,4)
    gene_gbM_segments@seqinfo@seqnames = levels(gene_gbM_segments@seqnames@values)
    
    gene_gbM_segments = gene_gbM_segments[unique(queryHits(findOverlaps(gene_gbM_segments, CG_site_ranges)))]
    
    gene_gbM_trans_segments = data.frame()
    for (this_gene in 1:length(gene_ranges)) {
      if(gene_locus_type[this_gene]=="protein_coding") {
        cat(paste0("cover: ", min_seg_cover, "  gene: ", this_gene, "\n"))
        doing_one = FALSE
        this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
        this_start = 0
        this_end = 0

        gene_segment_olaps = findOverlaps(gene_ranges[this_gene], gene_gbM_segments)

        save_this = FALSE
        if (length(gene_segment_olaps) > 0) {
          save_this = TRUE
          this_start = min(gene_gbM_segments[subjectHits(gene_segment_olaps)]@ranges@start)
          this_end = max(gene_gbM_segments[subjectHits(gene_segment_olaps)]@ranges@start + gene_gbM_segments[subjectHits(gene_segment_olaps)]@ranges@width -1)
          if (this_start<gene_ranges[this_gene]@ranges@start) {
            this_start = gene_ranges[this_gene]@ranges@start
          }
          if (this_end>(gene_ranges[this_gene]@ranges@start+gene_ranges[this_gene]@ranges@width-1)) {
            this_end = gene_ranges[this_gene]@ranges@start+gene_ranges[this_gene]@ranges@width-1
          }
        }
        
        if (save_this) {
          gene_gbM_trans_segments = rbind(gene_gbM_trans_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=this_start, "end"=this_end))
          colnames(gene_gbM_trans_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
          save_this=FALSE
        }
        
      }
    }
    write.table(gene_gbM_trans_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_Col-0_or_not_non_te_spacing_",min_site_spacing,"_trans_2.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    

    segment_site_stats2 = data.frame()
    gene_gbM_trans_segments.gr = makeGRangesFromDataFrame(gene_gbM_trans_segments, seqnames.field="Chromosome", start.field="start", end.field="end")
    
    segment_site_olaps = findOverlaps(CG_site_ranges, gene_gbM_trans_segments.gr)
    for (this_segment in 1:length(gene_gbM_trans_segments.gr)) {
      this_segment.gr = gene_gbM_trans_segments.gr[this_segment]
      these_sites = queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])
      gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps[subjectHits(segment_site_olaps)==this_segment])]
      if (length(these_sites)>0) {
        #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_trans_segments.gr[this_segment]@ranges@start, "end" = gene_gbM_trans_segments.gr[this_segment]@ranges@start+gene_gbM_trans_segments.gr[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
        colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
      }
    }
    #for (this_segment in 1:length(gene_gbM_segments)) {
    #  this_segment.gr = gene_gbM_segments[this_segment]
    #  segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    #  these_sites = queryHits(segment_site_olaps)
    #  gene_site_counts = CG_site_gbM_counts_Col0like[queryHits(segment_site_olaps)]
    #  if (length(these_sites)>0) {
    #    #segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    segment_site_stats2 = rbind(segment_site_stats2, c("start" = gene_gbM_segments[this_segment]@ranges@start, "end" = gene_gbM_segments[this_segment]@ranges@start+gene_gbM_segments[this_segment]@ranges@width-1, "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
    #    colnames(segment_site_stats2) = c("start", "end", "no_sites", "mean_coverage")
    #  }
    #}
    segment_site_stats2$start = as.numeric(segment_site_stats2$start)
    segment_site_stats2$end = as.numeric(segment_site_stats2$end)
    segment_site_stats2$no_sites = as.numeric(segment_site_stats2$no_sites)
    segment_site_stats2$mean_coverage = as.numeric(segment_site_stats2$mean_coverage)
    
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))))
    #print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  min. dense quantile.cov.: ", min_dense_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900))
    print(ggplot(segment_site_stats2) +geom_point(aes(x=(1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), y=segment_site_stats2$mean_coverage, colour=log(1+segment_site_stats2$end-segment_site_stats2$start)))  +scale_colour_gradient2() +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("bp/CG threshold: ", min_site_spacing,"  min. quantile.cov.: ", min_quantile_cov,"  R=", cor((1+segment_site_stats2$end-segment_site_stats2$start)/(segment_site_stats2$no_sites-1), segment_site_stats2$mean_coverage))) +xlim(0,100) +ylim(0,900) + theme(legend.title = element_blank()))
    
    #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats2, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_min_non_Col-0_spacing_",min_site_spacing,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
    # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
    #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
    #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}


# Check for overlaps among gbM segments
# I needed to do this iteratively (three times) because the first run left about 20 entries behind that were still duplicated, and the second run left 3 behind

for (min_site_spacing in c(30, 33, 37, 43, 47, 50, 60)) {
  
these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_0.05_min_CG_3_min_Col-0_or_not_non_te_spacing_",min_site_spacing,"_trans_no_olaps_3.tsv"), sep="\t", header=TRUE)
these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="Chromosome", start.field="start", end.field="end")
length(these_gbM_segments.gr)
#[1] 13080
length(reduce(these_gbM_segments.gr))
#[1] 11948

# Find genes where segment in one completely subsumes segment in the other (by 90% or more)
# Drop the shorter of the pair

# Find genes where segments still overlap. Remove the overlapping part from both segments
length_tolerance = 0.9
discarded_segments = NULL
retained_segments = NULL
for (this_segment in 1:length(these_gbM_segments.gr)) {
  seg_olaps = findOverlaps(these_gbM_segments.gr[this_segment], these_gbM_segments.gr)
  if (length(seg_olaps) > 1) {
    cat(paste0(this_segment, "\n"))
    this_seg_done = FALSE
    for (this_olap in 1:length(seg_olaps)) {
      if (subjectHits(seg_olaps[this_olap])!=this_segment) {
        #find length of olap
        olap_segment = intersect(these_gbM_segments.gr[this_segment], these_gbM_segments.gr[subjectHits(seg_olaps[this_olap])])
        #find lengths of 2 segments overlapping
        if (these_gbM_segments.gr[this_segment]@ranges@width > these_gbM_segments.gr[subjectHits(seg_olaps[this_olap])]@ranges@width) {
          #this_segment is longer
          shorter_width = these_gbM_segments.gr[subjectHits(seg_olaps[this_olap])]@ranges@width
          shorter_seg = subjectHits(seg_olaps[this_olap])
        } else {
          if (these_gbM_segments.gr[this_segment]@ranges@width < these_gbM_segments.gr[subjectHits(seg_olaps[this_olap])]@ranges@width) {
            #this_segment is shorter
            shorter_width = these_gbM_segments.gr[this_segment]@ranges@width
            shorter_seg = this_segment
          } else {
            #they must be the same length, in which case, retain the first and discard the second
            if (this_segment>subjectHits(seg_olaps[this_olap])) {
              shorter_width = these_gbM_segments.gr[this_segment]@ranges@width
              shorter_seg = this_segment
            } else {
              shorter_width = these_gbM_segments.gr[subjectHits(seg_olaps[this_olap])]@ranges@width
              shorter_seg = subjectHits(seg_olaps[this_olap])
            }
          }
        }
        #if length of olap > 0.9*length of shorter segment
        if (shorter_width * length_tolerance < olap_segment@ranges@width) {
          # If we are doing the longer segment this time round, add it to the retained_segments
          if (shorter_seg != this_segment) {
            retained_segments = rbind(retained_segments, these_gbM_segments[this_segment,])
            #add shorter segment to discarded_segments
            discarded_segments = c(discarded_segments, shorter_seg)
          }
        } else {
          # Overlap is less than 90% of length of shorter segment, so we need to trim this_segment and save it
          if (!this_seg_done) {
            this_modified_segment = these_gbM_segments[this_segment,]
            this_new_segment = setdiff(these_gbM_segments.gr[this_segment], these_gbM_segments.gr[subjectHits(seg_olaps[this_olap])])
            this_modified_segment$start = this_new_segment@ranges@start
            this_modified_segment$end = this_new_segment@ranges@start + this_new_segment@ranges@width -1
            retained_segments = rbind(retained_segments, this_modified_segment)
            this_seg_done = TRUE
          }
        }
      }
    }
  } else {
    # there are no overlaps to other segments
    retained_segments = rbind(retained_segments, these_gbM_segments[this_segment,])
  }
}
write.table(retained_segments, file=paste0("gene_gbm_segments_min_acc_decile_0.05_min_CG_3_min_Col-0_or_not_non_te_spacing_",min_site_spacing,"_trans_no_olaps_4.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


# Remove anything shorter than 100bp

for (min_site_spacing in c(30, 33, 37, 43, 47, 50, 60)) {
  
  these_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_0.05_min_CG_3_min_Col-0_or_not_non_te_spacing_",min_site_spacing,"_trans_no_olaps_5.tsv"), sep="\t", header=TRUE)
  these_gbM_segments.gr = makeGRangesFromDataFrame(these_gbM_segments, seqnames.field="Chromosome", start.field="start", end.field="end")
  length(these_gbM_segments.gr)
  #[1] 13080
  length(reduce(these_gbM_segments.gr))
  #[1] 11948
  
  # Find genes where segment in one completely subsumes segment in the other (by 90% or more)
  # Drop the shorter of the pair
  
  # Find genes where segments still overlap. Remove the overlapping part from both segments
  min_seg_length = 100
  discarded_segments = NULL
  retained_segments = NULL
  for (this_segment in 1:length(these_gbM_segments.gr)) {
    if (these_gbM_segments.gr[this_segment]@ranges@width >= min_seg_length) {
      retained_segments = rbind(retained_segments, these_gbM_segments[this_segment,])
    } else {  
      discarded_segments = c(discarded_segments, this_segment)
    }
  }
  write.table(retained_segments, file=paste0("gene_gbm_segments_min_acc_decile_0.05_min_CG_3_min_Col-0_or_not_non_te_spacing_",min_site_spacing,"_trans_no_olaps_min_",min_seg_length,"_5.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}




# We now want to extract some segments as comparators for chromatin features - segments that are rarely methylated

# We traverse CG_site_gbM_counts one gene at a time (which already contains counts for each site of overlapping gbM segments longer than 3 sites) and identify contiguous segments within each gene which are covered by relevant gbM segments in less than 5 accessions

for (max_seg_cover in 1:10) {
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  gene_gbM_segments = data.frame()
  for (this_gene in 1:length(gene_ranges)) {
    cat(paste0("cover: ", max_seg_cover, "  gene: ", this_gene, "\n"))
    doing_one = FALSE
    this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
    this_start = 0
    this_end = 0
    gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
    gene_site_counts = CG_site_gbM_counts[subjectHits(gene_site_olaps)]
    no_sites = length(gene_site_olaps)
    if (no_sites>0) {
      if (max(gene_site_counts) <= max_seg_cover) {
        gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start" = gff.genes[this_gene,"V4"], "end" = gff.genes[this_gene,"V5"]))
        colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
      }
    }
  }
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_max_accesions_",max_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}




# Load in bigWig files, and export contents as .RDS.  This must be done on cluster as Windows cannot read .bigWig files
#library(rtracklayer)
#saveRDS(import("/media/sf_D_DRIVE/Reference_data/H2AZ.devin_plos_genetics.w50.Chr.bw", format="BigWig")Error in seqinfo(ranges) : UCSC library operation failed
#saveRDS(import("H2AZ.devin_plos_genetics.w50.Chr.bw", format="BigWig"), "H2AZ.devin_plos_genetics.w50.C> saveRDS(import("k27me1.col.GSM3040071.bw", format="BigWig"), "k27me1.col.GSM3040071.RDS")              >


# Load in the .RDS files as genomic ranges with scores
H2AZ_positions = readRDS("../0-reference/chipseq/H2AZ.devin_plos_genetics.w50.Chr.RDS")
H3K27me1_positions = readRDS("../0-reference/chipseq/k27me1.col.GSM3040071.RDS")

# Fix chromosome names in chipseq objects
levels(H2AZ_positions@seqnames@values) = substr(as.character(levels(H2AZ_positions@seqnames@values)),4,4)
H2AZ_positions@seqinfo@seqnames = levels(H2AZ_positions@seqnames@values)
levels(H3K27me1_positions@seqnames@values) = substr(as.character(levels(H3K27me1_positions@seqnames@values)),4,4)
H3K27me1_positions@seqinfo@seqnames = levels(H3K27me1_positions@seqnames@values)

gene_mean_h3k27me1 = NULL
gene_mean_h2az = NULL
for (this_gene in 1:nrow(gff.genes)) {
  cat(paste0(gff.genes[this_gene,"gene_ID"], "\n"))
  gene_mean_h3k27me1 = rbind(gene_mean_h3k27me1, c("gene_ID"=gff.genes$gene_ID[this_gene], "mean_h3k27me1"=mean(as.data.frame(H3K27me1_positions)[subjectHits(findOverlaps(gene_ranges[this_gene], H3K27me1_positions)),]$score)))
  gene_mean_h2az = rbind(gene_mean_h2az, c("gene_ID"=gff.genes$gene_ID[this_gene], "mean_h2az"=mean(as.data.frame(H2AZ_positions)[subjectHits(findOverlaps(gene_ranges[this_gene], H2AZ_positions)),]$score)))
}







# Let's look at the 100 most Col-0-like accessions
schmitz_gbM_density = sum(schmitz_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"]=="M", na.rm=TRUE) / sum(merged_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"] %in% c("M","U"), na.rm=TRUE)
#[1] 0.6124375

sorted_accessions = sort(accession_gbm_levels)

Col_0_like_100 = names(c(head(sorted_accessions[sorted_accessions>schmitz_gbM_density], 50), tail(sorted_accessions[sorted_accessions<schmitz_gbM_density], 50)))


# CG_site_gbM_counts contains counts for all accessions. We want it just for Col0-like ones
CG_site_gbM_counts_Col0like100 = rep(0,length(CG_site_ranges))

# at this point, check that CG_site_ranges have chromosomes identified just as numbers rather than ChrX
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if ((this_sample %in% Col_0_like_100)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)
    CG_site_gbM_counts_Col0like100[subjectHits(model_site_olaps)] = CG_site_gbM_counts_Col0like100[subjectHits(model_site_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}


# We traverse CG_site_gbM_counts_Col0like100 one gene at a time (which already contains counts for each site of overlapping gbM segments longer than 3 sites) and identify contiguous segments within each gene which are covered by relevant gbM segments in at least x% of accessions
# Make sure gene_ranges has numbers rather than ChrX

for (min_quantile_cov in 1:9) {
  min_seg_cover = round(min_quantile_cov*length(Col_0_like_100)/10)
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  gene_gbM_segments = data.frame()
  for (this_gene in 1:length(gene_ranges)) {
    if(gene_locus_type[this_gene]=="protein_coding") {
      cat(paste0("cover: ", min_seg_cover, "  gene: ", this_gene, "\n"))
      doing_one = FALSE
      this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
      this_start = 0
      this_end = 0
      gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
      gene_site_counts = CG_site_gbM_counts_Col0like100[subjectHits(gene_site_olaps)]
      no_sites = length(gene_site_olaps)
      if (no_sites>0) {
        prev_site = 0
        last_site = FALSE
        save_this = FALSE
        for (this_site in 1:no_sites) {
          if (this_site==no_sites) {}
          if (doing_one) {
            if (gene_site_counts[this_site]>=min_seg_cover) {
              if (this_site==no_sites) {
                this_end = this_site
                save_this = TRUE
              }
            } else {
              this_end = this_site -1
              doing_one = FALSE
              save_this = TRUE
            }
          } else {
            if (gene_site_counts[this_site]>=min_seg_cover) {
              doing_one=TRUE
              this_start = this_site
              if (this_site==no_sites) {
                this_end = this_site
                save_this = TRUE
              }
            }
          }
          if (save_this) {
            gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start, "end"=CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1))
            colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
            save_this=FALSE
          }
        }
      } else {
        #gene_firsts[this_gene] = NA
        #gene_lasts[this_gene] = NA
      }
    }
  }
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_Col0100.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])


pdf("annotation_summary_stats_Col0100.pdf")
# For each annotation - plot CG density vs. mCG in Col-0, average number of accessions each CG site is annotated as gbM
for (min_quantile_cov in 1:9) {
  min_seg_cover = round(min_quantile_cov*length(relevant_samples)/10)
  #for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  #gene_gbM_segments = data.frame()
  
  segment_site_stats = data.frame()
  
  gene_gbM_segments = read.table(file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_Col0100.tsv"), sep="\t", header=TRUE)
  
  for (this_segment in 1:nrow(gene_gbM_segments)) {
    this_segment.gr = makeGRangesFromDataFrame(gene_gbM_segments[this_segment,], seqnames.field="chromosome", start.field="start", end.field="end")
    segment_site_olaps = findOverlaps(CG_site_ranges, this_segment.gr)
    these_sites = queryHits(segment_site_olaps)
    gene_site_counts = CG_site_gbM_counts_Col0like100[queryHits(segment_site_olaps)]
    if (no_sites>0) {
      segment_site_stats = rbind(segment_site_stats, c("gene_ID"=gene_gbM_segments[this_segment,"gene_ID"], "chromosome"=gene_gbM_segments[this_segment,"chromosome"], "start" = gene_gbM_segments[this_segment,"start"], "end" = gene_gbM_segments[this_segment,"end"], "no_sites" = length(these_sites), "mean_coverage" = mean(gene_site_counts)))
      colnames(segment_site_stats) = c("gene_ID", "chromosome", "start", "end", "no_sites", "mean_coverage")
    }
  }
  segment_site_stats$start = as.numeric(segment_site_stats$start)
  segment_site_stats$end = as.numeric(segment_site_stats$end)
  segment_site_stats$no_sites = as.numeric(segment_site_stats$no_sites)
  segment_site_stats$mean_coverage = as.numeric(segment_site_stats$mean_coverage)
  
  print(ggplot(segment_site_stats) +geom_point(aes(x=(1+segment_site_stats$end-segment_site_stats$start)/segment_site_stats$no_sites, y=segment_site_stats$mean_coverage)) +xlab("segment CG site spacing") +ylab("No. accessions with gbM in segment") +ggtitle(paste0("R=", cor((1+segment_site_stats$end-segment_site_stats$start)/segment_site_stats$no_sites, segment_site_stats$mean_coverage))))
  
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(segment_site_stats, file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_Col0100.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  #segment_site_stats = read.table(file=paste0("segment_site_stats_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,".tsv"), sep="\t", header=TRUE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  #gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  #write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_acc_decile_",min_quantile_cov,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}
dev.off()









reference_gff = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff")  # AraPort11 annotation

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

gene_gbm_olaps = intersect(gene_ranges, overlap_model, ignore.strand=TRUE)

# genes with at least one gbM segment
candidate_genes = gff.genes[unique(queryHits(findOverlaps(gene_ranges, gene_gbm_olaps))),]
candidate_genes$gbm_start = rep(0, nrow(candidate_genes))
candidate_genes$gbm_end = rep(0, nrow(candidate_genes))

for (this_gene in 1:nrow(candidate_genes)) {
  cat(paste0(this_gene,"\n"))
  candidate_segments = intersect(makeGRangesFromDataFrame(candidate_genes[this_gene,c(1,4,5)], seqnames.field = "V1", start.field = "V4", end.field = "V5"), gene_gbm_olaps)
  candidate_genes[this_gene, "gbm_start"] = range(candidate_segments)@ranges@start
  candidate_genes[this_gene, "gbm_end"] = range(candidate_segments)@ranges@start + range(candidate_segments)@ranges@width -1
  
  #candidate_segments = unique(subjectHits(findOverlaps(makeGRangesFromDataFrame(candidate_genes[this_gene,c(1,4,5)], seqnames.field = "V1", start.field = "V4", end.field = "V5"), gene_gbm_olaps)))
  #candidate_genes[this_gene, "gbm_start"] = range(gene_gbm_olaps[candidate_segments])@ranges@start  
  #candidate_genes[this_gene, "gbm_end"] = range(gene_gbm_olaps[candidate_segments])@ranges@start + range(gene_gbm_olaps[candidate_segments])@ranges@width -1
}

candidate_genes = candidate_genes[,c(1,4,5,10,12,13)]
colnames(candidate_genes)[1:3] = c("chromosome","start","end")

write.table(candidate_genes, file="gene_trans_accession_gbm_segments_v2.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


# New version - mark up, for each CG site, the percentage of accessions in which it is methylated

all_samples_meth_status = readRDS(file="m1001_CG_all_samples_meth_status_0.005_0.05_10_25_SNPs.rds")

site_m_count = rowSums(all_samples_meth_status[,4:ncol(all_samples_meth_status)]=="M", na.rm=TRUE)
site_u_count = rowSums(all_samples_meth_status[,4:ncol(all_samples_meth_status)]=="U", na.rm=TRUE)

# make granges of all CG sites
CG_site_ranges =makeGRangesFromDataFrame(df = all_samples_meth_status, start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
gene_site_ranges = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, gene_ranges))]
gene_site_m_count = rep(0, length(gene_site_ranges))

# load in all the segmentation models, and create the union of all the gbM segments in each gene
no_models = 0
ignore_list=c()

for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 10,19)
  if (substr(this_sample,10,10)=="_") {this_sample=substr(this_sample,1,9)}
  if (!(this_sample %in% ignore_list)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model$Chromosome = paste0("Chr", as.character(this_model$Chromosome))
    this_model = this_model[substr(this_model$Type,1,4)=="GBM_",]
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    msite_ranges = makeGRangesFromDataFrame(df = all_samples_meth_status[(!is.na(all_samples_meth_status[,this_sample])) & (all_samples_meth_status[,this_sample]=="M"),], start.field = "Locus", end.field = "Locus", seqnames.field = "Chromosome")
    gbM_msites = msite_ranges[queryHits(findOverlaps(msite_ranges, this_model.gr))]
    gene_site_gbm_msite_olaps = findOverlaps(gene_site_ranges, gbM_msites)
    gene_site_m_count[queryHits(gene_site_gbm_msite_olaps)] = gene_site_m_count[queryHits(gene_site_gbm_msite_olaps)] +1

    #    if (no_models == 1) {
    #  overlap_model = this_model.gr
    #} else {
    #  overlap_model = union(overlap_model, this_model.gr, ignore.strand=TRUE)
    #}
  }
}

write.table(gene_site_m_count, file="gene_site_m_count.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


z=cbind(as.data.frame(CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, gene_ranges))]), gene_site_m_count, gff.genes$gene_ID[subjectHits(findOverlaps(CG_site_ranges, gene_ranges))], gff.genes$gene_name[subjectHits(findOverlaps(CG_site_ranges, gene_ranges))])
colnames(z)[c(7,8)] = c("gene_ID", "gene_name")
# This is a proxy for the protein coding subset (ones whose gene_name begins with AT)
View(z[substr(z$gene_name,1,2)=="AT",])
write.table(z[substr(z$gene_name,1,2)=="AT",], file="gene_site_m_count.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Add "Chr" to the chromosome names
levels(overlap_model@seqnames@values) = paste0("Chr",as.character(levels(overlap_model@seqnames@values)))
overlap_model@seqinfo@seqnames = levels(overlap_model@seqnames@values)


for each gene {
  first_meth_site = 0
  last_meth_site = 0
  for this_site in gene {
    if site_m_count[this_site] > min_m_threshold {
      if (first_meth_site == 0) {
        first_meth_site = this_site
      }
      last_meth_site = this_site
    }    
  }
  save first_site and last_site
}


for each accession {
  read in gbM segments
  overlap gbm Segments with CG sites
  update accession count for that site
}





for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 10,19)
  if (substr(this_sample,10,10)=="_") {this_sample=substr(this_sample,1,9)}
  if (!(this_sample %in% ignore_list)) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    if (no_models == 1) {
      start_sites = cbind(this_model$Chromosome, this_model$Start)
      end_sites = cbind(this_model$Chromosome, this_model$End)
    } else {
      #start_sites = rbind(start_sites, cbind(this_model$Chromosome, this_model$Start))
      end_sites = rbind(end_sites, cbind(this_model$Chromosome, this_model$End))
      
    }
  }
}

# We are going to construct a set of epihaplotype blocks, each with a unique pattern of uM/gbM/TEM across accessions in comparison to its neighbours

# Make a list of unique end sites for segments
#unique_start_sites = unique(start_sites)
unique_blocks = unique(end_sites)
unique_blocks = as.data.frame(unique_blocks)
colnames(unique_blocks) = c("Chromosome", "End")
unique_blocks = unique_blocks[order(unique_blocks$Chromosome,unique_blocks$End),]

# Add blank start sites for each end site
unique_blocks = cbind(unique_blocks$Chromosome, "Start" = rep(0,nrow(unique_blocks)), unique_blocks$End)
unique_blocks = as.data.frame(unique_blocks)
colnames(unique_blocks) = c("Chromosome", "Start", "End")

# use the lag() function from dplyr to assign the end+1 of each previous row as the start of each next row
#install.packages("dplyr")
library(dplyr)
unique_blocks$Start=lag(unique_blocks$End+1)

# Correct the first row, and the first rows of each subsequent chromosome
unique_blocks[1,"Start"] = 1
unique_blocks[unique_blocks$Start>unique_blocks$End,"Start"] = 1

# overlap the epihaplotype blocks with the segment model in each accession and create a set of segment types for each epihaplotype block for the accession

library(GenomicRanges)

# do one chromosome at a time to avoid over-running memory
for (this_chrom in 1:5) {
  this_chrom_name=this_chrom
  # Make a gr for the epihaplotype blocks
  unique_blocks.gr = makeGRangesFromDataFrame(df=unique_blocks[unique_blocks$Chromosome==this_chrom_name,], start.field="Start", end.field="End", seqnames.field="Chromosome")
  
  no_models = 0
  for (this_model_file in list.files(path=model_dir)) {
    this_sample = substr(this_model_file, 10,19)
    if (substr(this_sample,10,10)=="_") {this_sample=substr(this_sample,1,9)}
    if (!(this_sample %in% ignore_list)) {
      no_models = no_models + 1
      cat(paste0(this_sample,"\n"))
      this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
      model_haplo_olaps=findOverlaps(makeGRangesFromDataFrame(df=this_model, start.field="Start", end.field="End", seqnames.field="Chromosome"), unique_blocks.gr)
      queryHits(model_haplo_olaps)
      
      if (no_models == 1) {
        epihaplo_patterns = cbind("block"=subjectHits(model_haplo_olaps), unique_blocks[unique_blocks$Chromosome==this_chrom_name,][subjectHits(model_haplo_olaps),], this_sample=substr(this_model[queryHits(model_haplo_olaps),c("Type")],1,3))
      } else {
        epihaplo_patterns = merge(epihaplo_patterns, cbind("block"=subjectHits(model_haplo_olaps), this_sample=substr(this_model[queryHits(model_haplo_olaps),c("Type")],1,3)), by=c("block"))
      }
      colnames(epihaplo_patterns)[no_models+4] = this_sample
    }
  }
  
  write.table(epihaplo_patterns, file=paste0("epihaplo_patterns_",this_chrom,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

