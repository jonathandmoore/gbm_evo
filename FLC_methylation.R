# parse all segmentation models, and extract blocks overlapping FLC. Extract methylation level of such blocks

#setwd("W:/Daniel-Zilberman/Jay_1001_methylomes")
#setwd("X:/Daniel-Zilberman/Projects/Jay_1001_methylomes/5-analysis")
setwd("/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis")



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


FLC_ranges = AT5G10140


model_dir = "segmentation_models"
# make an empty structure to hold the resulting model
merged_model = data.frame()
no_models = 0
ignore_list=c()

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

