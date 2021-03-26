# Load in met1 and wt methylation data from Choi et al, 2021; compile mean values per gene

pathroot="X:/Daniel-Zilberman"

# Source directory containing tables downloaded from GEO
source_dir = paste0(pathroot,"/Projects/Jay-1001_methylomes/0-reference/supporting/met1")

# Working directory where outputs will be directed to
setwd(dir = paste0(pathroot,"/Projects/Jay-1001_methylomes/0-reference/supporting/met1/"))

reference_fasta = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta")  # TAIR10 genome assembly This never gets used at the moment
reference_CG_sites = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.CG_sites.tsv")   # We prepared this earlier using a perl script to parse the TAIR10 genome assembly
reference_gff = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff")  # AraPort11 annotation
reference_exons = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-exon.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_introns = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-intron.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_5UTR = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-five_prime_UTR.gff") # Prepared from AraPort11 annotation by a script in bs_sequel
reference_3UTR = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/GFF/Araport11_GFF3_genes_transposons.201606-three_prime_UTR.gff") # Prepared from AraPort11 annotation by a script in bs_sequel

library(stringr)

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
# fix chromosome IDs
levels(gene_ranges@seqnames@values) = substr(as.character(levels(gene_ranges@seqnames@values)),4,4)
gene_ranges@seqinfo@seqnames = levels(gene_ranges@seqnames@values)

# load in 50bp window methylation data
col0_CG = read.table(file="GSE122394_col.bisulfite.CG.w50.gff", sep="\t", header=TRUE)
col0_CHG = read.table(file="GSE122394_col.bisulfite.CHG.w50.gff", sep="\t", header=TRUE)
col0_CHH = read.table(file="GSE122394_col.bisulfite.CHH.w50.gff", sep="\t", header=TRUE)
met1_CG = read.table(file="GSE122394_met1.bisulfite.CG.w50.gff", sep="\t", header=TRUE)
met1_CHG = read.table(file="GSE122394_met1.bisulfite.CHG.w50.gff", sep="\t", header=TRUE)
met1_CHH = read.table(file="GSE122394_met1.bisulfite.CHH.w50.gff", sep="\t", header=TRUE)
colnames(col0_CG) = c("chr","dot","context","start","end","fraction","dot2","dot3","stuff")
colnames(col0_CHG) = c("chr","dot","context","start","end","fraction","dot2","dot3","stuff")
colnames(col0_CHH) = c("chr","dot","context","start","end","fraction","dot2","dot3","stuff")
colnames(met1_CG) = c("chr","dot","context","start","end","fraction","dot2","dot3","stuff")
colnames(met1_CHG) = c("chr","dot","context","start","end","fraction","dot2","dot3","stuff")
colnames(met1_CHH) = c("chr","dot","context","start","end","fraction","dot2","dot3","stuff")

col0_CG = col0_CG[,c(1,4,5,6)]
col0_CHG = col0_CHG[,c(1,4,5,6)]
col0_CHH = col0_CHH[,c(1,4,5,6)]
met1_CG = met1_CG[,c(1,4,5,6)]
met1_CHG = met1_CHG[,c(1,4,5,6)]
met1_CHH = met1_CHH[,c(1,4,5,6)]

col0_CG$fraction[col0_CG$fraction=="."] = NA
col0_CHG$fraction[col0_CHG$fraction=="."] = NA
col0_CHH$fraction[col0_CHH$fraction=="."] = NA
met1_CG$fraction[met1_CG$fraction=="."] = NA
met1_CHG$fraction[met1_CHG$fraction=="."] = NA
met1_CHH$fraction[met1_CHH$fraction=="."] = NA

col0_CG$fraction = as.numeric(col0_CG$fraction)
col0_CHG$fraction = as.numeric(col0_CHG$fraction)
col0_CHH$fraction = as.numeric(col0_CHH$fraction)
met1_CG$fraction = as.numeric(met1_CG$fraction)
met1_CHG$fraction = as.numeric(met1_CHG$fraction)
met1_CHH$fraction = as.numeric(met1_CHH$fraction)

all_samples_meth = col0_CG
#all_samples_meth = merge(col0_CG, col0_CHG, col0_CHH, met1_CG, met1_CHG, met1_CHH, by=c("chr","start","end"))
all_samples_meth = merge(col0_CG, col0_CHG, by=c(1:3))
all_samples_meth = merge(all_samples_meth, col0_CHH, by=c(1:3))
all_samples_meth = merge(all_samples_meth, met1_CG, by=c(1:3))
all_samples_meth = merge(all_samples_meth, met1_CHG, by=c(1:3))
all_samples_meth = merge(all_samples_meth, met1_CHH, by=c(1:3))

colnames(all_samples_meth)[4:9] = c("col0_CG","col0_CHG","col0_CHH","met1_CG","met1_CHG","met1_CHH")

meth_ranges = makeGRangesFromDataFrame(all_samples_meth, start.field = "start", end.field = "end", seqnames.field = "chr")
# fix chromosome IDs
levels(meth_ranges@seqnames@values) = substr(as.character(levels(meth_ranges@seqnames@values)),4,4)
meth_ranges@seqinfo@seqnames = levels(meth_ranges@seqnames@values)

gene_meth_olaps = findOverlaps(gene_ranges, meth_ranges)

for (this_col in 4:9) {
  gff.genes = cbind(gff.genes, rep(NA, nrow(gff.genes)))
  colnames(gff.genes)[11+this_col-3] = colnames(all_samples_meth)[this_col]
  for (this_gene in 1:nrow(gff.genes)) {
    gff.genes[this_gene, 11+this_col-3] = mean(all_samples_meth[subjectHits(gene_meth_olaps[queryHits(gene_meth_olaps)==this_gene]),this_col], na.rm=TRUE)
  }
}

write.table(gff.genes[,c(10,12:17)], file="gene_col0_met1_mC.txt", sep="\t",col.names=TRUE, row.names=FALSE)


library(rtracklayer)

#[3:53 PM] Jaemyung Choi (JIC)
# log2 ratios of output to input DNA; output from deeptools2 bamcompare
col0_H3K9me2 = rtracklayer::import("GSM3465913_col.H3K9me2.w50.bw", format="BigWig")
met1_H3K9me2 = rtracklayer::import("GSM3465914_met1.H3K9me2.w50.bw", format="BigWig")

# fix chromosome IDs
levels(col0_H3K9me2@seqnames@values) = substr(as.character(levels(col0_H3K9me2@seqnames@values)),4,4)
col0_H3K9me2@seqinfo@seqnames = levels(col0_H3K9me2@seqnames@values)
levels(met1_H3K9me2@seqnames@values) = substr(as.character(levels(met1_H3K9me2@seqnames@values)),4,4)
met1_H3K9me2@seqinfo@seqnames = levels(met1_H3K9me2@seqnames@values)


gene_chip_olaps = findOverlaps(gene_ranges, col0_H3K9me2)
gff.genes = cbind(gff.genes, rep(NA, nrow(gff.genes)))
colnames(gff.genes)[ncol(gff.genes)] = "col0_H3K9me2"
for (this_gene in 1:nrow(gff.genes)) {
  gff.genes[this_gene, ncol(gff.genes)] = mean(as.data.frame(col0_H3K9me2)[subjectHits(gene_chip_olaps[queryHits(gene_chip_olaps)==this_gene]),"score"], na.rm=TRUE)
}


gene_chip_olaps = findOverlaps(gene_ranges, met1_H3K9me2)
gff.genes = cbind(gff.genes, rep(NA, nrow(gff.genes)))
colnames(gff.genes)[ncol(gff.genes)] = "met1_H3K9me2"
for (this_gene in 1:nrow(gff.genes)) {
  gff.genes[this_gene, ncol(gff.genes)] = mean(as.data.frame(met1_H3K9me2)[subjectHits(gene_chip_olaps[queryHits(gene_chip_olaps)==this_gene]),"score"], na.rm=TRUE)
}

write.table(gff.genes[,c(10,18:19)], file="gene_col0_met1_h3k9ME2.txt", sep="\t",col.names=TRUE, row.names=FALSE)

