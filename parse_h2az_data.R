# load in H2A.Z data from CHIP-Seq or CHIP-Chip, load in an annotation and a modelling performance run. Use the H2AZ data to explore the performance split by triads or deciles, and set a smoothed H2A.Z threshold to trim the annotation segments.
 
# met1 Data is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12212 which was part of Zilberman et al, 2008 Nature https://pubmed.ncbi.nlm.nih.gov/18815594/

#Paper uses two colour microarrays to assay CHIP for wt and met1 (among other stuff)

#GSM307377.gff	Zilberman H2A.Z IP/total DNA 1823502 (Ab-2 WT)
#GSM307378.gff	Zilberman H2A.Z IP/total DNA 1824902 (Ab-1 WT)
#GSM307379.gff	Zilberman H2A.Z IP/total DNA 1825002 (Ab-2 met1-6)
#GSM307380.gff	Zilberman H2A.Z IP/total DNA 1825302 (Ab-1 met1-6)

# The rationale for considering using the met1 h2az chip data is that in met1, h2az distribution should be independent of gbM

# Ab-1 and Ab-2 are reps

#Array platform: GPL3371	FHCRC Arabidopsis Tiling Array
#GPL3371-3816.txt contains array design
#The array was designed for TAIR7, so the array probes need to be realigned to TAIR9/10 to correctly assign IP values to genomic loci 

#this needs executing first in bash to align probes to genome:
#makeblastdb -in ../../../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -dbtype nucl -out TAIR10
#blastn -db TAIR10 -query GPL3371-3816.fasta -outfmt 6 > GPL3371-3816_v_TAIR10.txt
#The best BLAST result for each probe will then be used to assign that probe's IP values to a locus

setwd("X:/Daniel-zilberman/Projects/Jay-1001_methylomes/5-analysis")

library(stringr)
library(GenomicRanges)
library(ggplot2)

probe_align_file = "../0-reference/chipseq/GPL3371-3816_v_TAIR10.txt"

probe_align = read.table(file=probe_align_file, sep="\t", header=FALSE)
colnames(probe_align) = c("query","subject","percent","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")


# here we set a flag for which version of H2AZ data we will use for the rest of the script - wt or met1
h2az_set = "wt_chipseq"
#h2az_set = "met1_chipchip"

# this flag can be set TRUE after the H2AZ data has been parsed once the first time - this is a slow process, but the results get saved to disk so can be used in subsequent runs of the script
h2az_loaded=TRUE


if (h2az_set=="met1_chipchip") {
  if (!h2az_loaded) {
    h2az_rep1 = read.table(file="../0-reference/chipseq/GSM307380.gff", sep="\t", header=FALSE)
    h2az_rep2 = read.table(file="../0-reference/chipseq/GSM307379.gff", sep="\t", header=FALSE)

    h2az_rep1$probe = str_split_fixed(str_split_fixed(h2az_rep1$V9, ';',3)[,2],"=",2)[,2]
    h2az_rep2$probe = str_split_fixed(str_split_fixed(h2az_rep2$V9, ';',3)[,2],"=",2)[,2]

    h2az_both=merge(h2az_rep1, h2az_rep2, by="probe")
    h2az_both=h2az_both[,c(1,2,5,6,7,16)]

    cor(h2az_both$V6.x, h2az_both$V6.y)
    #[1] 0.556123

    h2az_both$new_chrom = NA
    h2az_both$new_start = NA
    h2az_both$new_end = NA

    for (this_probe in 1:nrow(h2az_both)) {
      these_aligns = probe_align[probe_align$query==h2az_both[this_probe,"probe"],]
      best_align = these_aligns[these_aligns$evalue==min(these_aligns$evalue),][1,]
      h2az_both[this_probe, "new_chrom"] = best_align$subject
      h2az_both[this_probe, "new_start"] = best_align$sstart
      h2az_both[this_probe, "new_end"] = best_align$send
    }

    h2az_both$score = (h2az_both$V6.x + h2az_both$V6.y)/2

    h2az_table = h2az_both[,7:10]
    colnames(h2az_table) = c("seqnames","start","end","score")
    h2az_table$seqnames = substr(h2az_table$seqnames,4,4)
    # slow, and only needed once:
    write.table(h2az_table, file="h2az_met1.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    write.table(cbind(h2az_table$seqnames, "source"=rep("-", nrow(h2az_table)), "identifier"=rep("-", nrow(h2az_table)), h2az_table$start, h2az_table$end, h2az_table$score, "strand"=rep("-", nrow(h2az_table)), "frame"=rep("-", nrow(h2az_table)), "attribute"=rep("-", nrow(h2az_table))), file=paste0("h2az_met1.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    h2az_loaded=TRUE
  } else {
    h2az_table = read.table(file="h2az_met1.txt", sep="\t", header=TRUE)
  }
  
  h2az_table = h2az_table[!is.na(h2az_table$start),]

  # fix the probes that aligned backwards
  backwards_probes = h2az_table$end<=h2az_table$start
  backwards_starts = h2az_table[backwards_probes, "start"]
  backwards_ends = h2az_table[backwards_probes, "end"]
  h2az_table[backwards_probes, "start"] = backwards_ends
  h2az_table[backwards_probes, "end"] = backwards_starts

  h2az = makeGRangesFromDataFrame(df=h2az_table, seqnames.field="seqnames", start.field="start", end.field="end")
  h2az$score = h2az_table$score
} 

if (h2az_set=="wt_chipseq") {
  if (!h2az_loaded) {
    # This is the h2az data, which is for wt, originates in CHIP-Seq, and is at 50bp resolution
    # from Dave, probably: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002988
    h2az = readRDS(file="../0-reference/chipseq/H2AZ.devin_plos_genetics.w50.Chr.RDS")
    levels(h2az@seqnames@values) = substr(as.character(levels(h2az@seqnames@values)),4,4)
    h2az@seqinfo@seqnames = levels(h2az@seqnames@values)
    h2az_table = as.data.frame(h2az)
    # slow, and only needed once:
    write.table(h2az_table, file="h2az_wt.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    # write a version out with integer scores for visualisation in SeqMonk:
 	write.table(cbind(h2az_table, "int_score"=round(h2az_table$score*5,0)), file="h2az_wt_int.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(cbind(h2az_table$seqnames, "source"=rep("-", nrow(h2az_table)), "identifier"=rep("-", nrow(h2az_table)), h2az_table$start, h2az_table$end, h2az_table$score, "strand"=rep("-", nrow(h2az_table)), "frame"=rep("-", nrow(h2az_table)), "attribute"=rep("-", nrow(h2az_table))), file=paste0("h2az.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  } else {
    h2az_table = read.table(file="h2az_wt.txt", sep="\t", header=TRUE)
    h2az = makeGRangesFromDataFrame(df=h2az_table, seqnames.field="seqnames", start.field="start", end.field="end")
    h2az$score = h2az_table$score
  }
}



#hypotheses for why the model fits better to some segments than others:
#Is it the short segments?
#Is it the dense segments?
#Is it dense knobs on the ends?
#Is it density of methylation on the ends?
#Is it the overall density of methylation in the segment (weird undermethylated segments - why? H2AZ?)
#Is it H2aZ from the ends?
#Also plot frequency of methylation in accessions from the ends?
#Also plot simulated patterns of methylation from the ends

#Some where the model is over-methylating (model predicts mCG>.65?)
#Some where the data is weird (H2AZ?)

#window size 50-100?

#What proportion of sites are called cleanly? Is this related to how wrong the model is?



# are we doing the model or the data?  This flag allows us to process methylation states from data or model identically
#doing_model=FALSE
doing_model=TRUE


# These are several alternative file sets, depending which annotation we are using

# MAF 5% 2500 genes
annotation_file = "gene_gbm_segments_min_seg_cover_47_min_CG_3_trans.tsv"
annotation_name = "maf_5_percent"
#model_file="gene_gbm_segments_MAF_5_perc_trans_h2az_wt_chipseq_1p2_Sampled_D3D4D5_output_SimState_0.tsv"
model_file="gene_gbm_segments_MAF_5_perc_trans_h2az_wt_chipseq_1p2_AllLoci_D3D4D5_output_SimState_0.tsv"
#model_performance_file = "gene_gbm_segments_min_seg_cover_47_min_CG_3_trans_Sampled_D3D4D5_output_Sim30.tsv"
#model_performance_file = "gene_gbm_segments_MAF_5_perc_trans_h2az_wt_chipseq_1p2_Sampled_D3D4D5_output_Sim30.tsv"
#model_performance_file = "gene_gbm_segments_MAF_5_perc_trans_h2az_met1_chipchip_0_Sampled_D3D4D5_output_Sim30.tsv"
model_performance_file = "gene_gbm_segments_MAF_5_perc_trans_h2az_wt_chipseq_1p2_Sampled_D3D4D5_output_V1b_Sim30.tsv"
model_performance_file = "gene_gbm_segments_MAF_5_perc_trans_h2az_wt_chipseq_1p2_AllLoci_D3D4D5_output_Sim30.tsv"

# MAF 10% all genes
#annotation_file = "gene_gbm_segments_min_seg_cover_94_min_CG_3_trans.tsv"
#annotation_name = "maf_10_percent"
#model_performance_file = "gene_gbm_segments_min_seg_cover_94_min_CG_3_trans_output_Sim_All_v2.tsv"

# MinSpacing_30 2500 genes
#annotation_file = "AnnoTrans_MinAcc_0p005_MinCG_3_MinSpacing_30.tsv"
#annotation_name = "MinSpacing_30"
#model_file = "AnnoTrans_MinAcc_0p005_MinCG_3_MinSpacing_30_Sampled_output_SimState_0.tsv"
#model_performance_file = "AnnoTrans_MinAcc_0p005_MinCG_3_MinSpacing_30_Sampled_output_Sim30_D3D4D5.tsv"

# MAF 3 2500 genes Very permissive annotation (supported by at least 3 accessions)
#annotation_file = "gene_gbm_segments_min_acc_decile_2_min_CG_3_trans.tsv"
#annotation_name = "maf_3_accessions"
#model_performance_file = "gene_gbm_segments_min_acc_decile_2_min_CG_3_trans_h2az_1.2_Sampled_D3D4D5_output_Sim30.tsv"
#model_performance_file = "gene_gbm_segments_MAF_3_trans_h2az_wt_chipseq_1p2_Sampled_D3D4D5_output_Sim30.tsv"

this_accession_section = 4
this_accession_no = 1





CG_site_density_file = "cg.density.btNuc.gff"




# Read in the annotation, and the data/model methylation state. make genomicRanges
this_annotation = read.table(file=annotation_file, sep="\t", header=TRUE)
annotation.gr = makeGRangesFromDataFrame(this_annotation, seqnames.field="chromosome", start.field="start", end.field="end")

#this_accession = modelling_accessions[[this_accession_section]][this_accession_no]
this_accession = "SRX1664635"

this_data = read.table(file=paste0("m1001_samples_meth_status_mCG_density_",this_accession_section,"_",this_accession,"_3.tsv"), sep="\t", header=TRUE)
#this_data.gr = makeGRangesFromDataFrame(this_data, seqnames.field="Chromosome", start.field="Start", end.field="End")

this_model = read.table(file=model_file, sep="\t", header=TRUE)

if (doing_model) {
  # if we're doing the model, use model state subsequently, otherwise use the example data state
  this_data = this_model
}

this_data.gr = makeGRangesFromDataFrame(this_data, seqnames.field="Chromosome", start.field="Start", end.field="End")


# This loop was trying to work out metaplot histograms, but instead we moved to using the ends_analysis.pl method instead

#gene_gbm_segments_min_seg_cover_94_min_CG_3_trans_output_Sim_All_v2.tsv contains Amy's simuation results averaged over segments
#gene_gbm_segments_min_seg_cover_94_min_CG_3_trans_output_D7_All_v2.tsv

#target_segments = annotation.gr
#segment_CG_density_window = 50
#segment_mCG_density_window = 50
#H2AZ_density_window = 50
#step_size = 50
#end_size=1000

#segment_stats = NULL
#segment_locus_stats = NULL
#for (this_segment in 1:length(target_segments)) {
#  cat(paste0(this_segment," "))
#  this_gene_ID = this_annotation[this_segment, "gene_ID"]
#  this_chrom = substr(this_gene_ID, 3,3)
#  
#  #work out segment length
#  segment_length = target_segments[this_segment]@ranges@width
#  
#  #work out segment CG site density
#  segment_CG_site_olaps = findOverlaps(target_segments[this_segment], CG_site_ranges)
#  segment_no_CG_sites = length(unique(subjectHits(segment_CG_site_olaps)))
#  
#  #work out segment mCG density in data
#  #segment_data_site_olaps = findOverlaps(target_segments[this_segment], this_data.gr)
#  this_segment_data = this_data[this_data$gene_ID==this_gene_ID,]
#  segment_mCG_density_data = sum(this_segment_data$Status %in% c("M","2"), na.rm=TRUE)/sum(this_segment_data$Status %in% c("M","U","2","4"), na.rm=TRUE)
#  
#  #work out segment mCG density in simulation
#  #this_segment_sim = this_sim[this_sim$gene_ID==this_gene_ID,]
#  #segment_mCG_density_sim = sum(this_segment_sim$Status %in% c("M","2"), na.rm=TRUE)/sum(this_segment_sim$Status %in% c("M","U","2","4"), na.rm=TRUE)
#  segment_mCG_density_sim = 0
#  
#  #work out segment_call_rate
#  segment_call_rate = sum(this_segment_data$Status %in% c("M","U","2","4"), na.rm=TRUE)/length(this_segment_data$Status)
#  
#  # take convenience subsets of the main objects, so subsetting 1000 times doesn't take so long
#  segment_CG_site_ranges = CG_site_ranges[subjectHits(segment_CG_site_olaps)]
#  segment_data_olaps = findOverlaps(target_segments[this_segment], this_data.gr)
#  segment_h2az = h2az[subjectHits(findOverlaps(target_segments[this_segment], h2az))]
#  
#  #work out per-site stuff
#  for (this_bp in seq(0, end_size-1, step_size)) {
#    this_mid_point = target_segments[this_segment]@ranges@start + this_bp
#  
#    #work out local CG site density
#	this_CG_window = makeGRangesFromDataFrame(as.data.frame(t(c("seqnames"=this_chrom, "start"=this_mid_point - segment_CG_density_window/2, "end"= this_mid_point + segment_CG_density_window/2))), seqnames.field="seqnames", start.field="start", end.field="end")
#    segment_locus_CG_site_olaps = findOverlaps(this_CG_window, segment_CG_site_ranges)
#    segment_locus_no_CG_sites = length(unique(subjectHits(segment_locus_CG_site_olaps)))
#   
#    #work out local mCG density in data
# 	this_mCG_window = makeGRangesFromDataFrame(as.data.frame(t(c("seqnames"=this_chrom, "start"=this_mid_point - segment_mCG_density_window/2, "end"= this_mid_point + segment_mCG_density_window/2))), seqnames.field="seqnames", start.field="start", end.field="end")
#    segment_locus_mCG_site_olaps = findOverlaps(this_mCG_window, this_data.gr[subjectHits(segment_data_olaps)])
#	this_segment_locus_data = this_data[subjectHits(segment_data_olaps),][subjectHits(segment_locus_mCG_site_olaps),]
#    segment_locus_mCG_density_data = sum(this_segment_locus_data$Status %in% c("M","2"), na.rm=TRUE)/sum(this_segment_locus_data$Status %in% c("M","U","2","4"), na.rm=TRUE)
#  
#    #work out local mCG density in sim
# 	#this_mCG_window = makeGRangesFromDataFrame(c("seqnames"=this_chrom, "start"=this_mid_point - segment_mCG_density_window/2, "end"= this_mid_point + segment_mCG_density_window/2), seqnames.field="seqnames", start.field="start", end.field="end")
#    #segment_locus_CG_site_olaps = findOverlaps(this_mCG_window, this_data.gr)
#	#this_segment_locus_data = this_data[subjectHits{segment_locus_CG_site_olaps),]
#    #segment_locus_mCG_density_data = sum(this_segment_locus_data$Status %in% c("M","2"), na.rm=TRUE)/sum(this_segment_locus_data$Status %in% c("M","U","2","4"), na.rm=TRUE)
#    segment_locus_mCG_density_sim = 0
#
#	#work out local H2AZ density
#	this_H2AZ_window = makeGRangesFromDataFrame(as.data.frame(t(c("seqnames"=this_chrom, "start"=this_mid_point - H2AZ_density_window/2, "end"= this_mid_point + H2AZ_density_window/2))), seqnames.field="seqnames", start.field="start", end.field="end")
#    segment_locus_H2AZ_olaps = findOverlaps(this_H2AZ_window, segment_h2az)
#    segment_locus_H2AZ_density = mean(segment_h2az[subjectHits(segment_locus_H2AZ_olaps)]@elementMetadata@listData$score)
# 
#    #update segment_locus_stats
#    if ((this_segment==1) & (this_bp==0)) {
#      segment_locus_stats = c("gene_ID"=this_gene_ID, "locus"=this_bp, "segment_locus_no_CG_sites"=segment_locus_no_CG_sites, "segment_locus_mCG_density_data"=segment_locus_mCG_density_data, "segment_locus_mCG_density_sim"=segment_locus_mCG_density_sim, "segment_locus_H2AZ_density"=segment_locus_H2AZ_density) 
#    } else {
#      segment_locus_stats = rbind(segment_locus_stats, c("gene_ID"=this_gene_ID, "locus"=this_bp, "segment_locus_no_CG_sites"=segment_locus_no_CG_sites, "segment_locus_mCG_density_data"=segment_locus_mCG_density_data, "segment_locus_mCG_density_sim"=segment_locus_mCG_density_sim, "segment_locus_H2AZ_density"=segment_locus_H2AZ_density))
#    }
#  }
#  #update segment_stats
#  if (this_segment==1) {
#    segment_stats = c("gene_ID"=this_gene_ID, "segment_length"=segment_length, "segment_no_CG_sites"=segment_no_CG_sites, "segment_mCG_density_data"=segment_mCG_density_data, "segment_mCG_density_sim"=segment_mCG_density_sim, "segment_call_rate"=segment_call_rate) 
#  } else {
#    segment_stats = rbind(segment_stats, c("gene_ID"=this_gene_ID, "segment_length"=segment_length, "segment_no_CG_sites"=segment_no_CG_sites, "segment_mCG_density_data"=segment_mCG_density_data, "segment_mCG_density_sim"=segment_mCG_density_sim, "segment_call_rate"=segment_call_rate))
#  }
#}

#write.table(segment_stats, file="segment_stats_20bp.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(segment_locus_stats, file="segment_locus_stats_20bp.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Plot stuff based on quantiles of various things in model_performance



#ends_analysis.pl requires 9 column GFF files
#        'seqname'   => lc $seqname,
#        'source'    => $source,
#        'feature'   => $feature,
#        'start'     => $start,
#        'end'       => $end,
#        'score'     => $score,
#        'strand'    => $strand,
#        'frame'     => $frame,
#        'attribute' => $attribute


write.table(cbind(this_annotation$chromosome, "source"=rep("-", nrow(this_annotation)), this_annotation$gene_ID, this_annotation$start, this_annotation$end, "score"=rep("-", nrow(this_annotation)), "strand"=rep("-", nrow(this_annotation)), "frame"=rep("-", nrow(this_annotation)), "attribute"=rep("-", nrow(this_annotation))), file=paste0(annotation_name,"_annotation.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


# we want to work out model performance variance with CG site density, segment length, etc, so load in some density data and the model performance data

CG_site_density = read.table(file=CG_site_density_file, sep="\t", header=FALSE)
CG_site_density$V1 = substr(CG_site_density$V1,4,4)
CG_site_density_ranges = makeGRangesFromDataFrame(CG_site_density, seqnames.field="V1", start.field="V4", end.field="V5")
density_annotation_olaps = findOverlaps(CG_site_density_ranges, annotation.gr)
relevant_CG_site_density_ranges = CG_site_density_ranges[queryHits(density_annotation_olaps)]
relevant_CG_site_density = CG_site_density[queryHits(density_annotation_olaps),]
#relevant_CG_sites = CG_site_ranges[queryHits(findOverlaps(CG_site_ranges, annotation.gr))]
relevant_methylation_data = this_data[queryHits(findOverlaps(this_data.gr, annotation.gr)),]
#throw away non-M/U values
relevant_methylation_values = relevant_methylation_data[relevant_methylation_data$Status %in% c("M","U","2","4"),]
relevant_methylation_values$Status = ifelse(relevant_methylation_values$Status %in% c("M","2"),1,0)
relevant_mCG_sites = makeGRangesFromDataFrame(relevant_methylation_values, seqnames.field="Chromosome", start.field="Start", end.field="End")

model_performance = read.table(file=model_performance_file, sep="\t", header=TRUE)
model_performance$mu_meth_diff = model_performance$D3D4D5_mu_meth_diff
ggplot(model_performance) + geom_histogram(aes(x=mu_meth_diff)) + labs(title=annotation_name)

# use the histogram of model performance to set limits for what we like to call low, mid, high as in under-, mid- and over-methylation of segments
low_max = -0.3
high_min = 0.3
mid_min = -0.25
mid_max = 0.25

ggplot(model_performance) + geom_point(aes(x=CG_density, y=mu_meth_diff)) + geom_density_2d(aes(x=CG_density, y=mu_meth_diff)) + ggtitle(paste0(annotation_name, "  Rsquared=",cor(model_performance$mu_meth_diff, model_performance$CG_density)^2))
ggplot(model_performance) + geom_point(aes(x=Col0_meth_frac, y=mu_meth_diff)) + geom_density_2d(aes(x=Col0_meth_frac, y=mu_meth_diff)) + ggtitle(paste0(annotation_name, "  Rsquared=",cor(model_performance$mu_meth_diff, model_performance$Col0_meth_frac)^2))
ggplot(model_performance) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_name, "  Rsquared=",cor(model_performance$Sim_mu_methylation_fraction, model_performance$D3D4D5_mu_methylation_fraction)^2))

for (min_segment_length in 1:40) {
  #cat(paste0("Min L_CG: ",min_segment_length, "  Rsquared: ", cor(model_performance[(model_performance$mu_meth_diff<(-.025)) & (model_performance$N_CG>=min_segment_length),]$mu_meth_diff, model_performance[(model_performance$mu_meth_diff<(-.025)) & (model_performance$N_CG>=min_segment_length),]$CG_density)^2, "\n"))
  cat(paste0("Min L_CG: ",min_segment_length, "  Rsquared: ", cor(model_performance[ (model_performance$N_CG>=min_segment_length),]$Sim_mu_methylation_fraction, model_performance[ (model_performance$N_CG>=min_segment_length),]$D3D4D5_mu_methylation_fraction)^2, "\n"))
}

# write out low mid and high sets of data for CG site density, mCG density and H2A.Z levels for subsequent end_analysis.pl

low_genes = model_performance[model_performance$mu_meth_diff<low_max, "gene_ID"]
high_genes = model_performance[model_performance$mu_meth_diff>high_min, "gene_ID"]
mid_genes = model_performance[(model_performance$mu_meth_diff>mid_min) & (model_performance$mu_meth_diff<mid_max), "gene_ID"]

low_CG_sites_table = relevant_CG_site_density[queryHits(findOverlaps(relevant_CG_site_density_ranges, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% low_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]
mid_CG_sites_table = relevant_CG_site_density[queryHits(findOverlaps(relevant_CG_site_density_ranges, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% mid_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]
high_CG_sites_table = relevant_CG_site_density[queryHits(findOverlaps(relevant_CG_site_density_ranges, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% high_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]

#low_CG_sites_table = as.data.frame(low_CG_sites)
#mid_CG_sites_table = as.data.frame(mid_CG_sites)
#high_CG_sites_table = as.data.frame(high_CG_sites)

write.table(cbind(low_CG_sites_table$V1, "source"=rep("-", nrow(low_CG_sites_table)), "feature"=rep("-", nrow(low_CG_sites_table)), low_CG_sites_table$V4, low_CG_sites_table$V5, low_CG_sites_table$V6, "strand"=rep("-", nrow(low_CG_sites_table)), "frame"=rep("-", nrow(low_CG_sites_table)), "attribute"=rep("-", nrow(low_CG_sites_table))), file=paste0(annotation_name,"_low_CG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(cbind(mid_CG_sites_table$V1, "source"=rep("-", nrow(mid_CG_sites_table)), "feature"=rep("-", nrow(mid_CG_sites_table)), mid_CG_sites_table$V4, mid_CG_sites_table$V5, mid_CG_sites_table$V6, "strand"=rep("-", nrow(mid_CG_sites_table)), "frame"=rep("-", nrow(mid_CG_sites_table)), "attribute"=rep("-", nrow(mid_CG_sites_table))), file=paste0(annotation_name,"_mid_CG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(cbind(high_CG_sites_table$V1, "source"=rep("-", nrow(high_CG_sites_table)), "feature"=rep("-", nrow(high_CG_sites_table)), high_CG_sites_table$V4, high_CG_sites_table$V5, high_CG_sites_table$V6, "strand"=rep("-", nrow(high_CG_sites_table)), "frame"=rep("-", nrow(high_CG_sites_table)), "attribute"=rep("-", nrow(high_CG_sites_table))), file=paste0(annotation_name,"_high_CG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


low_mCG_sites_table = relevant_methylation_values[relevant_methylation_values$gene_ID %in% low_genes, ]
mid_mCG_sites_table = relevant_methylation_values[relevant_methylation_values$gene_ID %in% mid_genes, ]
high_mCG_sites_table = relevant_methylation_values[relevant_methylation_values$gene_ID %in% high_genes, ]

write.table(cbind(low_mCG_sites_table$Chromosome, "source"=rep("-", nrow(low_mCG_sites_table)), "feature"=rep("-", nrow(low_mCG_sites_table)), low_mCG_sites_table$Start, low_mCG_sites_table$End+1, low_mCG_sites_table$Status, "strand"=rep("-", nrow(low_mCG_sites_table)), "frame"=rep("-", nrow(low_mCG_sites_table)), "attribute"=rep("-", nrow(low_mCG_sites_table))), file=paste0(annotation_name,"_low_mCG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(cbind(mid_mCG_sites_table$Chromosome, "source"=rep("-", nrow(mid_mCG_sites_table)), "feature"=rep("-", nrow(mid_mCG_sites_table)), mid_mCG_sites_table$Start, mid_mCG_sites_table$End+1, mid_mCG_sites_table$Status, "strand"=rep("-", nrow(mid_mCG_sites_table)), "frame"=rep("-", nrow(mid_mCG_sites_table)), "attribute"=rep("-", nrow(mid_mCG_sites_table))), file=paste0(annotation_name,"_mid_mCG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(cbind(high_mCG_sites_table$Chromosome, "source"=rep("-", nrow(high_mCG_sites_table)), "feature"=rep("-", nrow(high_mCG_sites_table)), high_mCG_sites_table$Start, high_mCG_sites_table$End+1, high_mCG_sites_table$Status, "strand"=rep("-", nrow(high_mCG_sites_table)), "frame"=rep("-", nrow(high_mCG_sites_table)), "attribute"=rep("-", nrow(high_mCG_sites_table))), file=paste0(annotation_name,"_high_mCG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


low_h2az_table = h2az_table[queryHits(findOverlaps(h2az, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% low_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]
mid_h2az_table = h2az_table[queryHits(findOverlaps(h2az, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% mid_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]
high_h2az_table = h2az_table[queryHits(findOverlaps(h2az, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% high_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]

write.table(cbind(low_h2az_table$seqnames, "source"=rep("-", nrow(low_h2az_table)), "feature"=rep("-", nrow(low_h2az_table)), low_h2az_table$start, low_h2az_table$end, low_h2az_table$score, "strand"=rep("-", nrow(low_h2az_table)), "frame"=rep("-", nrow(low_h2az_table)), "attribute"=rep("-", nrow(low_h2az_table))), file=paste0(annotation_name,"_low_h2az_",h2az_set,".gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(cbind(mid_h2az_table$seqnames, "source"=rep("-", nrow(mid_h2az_table)), "feature"=rep("-", nrow(mid_h2az_table)), mid_h2az_table$start, mid_h2az_table$end, mid_h2az_table$score, "strand"=rep("-", nrow(mid_h2az_table)), "frame"=rep("-", nrow(mid_h2az_table)), "attribute"=rep("-", nrow(mid_h2az_table))), file=paste0(annotation_name,"_mid_h2az_",h2az_set,".gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(cbind(high_h2az_table$seqnames, "source"=rep("-", nrow(high_h2az_table)), "feature"=rep("-", nrow(high_h2az_table)), high_h2az_table$start, high_h2az_table$end, high_h2az_table$score, "strand"=rep("-", nrow(high_h2az_table)), "frame"=rep("-", nrow(high_h2az_table)), "attribute"=rep("-", nrow(high_h2az_table))), file=paste0(annotation_name,"_high_h2az_",h2az_set,".gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


#write.table(cbind(h2az_table$seqnames, "source"=rep("-", nrow(h2az_table)), "identifier"=rep("-", nrow(h2az_table)), h2az_table$start, h2az_table$end, h2az_table$score, "strand"=rep("-", nrow(h2az_table)), "frame"=rep("-", nrow(h2az_table)), "attribute"=rep("-", nrow(h2az_table))), file=paste0("h2az.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


# generate the relevant bash commands to do the ends analysis

cat(paste0("perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_low_CG_sites.gff > model_ends/", annotation_name, "_low_CG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_mid_CG_sites.gff > model_ends/", annotation_name, "_mid_CG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_high_CG_sites.gff > model_ends/", annotation_name, "_high_CG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_low_CG_sites.gff > model_ends/", annotation_name, "_low_CG_sites_3_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_mid_CG_sites.gff > model_ends/", annotation_name, "_mid_CG_sites_3_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_high_CG_sites.gff > model_ends/", annotation_name, "_high_CG_sites_3_",ifelse(doing_model,"model","data"),".txt\n",
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_low_mCG_sites.gff > model_ends/", annotation_name, "_low_mCG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_mid_mCG_sites.gff > model_ends/", annotation_name, "_mid_mCG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_high_mCG_sites.gff > model_ends/", annotation_name, "_high_mCG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_low_mCG_sites.gff > model_ends/", annotation_name, "_low_mCG_sites_3_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_mid_mCG_sites.gff > model_ends/", annotation_name, "_mid_mCG_sites_3_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_high_mCG_sites.gff > model_ends/", annotation_name, "_high_mCG_sites_3_",ifelse(doing_model,"model","data"),".txt\n",
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_low_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_low_h2az_",h2az_set,"_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_mid_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_mid_h2az_",h2az_set,"_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_high_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_high_h2az_",h2az_set,"_5_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_low_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_low_h2az_",h2az_set,"_3_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_mid_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_mid_h2az_",h2az_set,"_3_",ifelse(doing_model,"model","data"),".txt\n", 
"perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_high_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_high_h2az_",h2az_set,"_3_",ifelse(doing_model,"model","data"),".txt\n"
))


# repeat the whole thing but do deciles instead of low/mid/high

performance_deciles = quantile(model_performance$mu_meth_diff, prob = seq(0, 1, length = 11), type = 5)

performance_deciles
#          0%          10%          20%          30%          40%          50%          60%          70%          80%          90%         100% 
#-0.760606061 -0.213378917 -0.133237179 -0.077688889 -0.033333333  0.004954955  0.052380952  0.111111111  0.192378242  0.323782051  0.746717172 


for (this_decile in 1:10) {
  # adjust performance_deciles[1] so that >/<= works
  min_mmd = performance_deciles[this_decile] - ifelse(this_decile==1,0.000001,0)
  max_mmd = performance_deciles[this_decile+1]
  decile_genes = model_performance[(model_performance$mu_meth_diff>min_mmd) & (model_performance$mu_meth_diff<=max_mmd), "gene_ID"]

  decile_CG_sites_table = relevant_CG_site_density[queryHits(findOverlaps(relevant_CG_site_density_ranges, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% decile_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]
  write.table(cbind(decile_CG_sites_table$V1, "source"=rep("-", nrow(decile_CG_sites_table)), "feature"=rep("-", nrow(decile_CG_sites_table)), decile_CG_sites_table$V4, decile_CG_sites_table$V5, decile_CG_sites_table$V6, "strand"=rep("-", nrow(decile_CG_sites_table)), "frame"=rep("-", nrow(decile_CG_sites_table)), "attribute"=rep("-", nrow(decile_CG_sites_table))), file=paste0(annotation_name,"_decile", this_decile, "_CG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

  #decile_mCG_sites_table = relevant_methylation_values[relevant_methylation_values$gene_ID %in% decile_genes, ]
  #write.table(cbind(decile_mCG_sites_table$Chromosome, "source"=rep("-", nrow(decile_mCG_sites_table)), "feature"=rep("-", nrow(decile_mCG_sites_table)), decile_mCG_sites_table$Start, decile_mCG_sites_table$End+1, decile_mCG_sites_table$Status, "strand"=rep("-", nrow(decile_mCG_sites_table)), "frame"=rep("-", nrow(decile_mCG_sites_table)), "attribute"=rep("-", nrow(decile_mCG_sites_table))), file=paste0(annotation_name,"_decile",this_decile, "_mCG_sites.gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

  decile_h2az_table = h2az_table[queryHits(findOverlaps(h2az, makeGRangesFromDataFrame(this_annotation[this_annotation$gene_ID %in% decile_genes, ], seqnames.field="chromosome", start.field="start", end.field="end"))),]
  write.table(cbind(decile_h2az_table$seqnames, "source"=rep("-", nrow(decile_h2az_table)), "feature"=rep("-", nrow(decile_h2az_table)), decile_h2az_table$start, decile_h2az_table$end, decile_h2az_table$score, "strand"=rep("-", nrow(decile_h2az_table)), "frame"=rep("-", nrow(decile_h2az_table)), "attribute"=rep("-", nrow(decile_h2az_table))), file=paste0(annotation_name,"_decile",this_decile, "_h2az_",h2az_set,".gff"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

  cat(paste0("perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_decile",this_decile, "_CG_sites.gff > model_ends/", annotation_name, "_decile", this_decile, "_CG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
  "perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_decile", this_decile, "_CG_sites.gff > model_ends/", annotation_name, "_decile", this_decile, "_CG_sites_3_",ifelse(doing_model,"model","data"),".txt\n", 
  "perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_decile", this_decile, "_mCG_sites.gff > model_ends/", annotation_name, "_decile", this_decile, "_mCG_sites_5_",ifelse(doing_model,"model","data"),".txt\n", 
  "perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_decile", this_decile, "_mCG_sites.gff > model_ends/", annotation_name, "_decile", this_decile, "_mCG_sites_3_",ifelse(doing_model,"model","data"),".txt\n", 
  "perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -5 --zero -x ID ", annotation_name, "_decile", this_decile, "_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_decile", this_decile, "_h2az_",h2az_set,"_5_",ifelse(doing_model,"model","data"),".txt\n", 
  "perl ../../../Software/bs-sequel/ends_analysis.pl -g ",annotation_name, "_annotation.gff -b 10 -d 2000 -s 2 -3 --zero -x ID ", annotation_name, "_decile", this_decile, "_h2az_",h2az_set,".gff > model_ends/", annotation_name, "_decile", this_decile, "_h2az_",h2az_set,"_3_",ifelse(doing_model,"model","data"),".txt\n"
  ))
}

# wt CHIP-Seq: decile 9 and 10 of mu_meth_diff are way out on the H2AZ axis compared to the other deciles, and all deciles have some 5' slope. Cut off of H2AZ at 1.2 should take care of both issues
# met1 CHIP-chip: decile 9 and 10 of mu_meth_diff are somewhat out on the H2AZ axis compared to the other deciles, and all deciles have some 5' slope. Cut off of H2AZ at -0.2 should take care of both issues



# all the above was preparation to find a cutoff of H2A.Z to use to trim annotations. Now we can get on with that

### Trim annotation with H2AZ cutoff of 1.2 or 0:

library(zoo)

if (h2az_set=="wt_chipseq") {
  # wt CHIP-Seq
  h2az_threshold = 1.2
}
if (h2az_set=="met1_chipchip") {
  # met1 CHIP-chip
   h2az_threshold = 0
}

# smooth H2A.Z data over 5 adjacent windows to damp down noise
# don't really care about ends of chromosomes, as not part of any gene, so I'm going to roll the smoothing window over the ends of the chromosomes
# pad with zeros for ends of sequence
h2az$smoothed_score = c(0,0,rollmean(h2az_table$score, k=5),0,0)

#this_annotation$chromosome

updated_annotation = this_annotation

updated_annotation$start = NA
updated_annotation$end = NA

# If TRUE, this version starts with whole genes in this_anotation, and trims them from gene ends rather than from annotation ends
do_whole_gene = TRUE

for (this_segment in 1:length(annotation.gr)) {
  segment_start = ifelse(do_whole_gene, this_annotation[this_segment,"gene_start"], this_annotation[this_segment,"start"])
  segment_end = ifelse(do_whole_gene, this_annotation[this_segment,"gene_end"], this_annotation[this_segment,"end"])
  segment_h2az = h2az[queryHits(findOverlaps(h2az, annotation.gr[this_segment]))]
  #cat(paste0(this_segment, ": ", segment_h2az$smoothed_score, "\n"))

  new_start = 1
  while ((segment_h2az$smoothed_score[new_start]>h2az_threshold) & (new_start<length(segment_h2az))) {
    new_start = new_start + 1
  }
  
  new_end = length(segment_h2az)
  if (new_end>1) {
    while ((new_end>1) & (segment_h2az$smoothed_score[new_end]>h2az_threshold)) {
      new_end = new_end - 1
    }
  }
  
  if (new_start<new_end) {
    if (new_start==1) {
      updated_annotation[this_segment, "start"] = segment_start
	} else {
      updated_annotation[this_segment, "start"] = segment_h2az@ranges@start[new_start]
	}
	
	if (new_end==length(segment_h2az)) { 
      updated_annotation[this_segment, "end"] = segment_end
	} else {
      updated_annotation[this_segment, "end"] = segment_h2az@ranges@start[new_end]+segment_h2az@ranges@width[new_end]-1
    }
  }
}

write.table(updated_annotation[!is.na(updated_annotation$start), ], file=paste0(substr(annotation_file,1,str_length(annotation_file)-4),"_h2az_",h2az_set,"_",h2az_threshold,"_2.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# What is mean H2A.Z density per segment?

mean_h2az = rep(NA, nrow(updated_annotation))
max_h2az = rep(NA, nrow(updated_annotation))

for (this_segment in 1:nrow(updated_annotation)) {
  segment_start = updated_annotation[this_segment,"start"]
  segment_end = updated_annotation[this_segment,"end"]
  if ((!is.na(segment_end-segment_start)) & (segment_end-segment_start)>=1) {
    this_segment.gr = makeGRangesFromDataFrame(as.data.frame(t(c("seqnames"=updated_annotation[this_segment,"chromosome"],"start"=updated_annotation[this_segment,"start"],"end"=updated_annotation[this_segment,"end"]))), seqnames.field="seqnames", start.field="start", end.field="end") 
    mean_h2az[this_segment] = mean(h2az[queryHits(findOverlaps(h2az, this_segment.gr))]$score)
    max_h2az[this_segment] = max(h2az[queryHits(findOverlaps(h2az, this_segment.gr))]$score)
  }
}

write.table(cbind(updated_annotation, mean_h2az, max_h2az)[!is.na(updated_annotation$start), ], file=paste0(substr(annotation_file,1,str_length(annotation_file)-4),"_h2az_",h2az_set,"_",h2az_threshold,"_mean_h2az_2.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# Analyse performance by various factors

performance_analysis = merge(cbind(updated_annotation, mean_h2az, max_h2az)[!is.na(updated_annotation$start), ], model_performance[,], by="gene_ID")

# the very highly methylated segments (M>0.8) are very short and very under-methylated
ggplot(performance_analysis) + geom_point(aes(x=N_CG, y=mu_meth_diff, colour=Col0_meth_frac))

performance_analysis_set = performance_analysis[performance_analysis$Col0_meth_frac<=0.75,]

ggplot(performance_analysis_set) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_name, "  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2))

# Rsquared improves progressively as we filter M<=1 down to M<=0.75 for Col-0 values
for (col0_m_threshold in c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6)) {
  
  filt_annotation = cbind(updated_annotation, mean_h2az, max_h2az)
  filt_annotation = filt_annotation[!is.na(filt_annotation$Col0_meth_frac) & (filt_annotation$Col0_meth_frac<=col0_m_threshold),]
  performance_analysis_set = performance_analysis[!is.na(performance_analysis$Col0_meth_frac) & (performance_analysis$Col0_meth_frac<=col0_m_threshold),]
  cat(paste0(col0_m_threshold, " ", nrow(performance_analysis_set), " ", (cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2), "\n"))
  
  print(ggplot(performance_analysis_set) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_name, " Col-0 M<=",col0_m_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2)))

  #write.table(filt_annotation, file=paste0(substr(annotation_file,1,str_length(annotation_file)-4),"_h2az_",h2az_set,"_",h2az_threshold,"_",mean_h2az_threshold,"_mean_h2az.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}


# Make a new annotation where the mean H2A.Z density per segment is less than a different, more strict threshold
# Rsquared improves progressively as we use a more strict mean_h2az threshold

for (mean_h2az_threshold in c(1, 1.04, 1.08, 1.12, 1.16, 1.2)) {
  
  filt_annotation = cbind(updated_annotation, mean_h2az, max_h2az)
  filt_annotation = filt_annotation[!is.na(filt_annotation$mean_h2az) & (filt_annotation$mean_h2az<=mean_h2az_threshold),]
  performance_analysis_set = performance_analysis[!is.na(performance_analysis$mean_h2az) & (performance_analysis$mean_h2az<=mean_h2az_threshold),]
  cat(paste0(mean_h2az_threshold, " ", nrow(performance_analysis_set), " ", (cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2), "\n"))
  
  print(ggplot(performance_analysis_set) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_name, " mean_h2az<=",mean_h2az_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2)))

  write.table(filt_annotation, file=paste0(substr(annotation_file,1,str_length(annotation_file)-4),"_h2az_",h2az_set,"_",h2az_threshold,"_",mean_h2az_threshold,"_mean_h2az_2.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}


# Read in file of nn_te genes and see effect on performance
non_te_genes = read.table(file="non_te_genes_0.05.txt", sep="\t", header=FALSE)
filt_annotation = updated_annotation[updated_annotation$gene_ID %in% non_te_genes$V1,]
performance_analysis_set = performance_analysis[performance_analysis$gene_ID %in% non_te_genes$V1,]
print(ggplot(performance_analysis_set) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_name, " mean_h2az<=",mean_h2az_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2)))

non_te_genes_strict = read.table(file="non_te_genes_strict_0.05.txt", sep="\t", header=FALSE)
filt_annotation = updated_annotation[updated_annotation$gene_ID %in% non_te_genes_strict$V1,]
performance_analysis_set = performance_analysis[performance_analysis$gene_ID %in% non_te_genes_strict$V1,]
print(ggplot(performance_analysis_set) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_name, " mean_h2az<=",mean_h2az_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2)))


#Bit of both:
for (this_te_maf in c(0.01, 0.02, 0.03, 0.04, 0.05)) {
  annotation_nickname = "gbM_MAF: 0.05"
  
  non_te_genes = read.table(file=paste0("non_te_genes_", this_te_maf, ".txt"), sep="\t", header=FALSE)

   performance_analysis_set = performance_analysis[!is.na(performance_analysis$mean_h2az) & (performance_analysis$mean_h2az<=mean_h2az_threshold),]
   performance_analysis_set = performance_analysis_set[performance_analysis_set$gene_ID %in% non_te_genes$V1,]
   print(ggplot(performance_analysis_set) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_nickname, " TE_MAF: ", this_te_maf, "  mean_h2az<=",mean_h2az_threshold,"  N=", nrow(performance_analysis_set),"  Rsquared=",round(cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2, 2))))

   non_te_genes_strict = read.table(file=paste0("non_te_genes_strict_", this_te_maf, ".txt"), sep="\t", header=FALSE)

   performance_analysis_set = performance_analysis[!is.na(performance_analysis$mean_h2az) & (performance_analysis$mean_h2az<=mean_h2az_threshold),]
   performance_analysis_set = performance_analysis_set[performance_analysis_set$gene_ID %in% non_te_genes_strict$V1,]
   print(ggplot(performance_analysis_set) + geom_point(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + geom_density_2d(aes(x=D3D4D5_mu_methylation_fraction, y=Sim_mu_methylation_fraction)) + ggtitle(paste0(annotation_nickname, " TE_MAF: ", this_te_maf, " strict  mean_h2az<=",mean_h2az_threshold,"  N=", nrow(performance_analysis_set),"  Rsquared=",round(cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2, 2))))
} 

 
print(ggplot(performance_analysis_set) + geom_point(aes(x=L_locus, y=CG_density, colour=mu_meth_diff)) + ggtitle(paste0(annotation_name, " mean_h2az<=",mean_h2az_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2))+ scale_color_gradientn(colours=rainbow(7)))

print(ggplot(performance_analysis_set) + geom_point(aes(x=mean_h2az, y=Col0_meth_frac, colour=mu_meth_diff)) + ggtitle(paste0(annotation_name, " mean_h2az<=",mean_h2az_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2))+ scale_color_gradientn(colours=rainbow(7)))

print(ggplot(performance_analysis_set) + geom_point(aes(x=mean_h2az, y=CG_density, colour=mu_meth_diff)) + ggtitle(paste0(annotation_name, " mean_h2az<=",mean_h2az_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2))+ scale_color_gradientn(colours=rainbow(7))+ scale_x_log10())
 
print(ggplot(performance_analysis_set) + geom_point(aes(x=mean_h2az, y=CG_density, colour=mu_meth_diff)) + ggtitle(paste0(annotation_name, " mean_h2az<=",mean_h2az_threshold,"  Rsquared=",cor(performance_analysis_set$Sim_mu_methylation_fraction, performance_analysis_set$D3D4D5_mu_methylation_fraction)^2))+ scale_color_gradientn(colours=rainbow(7)))
 


# genes over- and under-methylated by the model
mean(performance_analysis_set$mu_meth_diff)+2*sd(performance_analysis_set$mu_meth_diff)
#[1] 0.3236266
mean(performance_analysis_set$mu_meth_diff)-2*sd(performance_analysis_set$mu_meth_diff)
#[1] -0.2664574
performance_analysis_set$performance_cat = ifelse(performance_analysis_set$mu_meth_diff<mean(performance_analysis_set$mu_meth_diff)-2*sd(performance_analysis_set$mu_meth_diff),-1,ifelse(performance_analysis_set$mu_meth_diff>mean(performance_analysis_set$mu_meth_diff)+2*sd(performance_analysis_set$mu_meth_diff),1,0))

ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=mu_meth_diff))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=CG_spacing))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=N_CG))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=mean_h2az))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=max_h2az))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=L_locus))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=CG_density))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=Sim_mu_methylation_fraction))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=Sim_sigma_methylation_fraction))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=D3D4D5_mu_methylation_fraction))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=D3D4D5_sigma_methylation_fraction))
ggplot(performance_analysis_set) + geom_boxplot(aes(group=performance_cat, y=Col0_meth_frac))


# Remove segments which are within 0, 500bp, 1kb of a TEM segment in COl-0

segmentation_model.gr = readRDS(file = "../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_segmentation_model_draft3.rds")

# Capitalise chromosome names in segmentation_model.gr
#levels(segmentation_model.gr@seqnames)=toupper(levels(segmentation_model.gr@seqnames))
#segmentation_model.gr@seqinfo@seqnames=levels(segmentation_model.gr@seqnames)

# Shorten chromosome names in segmentation_model.gr
levels(segmentation_model.gr@seqnames)=substr(levels(segmentation_model.gr@seqnames),4,4)
segmentation_model.gr@seqinfo@seqnames=levels(segmentation_model.gr@seqnames)

GBM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]
TEM_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="TEM"]
UMR_segments.gr = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="UMR"]
length(GBM_segments.gr)
length(TEM_segments.gr)
length(UMR_segments.gr)

# find segments in model overlapping or near to TEM segments
performance_analysis_set.gr = makeGRangesFromDataFrame(performance_analysis_set, seqnames.field="chromosome", start.field="start", end.field="end")

model_tem_olaps = findOverlaps(performance_analysis_set.gr, TEM_segments.gr)
new_performance_analysis_set = performance_analysis_set[-queryHits(model_tem_olaps),] 

# What if we try to fit a linear model to some of the data, and use that to predict the remaining data?
 
 
 
 

# Find genes which are overlapped by a TE in more than 5% of accessions and remove them
gene_accession_status = read.table(file="m1001_gene_gbm_tem_umr_calls_2.txt", sep="\t", header=TRUE)

for (tem_maf in c(0.01, 0.02, 0.03, 0.04, 0.05)) {
  
}




# We now make an annotation of mean H2A.Z at the closest point to each CG site, so it can be used for future modelling

# this step takes ages!
all_samples_meth_status = readRDS(file=paste0("m1001_CG_all_samples_meth_status_0.005_0.05_10_25_SNPs.rds"))

CG_site_ranges = makeGRangesFromDataFrame(all_samples_meth_status, seqnames.field="Chromosome", start.field="Locus", end.field = "Locus")
levels(CG_site_ranges@seqnames@values) = substr(as.character(levels(CG_site_ranges@seqnames@values)),4,4)
CG_site_ranges@seqinfo@seqnames = levels(CG_site_ranges@seqnames@values)

nuclear_CG_sites = CG_site_ranges[seqnames(CG_site_ranges) %in% c("1","2","3","4","5")]
nuclear_h2az=h2az[seqnames(h2az) %in% c("1","2","3","4","5")]
site_h2az = nuclear_h2az[nearest(nuclear_CG_sites, nuclear_h2az)]$score

write.table(cbind(as.data.frame(nuclear_CG_sites)[,1:2], "h2az"=site_h2az), file=paste0("nuclear_CG_site_nearest_h2az_",h2az_set,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)



# Load gff.genes and find a list of non-overlapping genes
 pathroot="X:/Daniel-Zilberman"
 
 reference_gff = paste0(pathroot,"/Reference_data/Genomes/Arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff")  # AraPort11 annotation
 
 
 library(stringr)
 
 # Load the annotation
 annot_gff = read.delim(reference_gff, header=F, comment.char="#")	
 
 gff.genes = annot_gff[annot_gff[,3]=="gene",]
 
 
 # Convert chromosome names to uppercase to match previous objects
 gff.genes$V1=toupper(gff.genes$V1)
 
 gff.genes$gene_ID=str_split_fixed(str_split_fixed(gff.genes$V9, ';',3)[,1],"=",2)[,2]
 
 gff.genes=within(gff.genes, {gene_name=str_split_fixed(str_split_fixed(V9, ';',3)[,2],"=",2)[,2]})
 
 gene_ranges=makeGRangesFromDataFrame(df = gff.genes, start.field = "V4", end.field = "V5",seqnames.field = "V1")
 # fix chromosome IDs
 levels(gene_ranges@seqnames@values) = substr(as.character(levels(gene_ranges@seqnames@values)),4,4)
 gene_ranges@seqinfo@seqnames = levels(gene_ranges@seqnames@values)

 gene_olaps = findOverlaps(gene_ranges, gene_ranges)
 overlapping_genes = unique(queryHits(gene_olaps[queryHits(gene_olaps)!=subjectHits(gene_olaps)]))
 non_overlapping_genes = gff.genes[-overlapping_genes,]
 
 write.table(non_overlapping_genes$gene_ID, file="non_overlapping_genes.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
 