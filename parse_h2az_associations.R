#setwd("X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results/")

#results_dir = "/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results"
results_dir = "X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results"
#project_dir = "/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis"
project_dir =  "X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis"



#go through all the genes filtered results.
#For each gene, find all the SNPs which are in the list of HTA-9 associated-ones that Zaigham posted
#make a matrix of SNP vs. gene_results they are found in

associated_snps = read.table(file="SNPs_associated_in_HTA9_region.txt", sep="\t", header=TRUE)


# This loop aims to go through all the accessions accumulating sum of -logs of all p-value

#Read in the previously prepared combined results set

no_genes = 0

for (this_threshold in c(8)) {
  combined_results = NULL
  no_genes = 0
  all_files = list.files(paste0(results_dir,"/",this_threshold))
  infiles = all_files[grep("_results_filtered.csv", all_files, fixed=T)]
  for (this_infile in infiles) {
    #if(no_genes>0) {
      cat(paste0(this_infile, "\n"))
	  this_gene = substr(this_infile,1,9)
	  this_dataset = read.table(file=paste0(results_dir,"/", this_threshold,"/",this_infile), sep=",", header=TRUE)
	  log_pvals = -log(this_dataset$pvals)
	  relevant_results = this_dataset[(this_dataset$chromosomes==1) & (this_dataset$positions %in% associated_snps$SNP.position),]
	  if (nrow(relevant_results)>0) {
        if (no_genes ==0) {
	      combined_results = cbind("gene_ID"=rep(this_gene, nrow(relevant_results)), "SNP_position" = relevant_results$positions)
        } else {
		  combined_results = rbind(combined_results, cbind("gene_ID"=rep(this_gene, nrow(relevant_results)), "SNP_position" = relevant_results$positions))
		}
      }                
      #colnames(combined_results)[ncol(combined_results)] = this_gene
	  # If the new accession has significant p-values, increment the totals column
    #}
    no_genes = no_genes + 1
  }
  #write.table(combined_results, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste0(this_threshold,"_genes_1-8000_log_pvals.txt"))
}

write.table(combined_results, file="HTA9_associated_SNP_associated_genes.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# These are the strongest associated genes with the nearest SNP to HTA9:
View(as.data.frame(combined_results)[as.data.frame(combined_results)$SNP_position=="19761854",])

# Add a column to performance indicating whether gene is associated with HTA9:

performance_analysis_set$HTA9_assoc = performance_analysis_set$gene_ID %in% as.data.frame(combined_results)[as.data.frame(combined_results)$SNP_position=="19761854","gene_ID"]


# We are looking to see the difference between mean h2az of associated and non-associated genes, before trimming for h2az

# Load in untrimmed annotation
annotation_file = "gene_gbm_segments_min_seg_cover_47_min_CG_3_trans.tsv"
annotation_name = "maf_5_percent"

this_annotation = read.table(file=annotation_file, sep="\t", header=TRUE)
annotation.gr = makeGRangesFromDataFrame(this_annotation, seqnames.field="chromosome", start.field="start", end.field="end")

#backup_updated_annotation = updated_annotation

# remove segments which are in the category of >5% MAF teM
updated_annotation = this_annotation[this_annotation$gene_ID %in% non_te_genes$V1,]

# this version does whole genes
#updated_annotation$start=updated_annotation$gene_start
#updated_annotation$end=updated_annotation$gene_end

# this version does all whole genes
#updated_annotation = gff.genes[,c(10,1,4,5,4,5)]
#colnames(updated_annotation) = colnames(this_annotation)
#updated_annotation$chromosome = substr(updated_annotation$chromosome,4,4)


# Score the mean_h2az of each segment in the annotation
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

write.table(cbind(updated_annotation, mean_h2az, max_h2az)[!is.na(updated_annotation$start), ], file=paste0(substr(annotation_file,1,str_length(annotation_file)-4),"_non_te_0.05.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(cbind(updated_annotation, mean_h2az, max_h2az)[!is.na(updated_annotation$start), ], file=paste0(substr(annotation_file,1,str_length(annotation_file)-4),"_non_te_0.05_whole_genes.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(cbind(updated_annotation, mean_h2az, max_h2az)[!is.na(updated_annotation$start), ], file=paste0(substr(annotation_file,1,str_length(annotation_file)-4),"_non_te_0.05_all_whole_genes.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

