# This code assigns methylation status to each methylated gene - gbM, TEM, or indeterminate.  The results are output in a single table for all methylated genes.
# requires gff.genes


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