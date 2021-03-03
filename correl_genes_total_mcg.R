setwd("X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-reference/supporting/correlation")

#Gene mCG leveles in each sample
gene_mcg_by_sample = read.table(file="bona_fide_gbm_8_reads_per_gene_cutoff_coverage_948_accessions.txt", sep="\t", header=TRUE)
library(stringr)
colnames(gene_mcg_by_sample) = substr(colnames(gene_mcg_by_sample), 2, str_length(colnames(gene_mcg_by_sample)))
colnames(gene_mcg_by_sample)[1] = "gene_ID"

#global mCG in each sample
genome_mcg_by_sample = read.table(file="global_gbM_mCG_and_col_gbm.txt", sep="\t", header=TRUE)

#conservation of methylation by gene
mcg_conservation_by_gene = read.table(file="methylation_conservation_of_genes_as_per_available_calls.txt", sep="\t", header=TRUE)

gene_genome_corrs = NULL
for (this_gene in gene_mcg_by_sample$gene_ID) {
  this_gene_mcg = cbind("accessionid"=as.numeric(colnames(gene_mcg_by_sample)[2:length(colnames(gene_mcg_by_sample))]), "gene_mcg"=t(gene_mcg_by_sample[gene_mcg_by_sample$gene_ID==this_gene,2:ncol(gene_mcg_by_sample)]))
  this_gene_mcg = merge(this_gene_mcg, genome_mcg_by_sample[,1:2], by="accessionid")
  if (!is.na(sd(this_gene_mcg[,2], na.rm=TRUE))) {
    gene_genome_corrs = rbind(gene_genome_corrs, c("gene"=this_gene, "corr"=cor(this_gene_mcg[,2], this_gene_mcg[,3], use="complete.obs", method="pearson")))
  }
}

gene_genome_corrs = as.data.frame(gene_genome_corrs)
gene_genome_corrs$corr = as.numeric(gene_genome_corrs$corr)
library(ggplot2)
ggplot(as.data.frame(gene_genome_corrs)) + geom_density(aes(x=corr))

 write.table(gene_genome_corrs, file="gene_genome_mcg_corrs.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
 
 colnames(gene_genome_corrs)[1]="gene_ID"
 z = merge(mcg_conservation_by_gene, gene_genome_corrs, by="gene_ID")
 