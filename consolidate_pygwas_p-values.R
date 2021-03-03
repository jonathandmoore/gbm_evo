#setwd("X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results/")

#results_dir = "/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results"
results_dir = "X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results"
#project_dir = "/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis"
project_dir =  "X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis"


# This loop aims to go through all the accessions accumulating sum of -logs of all p-value

#Read in the previously prepared combined results set
#combined_results = NULL
combined_results = read.table(file=paste0(this_threshold,"_genes_1-7801_log_pvals.txt"), sep="\t", header=TRUE)
no_genes = 2498
no_genes = 5302
1326

#for (this_threshold in c(8,10,12,14)) {
for (this_threshold in c(8)) {
  #combined_results = NULL
  no_genes = 0
  all_files = list.files(paste0(results_dir,"/",this_threshold,"/full_results"))
  infiles = all_files[grep("_results.csv", all_files, fixed=T)]
  for (this_infile in infiles) {
    if(no_genes<4237) {
      cat(paste0(this_infile, "\n"))
	  this_gene = substr(this_infile,1,9)
	  this_dataset = read.table(file=paste0(results_dir,"/", this_threshold,"/full_results","/",this_infile), sep=",", header=TRUE)
	  log_pvals = -log(this_dataset$pvals)
      if (no_genes ==0) {
        combined_results = cbind(this_dataset[,c(1,2)], "totals" = log_pvals)
      } else {
        combined_results = merge(combined_results, cbind(this_dataset[,c(1,2)], "log_pvals" = log_pvals), by=c("chromosomes", "positions"), all=TRUE)
	    combined_results$totals = ifelse(is.na(combined_results$totals),0,combined_results$totals)
        combined_results$totals = combined_results$totals + ifelse(is.na(combined_results$log_pvals), 0, combined_results$log_pvals)
	    combined_results = combined_results[,1:3]
      }                
      #colnames(combined_results)[ncol(combined_results)] = this_gene
	  # If the new accession has significant p-values, increment the totals column
    }
    no_genes = no_genes + 1
  }
  write.table(combined_results, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste0(this_threshold,"_genes_1-8000_log_pvals.txt"))
}

# PCA of gene-accession mean mCG matrix

 
min_coverage = 8
pheno_values_gbm = read.table(file=paste0(project_dir, "/mean_mCG_per_gbm_gene_per_sample_cov_",min_coverage,".txt"), sep="\t", header=TRUE)
log_pheno_values_gbm = cbind(pheno_values_gbm[,1],-log(pheno_values_gbm[,2:ncol(pheno_values_gbm)]))



#Need to filter out columns/rows with no variation

library(matrixStats)
gbm_mat = data.matrix(pheno_values_gbm)[,2:ncol(pheno_values_gbm)]
gbm_row_vars = rowVars(gbm_mat, na.rm=TRUE)
gbm_col_vars = colVars(gbm_mat, na.rm=TRUE)
var_pheno_values_gbm = data.matrix(pheno_values_gbm[(!is.na(gbm_row_vars)) & (gbm_row_vars>0),(!is.na(gbm_col_vars)) & (gbm_col_vars>0)])

# this still did not quite get rid of every column that pca considers to have no variance, so found this extra step. problem is, na.omit gets rid of every column with an NA
#pca = prcomp(na.omit(var_pheno_values_gbm), center=TRUE, scale=TRUE)
pca = prcomp(na.omit(var_pheno_values_gbm[,which(apply(var_pheno_values_gbm, 2, var) != 0)]), center=TRUE, scale=TRUE)


#Let's try instead the imputation method, which imputes missing values first
#http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html
# 

install.packages("missMDA")
library(missMDA)
summary(pca)

 

# factoextra for visualisation
install.packages("factoextra")
library(factoextra)

