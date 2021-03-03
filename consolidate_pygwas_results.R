#results_dir = "/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results"
results_dir = "X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/pygwas/results"
#project_dir = "/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis"
project_dir =  "X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/5-analysis"
#project_dir = "~/Desktop/1001_methylomes"


# This loop aims to make a matrix of all the significant SNPs in all the accessions
# That works OK for the first 2000 but might get unwieldy with more

#for (this_threshold in c(8,10,12,14)) {
for (this_threshold in c(14)) {
  combined_results = NULL
  no_genes = 0
  all_files = list.files(paste0(results_dir,"/",this_threshold))
  infiles = all_files[grep("_results_filtered.csv", all_files, fixed=T)]
  for (this_infile in infiles) {
    cat(paste0(this_infile, "\n"))
	this_gene = substr(this_infile,1,9)
	this_dataset = read.table(file=paste0(results_dir,"/", this_threshold,"/",this_infile), sep=",", header=TRUE)
    if (no_genes ==0) {
      combined_results = this_dataset[,c(1,2,3)]
    } else {
      combined_results = merge(combined_results, this_dataset[,c(1,2,3)], by=c("chromosomes", "positions"), all=TRUE)
    }                
    colnames(combined_results)[ncol(combined_results)] = this_gene
    no_genes = no_genes + 1
  }
  write.table(combined_results, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste0(this_threshold,"_genes_1-2499_combined_results.txt"))
}


# This loop aims to go through all the accessions accumulating counts of significant p-values (and accumulating sum of -logs of all p-values?)

#for (this_threshold in c(8,10,12,14)) {
for (this_threshold in c(14)) {
for (this_threshold in c(8)) {
  #combined_results = NULL
  no_genes = 0
  all_files = list.files(paste0(results_dir,"/",this_threshold))
  infiles = all_files[grep("_results_filtered.csv", all_files, fixed=T)]
  for (this_infile in infiles) {
    #if (no_genes>18473) {
      cat(paste0(this_infile, "\n"))
	  this_gene = substr(this_infile,1,9)
	  this_dataset = read.table(file=paste0(results_dir,"/", this_threshold,"/",this_infile), sep=",", header=TRUE)
      if (no_genes ==0) {
        combined_results = cbind(this_dataset[,c(1,2)], "totals" = rep(1, nrow(this_dataset)))
      } else {
        combined_results = merge(combined_results, this_dataset[,c(1,2,3)], by=c("chromosomes", "positions"), all=TRUE)
	    combined_results$totals = ifelse(is.na(combined_results$totals),0,combined_results$totals)
        combined_results$totals = combined_results$totals + ifelse(is.na(combined_results$pvals), 0, 1)
	    combined_results = combined_results[,1:3]
      }  
	  #}	  
    #colnames(combined_results)[ncol(combined_results)] = this_gene
	# If the new accession has significant p-values, increment the totals column
    no_genes = no_genes + 1
  }
  write.table(combined_results, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste0(this_threshold,"_genes_tem_sig_gene_counts.txt"))
}

# PCA of gene-accession mean mCG matrix

 
min_coverage = 8
# all samples all genes
pheno_values_gbm = read.table(file=paste0(project_dir, "/mean_mCG_per_gbm_gene_per_sample_cov_",min_coverage,".txt"), sep="\t", header=TRUE)
#pheno_values_gbm = read.table(file="Desktop/1001 methylomes/mean_mCG_per_gene_per_sample_cov_8.txt", sep="\t", header=TRUE)

# Zaigham's set of bona fide gbm (TE-like genes marked as NA), 1 sample per accession, leaf tissue chosen first then at random in Excel. 32000 genes, 948 accessions
pheno_values_gbm = read.table(file=paste0(project_dir, "/bona_fide_gbm_",min_coverage,"_reads_per_gene_cutoff_coverage_948_accessions.txt"), sep="\t", header=TRUE)
row.names(pheno_values_gbm) = pheno_values_gbm[,1]
pheno_values_gbm = pheno_values_gbm[,-1]

#might want to log this sometime (if bringing in results values as p-values for instance)                              
#log_pheno_values_gbm = cbind(pheno_values_gbm[,1],-log(pheno_values_gbm[,2:ncol(pheno_values_gbm)]))


#Need to filter out columns/rows with no variation

library(matrixStats)
gbm_mat = data.matrix(pheno_values_gbm)[,2:ncol(pheno_values_gbm)]
gbm_row_vars = rowVars(gbm_mat, na.rm=TRUE)
gbm_col_vars = colVars(gbm_mat, na.rm=TRUE)
var_pheno_values_gbm = data.matrix(pheno_values_gbm[(!is.na(gbm_row_vars)) & (gbm_row_vars>0),(!is.na(gbm_col_vars)) & (gbm_col_vars>0)])

# SVD with imputation
install.packages("bcv")
library(bcv) 
SVDimputation = impute.svd(var_pheno_values_gbm[,2:ncol(var_pheno_values_gbm)], k = 1, maxiter=100)$x

# this still did not quite get rid of every column that pca considers to have no variance, so found this extra step. problem is, na.omit gets rid of every column with an NA
#pca = prcomp(na.omit(var_pheno_values_gbm), center=TRUE, scale=TRUE)
#pca = prcomp(na.omit(var_pheno_values_gbm[,which(apply(var_pheno_values_gbm, 2, var) != 0)]), center=TRUE, scale=TRUE)
#pca = prcomp(SVDimputation, center=TRUE, scale=TRUE)
#pca1 = prcomp(t(SVDimputation), center=TRUE, scale=TRUE)


# Let's try to get rid of columns with mostly NAs (21 of)
valid_rows = rowSums(!is.na(var_pheno_values_gbm))>400
var_pheno_values_gbm = var_pheno_values_gbm[valid_rows,]

valid_cols = colSums(!is.na(var_pheno_values_gbm))>10000
var_pheno_values_gbm = var_pheno_values_gbm[,valid_cols]


#Let's try instead the imputation method, which imputes missing values first
#http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html
# 

install.packages("missMDA")
library(missMDA)

# This part took a day and still didn't finish
nb <- estim_ncpPCA(var_pheno_values_gbm,method.cv = "gcv", verbose = FALSE, scale=TRUE) # estimate the number of components from incomplete data, can use Kfold instead
#(available methods include GCV to approximate CV - gcv is the 'quick' approximation of loo x-v)

nb$ncp

var_pheno_values_gbm_imputed = imputePCA(var_pheno_values_gbm, ncp=nb$ncp)


summary(pca)

 

# factoextra for visualisation
install.packages("factoextra")
library(factoextra)

fviz_screeplot(pca, addlabels = TRUE, ylim = c(0, 50))
var = get_pca_var(pca)

fviz_screeplot(pca1, addlabels = TRUE, ylim = c(0, 50))
var1 = get_pca_var(pca1)

#visualise circular plot of contributions of each factor to components 1 and 2
fviz_pca_var(pca, col.var = "black")
fviz_pca_var(pca1, col.var = "black")

#visualise to contributing factors to each component, and their contributions
fviz_contrib(pca1, choice = "var", axes = 1, top = 100)
pheno_values_gbm[c(4901,6175,6178,4748,4744,4746,18131,26335,3643,18107),1]
#[1] "AT1G35860" "AT1G56418" "AT1G56430" "AT1G34050" "AT1G34042" "AT1G34047" "AT3G29075" "AT5G04165" "AT1G23037" "AT3G28918"
# It is interesting that in the top 6 genes, there are two entries which contain nearby genes - is there some local effect going on here?
# Is there a separate pathway devoted to mCG of this situation?

fviz_contrib(pca1, choice = "var", axes = 2, top = 100)
pheno_values_gbm[c(1546, 30354,20457,23716,29358,27373,5936,25129,28237,28126),1]
#[1] "AT1G08123" "AT5G43570" "AT3G62450"  "AT4G23690" "AT5G28820" "AT5G09565" "AT1G54320" "AT4G36070" "AT5G17540" "AT5G16490" - no. 3 is a DNA repair enzyme

fviz_contrib(pca1, choice = "var", axes = 3, top = 100)
pheno_values_gbm[c(15759,2661,4066,3009,11099,17962,20433,2135,4444,29256),1]
#[1] "AT3G09250" "AT1G14040" "AT1G27570" "AT1G17220" "AT2G21980" "AT3G27650" "AT3G62230" "AT1G09833" "AT1G31050" "AT5G27700"


fviz_contrib(pca1, choice = "var", axes = 4, top = 100)
pheno_values_gbm[c(13880,28648,31424,7019,5154,19446,4884,22937,10351,19424),1]
#[1] "AT2G46700" "AT5G21930" "AT5G52990" "AT1G65032" "AT1G47056" "AT3G53220" "AT1G35610" "AT4G16880" "AT2G14750" "AT3G53020"


pheno_combo = c(4901,6175,6178,4748,4744,4746,18131,26335,3643,18107, 1546,30354,20457,23716,29358,27373,5936,25129,28237,28126, 15759,2661,4066,3009,11099,17962,20433,2135,4444,29256, 13880,28648,31424,7019,5154,19446,4884,22937,10351,19424)

install.packages("pheatmap")
library(pheatmap)
library("RColorBrewer")

# scale each gene so that it has a mean of 0 and SD of 1, then cluster these values
# scale scales column-wise, so transpose matrix first, and remove first column containing gene nos
scaled_var_pheno_values_gbm = scale(t(var_pheno_values_gbm[,2:ncol(var_pheno_values_gbm)]))
trunc_val = 3
trunc_scaled_var_pheno_values_gbm = ifelse(scaled_var_pheno_values_gbm >trunc_val,trunc_val,ifelse(scaled_var_pheno_values_gbm <(0-trunc_val),(0-trunc_val), scaled_var_pheno_values_gbm ))

# remove the not very variable genes
variable_enough_genes = (colMaxs(trunc_scaled_var_pheno_values_gbm, na.rm=TRUE)>0.5) & (colMins(trunc_scaled_var_pheno_values_gbm, na.rm=TRUE) < -0.5)

this_heatmap = pheatmap(trunc_scaled_var_pheno_values_gbm[,variable_enough_genes][,1:1000])


SVDimputation = impute.svd(trunc_scaled_var_pheno_values_gbm, k = 1, maxiter=100)$x
rownames(SVDimputation) = rownames(trunc_scaled_var_pheno_values_gbm)
colnames(SVDimputation) = colnames(trunc_scaled_var_pheno_values_gbm)

# sort accessions by mean methylation level
accession_totals = rowSums(SVDimputation)
names(accession_totals) = rownames(SVDimputation)
accession_totals = sort(accession_totals)

head()


install.packages("dendsort")
library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(dist(SVDimputation[,variable_enough_genes]))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")


# Let's make colours for each subpopulation
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(mat)

# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(col_groups)


that_heatmap = pheatmap(t(SVDimputation[,variable_enough_genes][,]), cluster_cols = mat_cluster_cols)
this_heatmap = pheatmap(t(SVDimputation[,variable_enough_genes][,]), kmeans=9, cluster_cols = mat_cluster_cols)
#this_heatmap = pheatmap(t(SVDimputation[match(names(accession_totals), rownames(SVDimputation)),][,variable_enough_genes][,1:100]), kmeans=9)
# need to add in gene names above, so we can tell which are in which cluster
# need to add in accession numbers instead of sample IDs

write.table(this_heatmap$kmeans$centers, file="mcg_k-means_9_cluster_centres.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(this_heatmap$kmeans$cluster, file="mcg_k-means_9_cluster_genes.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=FALSE)
pdf("mcg_clustered_genes.pdf")
that_heatmap
dev.off()
pdf("mcg_k-means_9_cluster_genes.pdf")
this_heatmap
dev.off()


# What is an appropriate value for k?
install.packages("mclust")
library(mclust)
d_clust <- Mclust(as.matrix(t(SVDimputation[,variable_enough_genes][,]), G=1:15, modelNames = mclust.options("emModelNames")))
d_clust$BIC
plot(d_clust)




# pheatmap can handle NAs in data, but not NAs in the distance matrix. Need to track these down and exclude first
# make a version of the distance matrix where 0s are NAs
na_scaled_var_pheno_values_gbm = scaled_var_pheno_values_gbm
na_scaled_var_pheno_values_gbm[na_scaled_var_pheno_values_gbm==0] = NA

giveNAs = which(is.na(as.matrix(dist(na_scaled_var_pheno_values_gbm))),arr.ind=TRUE)
head(giveNAs)




#pheatmap(data.matrix(na.omit(methylation_per_gene_per_sample)[1:1010,2:ncol(methylation_per_gene_per_sample)]), color=brewer.pal(9,"Blues"))
pheatmap(data.matrix(methylation_per_gene_per_sample[as.logical((rowSums(is.na(methylation_per_gene_per_sample))-ncol(methylation_per_gene_per_sample)+1)),][1:500,2:ncol(methylation_per_gene_per_sample)]), color=brewer.pal(9,"Blues"))
  
pca = prcomp(na.omit(data.matrix(methylation_per_gene_per_sample[methylation_per_gene_per_sample$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,][as.logical((rowSums(is.na(methylation_per_gene_per_sample[methylation_per_gene_per_sample$gene_ID %in% gff.genes[gff.genes$m_class=="Gene-body Methylated",]$gene_ID,]))-ncol(methylation_per_gene_per_sample)+1)),][,2:ncol(methylation_per_gene_per_sample)])), center=TRUE, scale=TRUE)
library(scatterplot3d)
scatterplot3d(pca$x[,1:3], pch=20, color="blue")