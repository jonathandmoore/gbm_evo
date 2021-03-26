

trans_data = read.table(file="../../../pygwas/phenotype/mCG_fitness_epiQTLs_for_trans_effects_NA_format.csv", sep=",", header=TRUE)
trans_data = read.table(file="../../../pygwas/phenotype/gbM_retained_genes_eQTLs.csv", sep=",", header=TRUE)
trans_data = read.table(file="../../../pygwas/phenotype/teM_retained_genes_eQTLs.csv", sep=",", header=TRUE)

library(MASS)

for (this_column in 2:ncol(trans_data)) {
  y = trans_data[,this_column]
  #y[is.na(y)] = 0
  y[y==0] = min(y[y>0], na.rm=TRUE)/100
  hist(y,breaks = 12)
  result = boxcox(y~1, lambda = seq(-5,5,0.5))
  mylambda = result$x[which.max(result$y)]
  mylambda
  y2 = (y^mylambda-1)/mylambda
  hist(y2)
  trans_data[,this_column] = y2
}

write.table(trans_data, file="../../../pygwas/phenotype/mCG_fitness_epiQTLs_for_trans_effects_NA_format_boxcox.csv", sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(trans_data, file="../../../pygwas/phenotype/gbM_retained_genes_eQTLs_boxcox.csv", sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(trans_data, file="../../../pygwas/phenotype/teM_retained_genes_eQTLs_boxcox.csv", sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

