setwd("X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-reference/supporting/haplo")

#Zaigham's list of gene regions of interest
gene_regions_of_interest = read.table(file="Retained_gbM_and_tem_eQTLs_for_haplogrouping.txt", sep="\t", header=TRUE)

#Zaigham's list of mCG gbM retained eQTL
mCG_gbm_retained_eqtl = read.table(file="mCG_gbm_retained_eqtl_tem_625_accessions.txt", sep="\t", header=TRUE)


#Umit's full 1001 genomes imputed SNP set
#BiocManager::install("zlib")
#BiocManager::install("rhdf5")
#library(zlibbioc)
library(rhdf5)
#this one is the lzh-encoded one - only works in python:
#h5f = H5Fopen("5/all_chromosomes_binary.hdf5")
h5f = H5Fopen("imputed_snps_binary_gzip.hdf5")
h5f
h5f$positions
h5f$accessions
h5f$snps

#find chromosome boundaries in snp positions
prev_position = 99999999999
for (this_position in 1:length(h5f$positions)) {
  if (h5f$positions[this_position]<prev_position) {
    cat(paste(prev_position, this_position, h5f$positions[this_position], "\n"))
  }
  prev_position = h5f$positions[this_position]
}

#find chromosome boundaries in snp positions 10k chunks
prev_position = 99999999999
for (this_position in seq(1, length(h5f$positions), by=10000)) {
     if (h5f$positions[this_position]<prev_position) {
         cat(paste(prev_position, this_position, h5f$positions[this_position], "\n"))
     }
     prev_position = h5f$positions[this_position]
 }
#99999999999 1 55 
#30287457 2600001 64192 
#19556477 4470001 56145 
#23441023 6670001 125162 
#18465705 8430001 28634 

prev_position = 99999999999
for (this_position in seq(2590001, 2600001, by=1)) {
     if (h5f$positions[this_position]<prev_position) {
         cat(paste(prev_position, this_position, h5f$positions[this_position], "\n"))
     }
     prev_position = h5f$positions[this_position]
 }
#30427620 2597736 1003 

prev_position = 99999999999
for (this_position in seq(4460001, 4470001, by=1)) {
     if (h5f$positions[this_position]<prev_position) {
         cat(paste(prev_position, this_position, h5f$positions[this_position], "\n"))
     }
     prev_position = h5f$positions[this_position]
 }
#19697934 4466531 115 

prev_position = 99999999999
for (this_position in seq(6660001, 6670001, by=1)) {
     if (h5f$positions[this_position]<prev_position) {
         cat(paste(prev_position, this_position, h5f$positions[this_position], "\n"))
     }
     prev_position = h5f$positions[this_position]
 }
#23459766 6660783 1005 

prev_position = 99999999999
for (this_position in seq(8420001, 8430001, by=1)) {
     if (h5f$positions[this_position]<prev_position) {
         cat(paste(prev_position, this_position, h5f$positions[this_position], "\n"))
     }
     prev_position = h5f$positions[this_position]
 }
#18585006 8427787 53
 
chrom_start_snps = c(1,2597736,4466531,6660783,8427787)
chrom_end_snps = c(2597735,4466530,6660782,8427786,10709466)


#write out a file for each gene
window_size = 4000
for (this_gene in 1:nrow(gene_regions_of_interest)) {
  this_chromosome = as.numeric(substr(gene_regions_of_interest[this_gene, "Chromosome"],4,4))
  positions_of_interest = h5f$positions[chrom_start_snps[this_chromosome]:chrom_end_snps[this_chromosome]]

  # make TRUE/FALSE list of accessions for inclusion
  accessions_of_interest = h5f$accessions %in% mCG_gbm_retained_eqtl$accessiondi
  accession_names = h5f$accessions[accessions_of_interest]
  
  #get out all the snps in the chromosome of interest, find the interesting region within that, and get the relevant subset of snps
  snps_of_chromosome = h5f$snps[accessions_of_interest,chrom_start_snps[this_chromosome]:chrom_end_snps[this_chromosome]]
  this_region = positions_of_interest[(positions_of_interest>=(gene_regions_of_interest[this_gene, "Start"]-window_size)) & (positions_of_interest<=(gene_regions_of_interest[this_gene, "End"]+ window_size))]
  region_snps = snps_of_chromosome[,(positions_of_interest>=(gene_regions_of_interest[this_gene, "Start"]-window_size)) & (positions_of_interest<=(gene_regions_of_interest[this_gene, "End"]+ window_size))]

  row.names(region_snps) = accession_names
  
  region_snps = cbind(region_snps, "combined" = rep(0, nrow(region_snps)))
  for (this_column in 1:ncol(region_snps)-1) {
     region_snps[,ncol(region_snps)]=paste0(region_snps[,ncol(region_snps)], region_snps[,this_column])
  }
  write.table(sort(region_snps[,ncol(region_snps)]), file=paste0(gene_regions_of_interest[this_gene,"Gene_ID"],"_snps.txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
} 

# read in the file for each gene, strip out rare SNPs, and regenerate the sorted files
rarity_threshold = 0.05

for (this_gene in 1:nrow(gene_regions_of_interest)) {
  region_snps = read.table(file=paste0(gene_regions_of_interest[this_gene,"Gene_ID"],"_snps.txt"), sep="\t", header=FALSE, colClasses = c("numeric","character"))
  for (this_column in 1:nchar(region_snps[1,2])) {
     region_snps=cbind(region_snps, substr(region_snps[,2],this_column,this_column))
	 region_snps[,2+this_column] = as.numeric(region_snps[,2+this_column])
  }
  row.names(region_snps) = region_snps[,1]
  region_snps = region_snps[,3:ncol(region_snps)]

  unrare_snps=colSums(region_snps)>=rarity_threshold*nrow(region_snps)
  
  region_snps = region_snps[,unrare_snps]
  region_snps = cbind(region_snps, "combined" = rep(0, nrow(region_snps)))
  for (this_column in 1:ncol(region_snps)-1) {
     region_snps[,ncol(region_snps)]=paste0(region_snps[,ncol(region_snps)], region_snps[,this_column])
  }
  final_patterns = region_snps[,ncol(region_snps)]
  names(final_patterns) = rownames(region_snps)
  write.table(sort(final_patterns), file=paste0(gene_regions_of_interest[this_gene,"Gene_ID"],"_snps_filtered.txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
}

#sample genes:
#AT1G09420, AT1G09910, AT4G21326, and AT5G10140 


# read in the file for each gene, figure out the haplotype groups
big_table_1 = matrix(ncol=nrow(region_snps), nrow=nrow(gene_regions_of_interest))
rownames(big_table_1) = gene_regions_of_interest$Gene_ID
colnames(big_table_1) = accession_names

min_hap_group_size = 15
max_group_size = 0
max_group_count = 0
for (this_gene in 1:nrow(gene_regions_of_interest)) {
  region_snps = read.table(file=paste0(gene_regions_of_interest[this_gene,"Gene_ID"],"_snps_filtered.txt"), sep="\t", header=FALSE, colClasses = c("numeric","character"))
  this_group = 0
  this_group_size = 0
  valid_hap_groups = NULL
  prev_pattern = "99"
  for (this_accession in 1:nrow(region_snps)) {
    if (region_snps[this_accession,2]==prev_pattern) {
	  this_group_size = this_group_size + 1
	  if (this_group_size>max_group_size) {
	    max_group_size = this_group_size
	  }
	} else {
	  if (this_group_size>=min_hap_group_size) {
	    valid_hap_groups = c(valid_hap_groups, this_group)
	  }
	  this_group_size = 0
	  this_group = this_group + 1
	  prev_pattern=region_snps[this_accession,2]
	}
    big_table_1[rownames(big_table_1)==gene_regions_of_interest[this_gene,"Gene_ID"], colnames(big_table_1)==region_snps[this_accession,1]] = this_group
  }
  # deal with the last group we were doing
  if (this_group_size>=min_hap_group_size) {
   valid_hap_groups = c(valid_hap_groups, this_group)
  }
  if (max_group_count<this_group) {
    max_group_count = this_group
  }
  cat(paste(this_gene, valid_hap_groups, "\n"))
}
write.table(big_table_1, file="haplogroups.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)

#max_group_size=268
#max_group_count=622

accession_gene_gbm = read.table(file="mCG_gbm_retained_eqtl_gbm_625_accessions.txt", sep="\t", header=TRUE)
accession_gene_expression_gbm = read.table(file="mRNA_retained_eqtl_gbm_625_accessions.txt", sep="\t", header=TRUE)
accession_gene_tem = read.table(file="mCG_gbm_retained_eqtl_tem_625_accessions.txt", sep="\t", header=TRUE)
accession_gene_expression_tem = read.table(file="mRNA_retained_eqtl_tem_625_accessions.txt", sep="\t", header=TRUE)

# make some empty tables for numbers of accessions, R and P, and no_valid_obs
big_table_2a = matrix(ncol=max_group_count, nrow=nrow(gene_regions_of_interest))
rownames(big_table_2a) = gene_regions_of_interest$Gene_ID
colnames(big_table_2a) = seq(1, max_group_count, by=1)

# ones for gbm
big_table_2b = matrix(ncol=max_group_count, nrow=nrow(gene_regions_of_interest))
rownames(big_table_2b) = gene_regions_of_interest$Gene_ID
colnames(big_table_2b) = seq(1, max_group_count, by=1)

big_table_2c = matrix(ncol=max_group_count, nrow=nrow(gene_regions_of_interest))
rownames(big_table_2c) = gene_regions_of_interest$Gene_ID
colnames(big_table_2c) = seq(1, max_group_count, by=1)

big_table_2d= matrix(ncol=max_group_count, nrow=nrow(gene_regions_of_interest))
rownames(big_table_2d) = gene_regions_of_interest$Gene_ID
colnames(big_table_2d) = seq(1, max_group_count, by=1)

# ones for tem
big_table_2e = matrix(ncol=max_group_count, nrow=nrow(gene_regions_of_interest))
rownames(big_table_2e) = gene_regions_of_interest$Gene_ID
colnames(big_table_2e) = seq(1, max_group_count, by=1)

big_table_2f = matrix(ncol=max_group_count, nrow=nrow(gene_regions_of_interest))
rownames(big_table_2f) = gene_regions_of_interest$Gene_ID
colnames(big_table_2f) = seq(1, max_group_count, by=1)

big_table_2g= matrix(ncol=max_group_count, nrow=nrow(gene_regions_of_interest))
rownames(big_table_2g) = gene_regions_of_interest$Gene_ID
colnames(big_table_2g) = seq(1, max_group_count, by=1)

# read in haplogroups table, for each gene find all valid haplogroups and do correlations 
min_hap_group_size = 15
min_obs_pairs = 10
for (this_gene in 1:nrow(gene_regions_of_interest)) {
  for (this_group in 1:max(big_table_1[rownames(big_table_1)==gene_regions_of_interest[this_gene,"Gene_ID"],])) {
    this_gene_id = gene_regions_of_interest[this_gene,"Gene_ID"]
    num_members = sum(big_table_1[rownames(big_table_1)==this_gene_id,]==this_group)
	big_table_2a[rownames(big_table_2a)==this_gene_id, colnames(big_table_2a)==this_group] = num_members
	
	if (num_members>=min_hap_group_size) {
	
      relevant_accessions = colnames(big_table_1)[big_table_1[rownames(big_table_1)==this_gene_id,]==this_group]
	
	  relevant_gbm = accession_gene_gbm[accession_gene_gbm$accessionid %in% relevant_accessions, colnames(accession_gene_gbm)==this_gene_id]
  	  relevant_expression_gbm = accession_gene_expression_gbm[accession_gene_expression_gbm$accessionid %in% relevant_accessions, colnames(accession_gene_expression_gbm)==this_gene_id]
	  
	  good_overall = !(is.na(relevant_gbm) | is.na(relevant_expression_gbm))
	  
	  if (sum(good_overall)>=min_obs_pairs) {
	
  	    test_result = cor.test(relevant_gbm, relevant_expression_gbm, use="complete.obs", method="pearson")
	    big_table_2b[rownames(big_table_2a)==this_gene_id, colnames(big_table_2a)==this_group] = test_result$p.value
	    big_table_2c[rownames(big_table_2a)==this_gene_id, colnames(big_table_2a)==this_group] = test_result$estimate
	  }
      big_table_2d[rownames(big_table_2d)==this_gene_id, colnames(big_table_2d)==this_group] = sum(good_overall)
	  
	  relevant_tem = accession_gene_tem[accession_gene_tem$accessiondi %in% relevant_accessions, colnames(accession_gene_tem)==this_gene_id]
  	  relevant_expression_tem = accession_gene_expression_tem[accession_gene_expression_tem$accessionid %in% relevant_accessions, colnames(accession_gene_expression_tem)==this_gene_id]
	  
	  if ((length(relevant_expression_tem)>0) & (length(relevant_tem) > 0)) {
  	    good_overall = !(is.na(relevant_tem) | is.na(relevant_expression_tem))
	  
	    if (sum(good_overall)>=min_obs_pairs) {
	
  	      test_result = cor.test(relevant_tem, relevant_expression_tem, use="complete.obs", method="pearson")
	      big_table_2e[rownames(big_table_2a)==this_gene_id, colnames(big_table_2a)==this_group] = test_result$p.value
	      big_table_2f[rownames(big_table_2a)==this_gene_id, colnames(big_table_2a)==this_group] = test_result$estimate
	    }
        big_table_2g[rownames(big_table_2g)==this_gene_id, colnames(big_table_2g)==this_group] = sum(good_overall)
      }
    }
  }
}  
write.table(big_table_2a, file="big_table_2a.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(big_table_2b, file="big_table_2b.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(big_table_2c, file="big_table_2c.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(big_table_2d, file="big_table_2d.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(big_table_2e, file="big_table_2e.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(big_table_2f, file="big_table_2f.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(big_table_2g, file="big_table_2g.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)


#Correct the p-values for multiple tests per gene
big_table_2b_corrected = big_table_2b
big_table_2e_corrected = big_table_2e

for (this_gene in 1:nrow(gene_regions_of_interest)) {
  these_tests = !is.na(big_table_2b[this_gene,])
  no_tests = sum(these_tests)
  if (no_tests>0) {
    big_table_2b_corrected[this_gene, these_tests] = p.adjust(big_table_2b[this_gene, these_tests], method="fdr", n = no_tests)
  }
  cat(paste0(min(big_table_2b[this_gene, these_tests]), " ", min(big_table_2b_corrected[this_gene, these_tests]), "\n"))

  these_tests = !is.na(big_table_2e[this_gene,])
  no_tests = sum(these_tests)
  if (no_tests>0) {
    big_table_2e_corrected[this_gene, these_tests] = p.adjust(big_table_2e[this_gene, these_tests], method="fdr", n = no_tests)
  }
}
write.table(big_table_2b_corrected, file="big_table_2b_corrected.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
write.table(big_table_2e_corrected, file="big_table_2e_corrected.txt", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)

# Extract the SNPs for a single locus, as a test
this_gene=98 #FLC At5g10140
region_snps = read.table(file=paste0(gene_regions_of_interest[this_gene,"Gene_ID"],"_snps.txt"), sep="\t", header=FALSE, colClasses = c("numeric","character"))
this_chromosome = as.numeric(substr(gene_regions_of_interest[this_gene, "Chromosome"],4,4))
positions_of_interest = h5f$positions[chrom_start_snps[this_chromosome]:chrom_end_snps[this_chromosome]]
# make TRUE/FALSE list of accessions for inclusion
accessions_of_interest = h5f$accessions %in% mCG_gbm_retained_eqtl$accessiondi
accession_names = h5f$accessions[accessions_of_interest]
#get out all the snps in the chromosome of interest, find the interesting region within that, and get the relevant subset of snps
snps_of_chromosome = h5f$snps[accessions_of_interest,chrom_start_snps[this_chromosome]:chrom_end_snps[this_chromosome]]
this_region = positions_of_interest[(positions_of_interest>=(gene_regions_of_interest[this_gene, "Start"]-window_size)) & (positions_of_interest<=(gene_regions_of_interest[this_gene, "End"]+ window_size))]
region_snps = snps_of_chromosome[,(positions_of_interest>=(gene_regions_of_interest[this_gene, "Start"]-window_size)) & (positions_of_interest<=(gene_regions_of_interest[this_gene, "End"]+ window_size))]

View(region_snps)
this_region
length(this_region)
[1] 930
ncol(region_snps)
[1] 930
colnames(region_snps) = this_region
View(region_snps)
unrare_snps=colSums(region_snps)>=rarity_threshold*nrow(region_snps)
 
sum(unrare_snps)
[1] 173
region_snps = region_snps[,unrare_snps]

row.names(region_snps) = accession_names
View(region_snps)
region_snps = cbind(region_snps, "combined" = rep(0, nrow(region_snps)))
for (this_column in 1:ncol(region_snps)-1) {
  region_snps[,ncol(region_snps)]=paste0(region_snps[,ncol(region_snps)], region_snps[,this_column])
}
 
region_snps[1:10,ncol(region_snps)]

write.table(region_snps, file="FLC_patterns.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


# Define haplogroups for some flowering and fitness genes
# AT5G10140 Chr5	3173382	3179448
# AT5G49440 Chr5	20048239	20049776
# AT2G39890 Chr2	16655577	16658412
# AT1G19410 Chr1	6714492	6716439

fandf_genes = c("AT5G10140", "AT5G49440", "AT2G39890", "AT1G19410")
fandf_chroms = c(5,5,2,1)
fandf_starts = c(3173382, 20048239, 16655577, 6714492)
fandf_ends = c(3179448, 20049776, 16658412, 6716439)

window_size = 4000
for (this_gene in 1:length(fandf_genes)) {
  this_chromosome = fandf_chroms[this_gene]
  positions_of_interest = h5f$positions[chrom_start_snps[this_chromosome]:chrom_end_snps[this_chromosome]]

  # make TRUE/FALSE list of accessions for inclusion
  #accessions_of_interest = h5f$accessions %in% mCG_gbm_retained_eqtl$accessiondi
  #accession_names = h5f$accessions[accessions_of_interest]
  
  accessions_of_interest = h5f$accessions %in% read.table(file=paste0("haplogroups_flowering_and_fitness_genes_", fandf_genes[this_gene],".txt"), sep="\t", header=TRUE, colClasses="character")[,1]
  accession_names = h5f$accessions[accessions_of_interest]
  
  #get out all the snps in the chromosome of interest, find the interesting region within that, and get the relevant subset of snps
  snps_of_chromosome = h5f$snps[accessions_of_interest,chrom_start_snps[this_chromosome]:chrom_end_snps[this_chromosome]]
  this_region = positions_of_interest[(positions_of_interest>=(fandf_starts[this_gene]-window_size)) & (positions_of_interest<=(fandf_ends[this_gene]+window_size))]
  region_snps = snps_of_chromosome[,(positions_of_interest>=(fandf_starts[this_gene]-window_size)) & (positions_of_interest<=(fandf_ends[this_gene]+ window_size))]
  rm(positions_of_interest, snps_of_chromosome)
  
  row.names(region_snps) = accession_names
  
  region_snps = cbind(region_snps, "combined" = rep(0, nrow(region_snps)))
  for (this_column in 1:ncol(region_snps)-1) {
     region_snps[,ncol(region_snps)]=paste0(region_snps[,ncol(region_snps)], region_snps[,this_column])
  }
  write.table(sort(region_snps[,ncol(region_snps)]), file=paste0(fandf_genes[this_gene],"_fandf_snps.txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
  rm(region_snps)
} 


# read in the file for each gene, strip out rare SNPs, and regenerate the sorted files
rarity_threshold = 0.05

for (this_gene in 1:length(fandf_genes)) {
  region_snps = read.table(file=paste0(fandf_genes[this_gene],"_fandf_snps.txt"), sep="\t", header=FALSE, colClasses = c("numeric","character"))
  for (this_column in 1:nchar(region_snps[1,2])) {
     region_snps=cbind(region_snps, substr(region_snps[,2],this_column,this_column))
	 region_snps[,2+this_column] = as.numeric(region_snps[,2+this_column])
  }
  row.names(region_snps) = region_snps[,1]
  region_snps = region_snps[,3:ncol(region_snps)]

  unrare_snps=colSums(region_snps)>=rarity_threshold*nrow(region_snps)
  
  region_snps = region_snps[,unrare_snps]
  region_snps = cbind(region_snps, "combined" = rep(0, nrow(region_snps)))
  for (this_column in 1:ncol(region_snps)-1) {
     region_snps[,ncol(region_snps)]=paste0(region_snps[,ncol(region_snps)], region_snps[,this_column])
  }
  final_patterns = region_snps[,ncol(region_snps)]
  names(final_patterns) = rownames(region_snps)
  write.table(sort(final_patterns), file=paste0(fandf_genes[this_gene],"_fandf_snps_filtered.txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
}

# read in the file for each fandf gene, figure out the haplotype groups
#big_table_1_fandf = matrix(ncol=nrow(region_snps), nrow=nrow(gene_regions_of_interest))
#rownames(big_table_1) = gene_regions_of_interest$Gene_ID
#colnames(big_table_1) = accession_names

min_hap_group_size = 15
max_group_size = 0
max_group_count = 0
for (this_gene in 1:length(fandf_genes)) {
  region_snps = read.table(file=paste0(fandf_genes[this_gene],"_fandf_snps_filtered.txt"), sep="\t", header=FALSE, colClasses = c("numeric","character"))
  big_table_1_fandf = matrix(ncol=nrow(region_snps), nrow=1)
  rownames(big_table_1_fandf) = fandf_genes[this_gene]
  colnames(big_table_1_fandf) = region_snps[,1]

  this_group = 0
  this_group_size = 0
  valid_hap_groups = NULL
  prev_pattern = "99"
  for (this_accession in 1:nrow(region_snps)) {
    if (region_snps[this_accession,2]==prev_pattern) {
	  this_group_size = this_group_size + 1
	  if (this_group_size>max_group_size) {
	    max_group_size = this_group_size
	  }
	} else {
	  if (this_group_size>=min_hap_group_size) {
	    valid_hap_groups = c(valid_hap_groups, this_group)
	  }
	  this_group_size = 0
	  this_group = this_group + 1
	  prev_pattern=region_snps[this_accession,2]
	}
    big_table_1_fandf[rownames(big_table_1_fandf)==fandf_genes[this_gene], colnames(big_table_1_fandf)==region_snps[this_accession,1]] = this_group
  }
  # deal with the last group we were doing
  if (this_group_size>=min_hap_group_size) {
   valid_hap_groups = c(valid_hap_groups, this_group)
  }
  if (max_group_count<this_group) {
    max_group_count = this_group
  }
  cat(paste(this_gene, valid_hap_groups, "\n"))
  write.table(big_table_1_fandf, file=paste0("haplogroups_",fandf_genes[this_gene],".txt"), sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
}


# extract the mean, standard deviation of mCG and expression, R values and mCG range for lowest p-value haplogroups
cat(paste("gene_id", "mean_mcg_gbm", "stdev_mcg_gbm", "mean_exp_gbm", "stdev_exp_gbm", "r_value_gbm", "mcg_min_gbm", "mcg_max_gbm", "\n", sep="\t"))
#, "mean_mcg_tem", "stdev_mcg_tem", "mean_exp_tem", "stdev_exp_tem", "r_value_tem", "mcg_min_tem", "mcg_max_tem",
for (this_gene in 1:nrow(gene_regions_of_interest)) {
  this_gene_id = gene_regions_of_interest$Gene_ID[this_gene]
  these_pvalues = big_table_2b[rownames(big_table_2b)==this_gene_id,]
  best_pvalue = min(these_pvalues, na.rm=TRUE)
  best_group = colnames(big_table_2b)[these_pvalues==best_pvalue][!is.na(colnames(big_table_2b)[these_pvalues==best_pvalue])]
  if (length(best_group)==1) {
    relevant_accessions = as.numeric(colnames(big_table_1)[big_table_1[rownames(big_table_1)==this_gene_id,]==as.numeric(best_group)])
	
    relevant_gbm = accession_gene_gbm[accession_gene_gbm$accessionid %in% relevant_accessions, colnames(accession_gene_gbm)==this_gene_id]
    relevant_expression_gbm = accession_gene_expression_gbm[accession_gene_expression_gbm$accessionid %in% relevant_accessions, colnames(accession_gene_expression_gbm)==this_gene_id]
    #relevant_tem = accession_gene_tem[accession_gene_tem$accessiondi %in% relevant_accessions, colnames(accession_gene_tem)==this_gene_id]
    #relevant_expression_tem = accession_gene_expression_tem[accession_gene_expression_tem$accessionid %in% relevant_accessions, colnames(accession_gene_expression_tem)==this_gene_id]

	
	mean_mcg_gbm = mean(relevant_gbm, na.rm=TRUE)
	stdev_mcg_gbm = sd(relevant_gbm, na.rm=TRUE)
	mean_exp_gbm = mean(relevant_expression_gbm, na.rm=TRUE)
	stdev_exp_gbm = sd(relevant_expression_gbm, na.rm=TRUE)
	
	r_value_gbm = cor(relevant_gbm, relevant_expression_gbm, method="pearson", use="pairwise.complete.obs")
	mcg_min_gbm = min(relevant_gbm, na.rm=TRUE)
	mcg_max_gbm = max(relevant_gbm, na.rm=TRUE)
	#mean_mcg_tem = mean(relevant_tem, na.rm=TRUE)
	#stdev_mcg_tem= sd(relevant_tem, na.rm=TRUE)
	#mean_exp_tem = mean(relevant_expression_tem, na.rm=TRUE)
	#stdev_exp_tem = sd(relevant_expression_tem, na.rm=TRUE)
	#r_value_tem = cor(relevant_tem, relevant_expression_tem, method="pearson", use="pairwise.complete.obs")
	#mcg_min_tem = min(relevant_tem, na.rm=TRUE)
	#mcg_max_tem = max(relevant_tem, na.rm=TRUE)
	cat(paste(this_gene_id, mean_mcg_gbm, stdev_mcg_gbm, mean_exp_gbm, stdev_exp_gbm, r_value_gbm, mcg_min_gbm, mcg_max_gbm, "\n", sep="\t"))
#, mean_mcg_tem, stdev_mcg_tem, mean_exp_tem, stdev_exp_tem, r_value_tem, mcg_min_tem, mcg_max_tem
  } else {
    cat(paste0(this_gene_id," has ", length(best_group), " best haplogroups.\n"))
  }
}

# extract the mean, standard deviation of mCG and expression, R values and mCG range for lowest p-value haplogroups
cat(paste("gene_id", "mean_mcg_tem", "stdev_mcg_tem", "mean_exp_tem", "stdev_exp_tem", "r_value_tem", "mcg_min_tem", "mcg_max_tem", "\n", sep="\t"))
#, "mean_mcg_gbm", "stdev_mcg_gbm", "mean_exp_gbm", "stdev_exp_gbm", "r_value_gbm", "mcg_min_gbm", "mcg_max_gbm",
for (this_gene in 1:nrow(gene_regions_of_interest)) {
  this_gene_id = gene_regions_of_interest$Gene_ID[this_gene]
  these_pvalues = big_table_2e[rownames(big_table_2e)==this_gene_id,]
  best_pvalue = min(these_pvalues, na.rm=TRUE)
  best_group = colnames(big_table_2e)[these_pvalues==best_pvalue][!is.na(colnames(big_table_2e)[these_pvalues==best_pvalue])]
  if (length(best_group)==1) {
    relevant_accessions = as.numeric(colnames(big_table_1)[big_table_1[rownames(big_table_1)==this_gene_id,]==as.numeric(best_group)])
	
    #relevant_gbm = accession_gene_gbm[accession_gene_gbm$accessionid %in% relevant_accessions, colnames(accession_gene_gbm)==this_gene_id]
    #relevant_expression_gbm = accession_gene_expression_gbm[accession_gene_expression_gbm$accessionid %in% relevant_accessions, colnames(accession_gene_expression_gbm)==this_gene_id]
    relevant_tem = accession_gene_tem[accession_gene_tem$accessiondi %in% relevant_accessions, colnames(accession_gene_tem)==this_gene_id]
    relevant_expression_tem = accession_gene_expression_tem[accession_gene_expression_tem$accessionid %in% relevant_accessions, colnames(accession_gene_expression_tem)==this_gene_id]

	
	#mean_mcg_gbm = mean(relevant_gbm, na.rm=TRUE)
	#stdev_mcg_gbm = sd(relevant_gbm, na.rm=TRUE)
	#mean_exp_gbm = mean(relevant_expression_gbm, na.rm=TRUE)
	#stdev_exp_gbm = sd(relevant_expression_gbm, na.rm=TRUE)
	
	#r_value_gbm = cor(relevant_gbm, relevant_expression_gbm, method="pearson", use="pairwise.complete.obs")
	#mcg_min_gbm = min(relevant_gbm, na.rm=TRUE)
	#mcg_max_gbm = max(relevant_gbm, na.rm=TRUE)
	mean_mcg_tem = mean(relevant_tem, na.rm=TRUE)
	stdev_mcg_tem= sd(relevant_tem, na.rm=TRUE)
	mean_exp_tem = mean(relevant_expression_tem, na.rm=TRUE)
	stdev_exp_tem = sd(relevant_expression_tem, na.rm=TRUE)
	r_value_tem = cor(relevant_tem, relevant_expression_tem, method="pearson", use="pairwise.complete.obs")
	mcg_min_tem = min(relevant_tem, na.rm=TRUE)
	mcg_max_tem = max(relevant_tem, na.rm=TRUE)
	cat(paste(this_gene_id, mean_mcg_tem, stdev_mcg_tem, mean_exp_tem, stdev_exp_tem, r_value_tem, mcg_min_tem, mcg_max_tem, "\n", sep="\t"))
#, mean_mcg_tem, stdev_mcg_tem, mean_exp_tem, stdev_exp_tem, r_value_tem, mcg_min_tem, mcg_max_tem
#, mean_mcg_gbm, stdev_mcg_gbm, mean_exp_gbm, stdev_exp_gbm, r_value_gbm, mcg_min_gbm, mcg_max_gbm
  } else {
    cat(paste0(this_gene_id," has ", length(best_group), " best haplogroups.\n"))
  }
}





install.packages("Rfast")
library(Rfast)

colMins(big_table_2b)
