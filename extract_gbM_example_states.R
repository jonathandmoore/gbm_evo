
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    model_site_olaps = findOverlaps(this_model_trimmed, CG_site_ranges)

	
	
# Now identify accessions with Col-0-like gbM levels
relevant_gbM_loci = cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",]
colnames(relevant_gbM_loci)[1]="gene_ID"
relevant_gbM_loci = relevant_gbM_loci[!is.na(relevant_gbM_loci$gene_firsts),]
relevant_gbM_loci.gr = makeGRangesFromDataFrame(relevant_gbM_loci, seqnames.field = "seqnames", start.field = "gene_firsts", end.field = "gene_lasts")
relevant_CG_sites = queryHits(findOverlaps(CG_site_ranges, relevant_gbM_loci.gr))

accession_gbm_levels = rep(0, length(list.files(path=model_dir)))
this_model_no = 0
for (this_model_file in list.files(path=model_dir)) {
  this_model_no = this_model_no + 1
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  names(accession_gbm_levels)[this_model_no] = this_sample
  this_model = read.table(paste0("segmentation_models/m1001_CG_", this_sample,"_segmentation_model_draft3.tsv"), sep="\t", header=TRUE)
  these_tem_loci = makeGRangesFromDataFrame(df=this_model[substr(this_model$Type,1,3)=="TEM",], seqnames.field="Chromosome", start.field="Start", end.field="End")
  # Exclude loci which are teM in the accession from the list of relevant CG sites to assess accession global gbM
  these_relevant_CG_sites = queryHits(findOverlaps(CG_site_ranges, setdiff(relevant_gbM_loci.gr, these_tem_loci)))
 
  accession_gbm_levels[this_model_no] = sum(all_samples_meth_status[these_relevant_CG_sites,this_sample]=="M", na.rm=TRUE) / sum(all_samples_meth_status[these_relevant_CG_sites,this_sample] %in% c("M","U"), na.rm=TRUE)
}

hist(accession_gbm_levels, breaks = seq(0,1, 0.01))

mean(accession_gbm_levels)
#[1] 0.2993232
sd(accession_gbm_levels)
#[1] 0.04787389

mean(accession_gbm_levels)-sd(accession_gbm_levels)
#[1] 0.2514493
mean(accession_gbm_levels)+sd(accession_gbm_levels)
#[1] 0.3471971

# Col-0 accession levels:
accession_gbm_levels[c("SRX446038", "SRX446039", "SRX446040", "SRX248644")]
#SRX446038 SRX446039 SRX446040 SRX248644 
#0.2308716 0.2850365 0.3006594 0.3499154

# all 4 Col-0 accessions have weirdly differing levels of mCG in gbM. Need to load in Schmitz data as a check

schmitz_all_samples_meth_status = readRDS("../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_all_samples_meth_status.rds")
schmitz_parental_consensus = read.table(file="../../Jay-SMP_variation/5-analysis/SRA035939/parental_consensus.txt")
schmitz_samples_meth_status = cbind(schmitz_all_samples_meth_status, schmitz_parental_consensus)[,c("Chromosome","Locus","V1")]
colnames(schmitz_samples_meth_status)[3]="Schmitz_parental"
schmitz_samples_meth_status$Chromosome = paste0("Chr",substr(schmitz_samples_meth_status$Chromosome,4,4))
merged_samples_meth_status = merge(all_samples_meth_status[,1:2], schmitz_samples_meth_status, by=c("Chromosome","Locus"), all.x=TRUE)
schmitz_samples_meth_status = schmitz_samples_meth_status[schmitz_samples_meth_status$Chromosome %in% c("Chr1","Chr2","Chr3","Chr4","Chr5"),]
schmitz_samples_meth_status$Chromosome = substr(schmitz_samples_meth_status$Chromosome,4,4)

segmentation_model.gr = readRDS(file = "../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_segmentation_model_draft3.rds")
these_tem_loci = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="TEM"]

levels(these_tem_loci@seqnames@values) = substr(as.character(levels(these_tem_loci@seqnames@values)),4,4)
these_tem_loci@seqinfo@seqnames = levels(these_tem_loci@seqnames@values)

schmitz_relevant_CG_sites = queryHits(findOverlaps(makeGRangesFromDataFrame(df=schmitz_samples_meth_status, seqnames.field = "Chromosome", start.field = "Locus", end.field="Locus"), setdiff(relevant_gbM_loci.gr, these_tem_loci)))
sum(schmitz_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"]=="M", na.rm=TRUE) / sum(schmitz_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"] %in% c("M","U"), na.rm=TRUE)
#[1] 0.2707315




#*** from here we don't run


# check out the values for the Schmitz GBM segments
#segmentation_model.gr = readRDS(file = "../../Jay-SMP_variation/5-analysis/SRA035939/SRA035939_CG_segmentation_model_draft3.rds")
GBM_segments = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="GBM"]
these_tem_loci = segmentation_model.gr[segmentation_model.gr$segment.mStatus=="TEM"]
relevant_CG_sites = queryHits(findOverlaps(CG_site_ranges, GBM_segments))
these_relevant_CG_sites = queryHits(findOverlaps(CG_site_ranges, setdiff(relevant_gbM_loci.gr, these_tem_loci)))
schmitz_relevant_CG_sites = queryHits(findOverlaps(makeGRangesFromDataFrame(df=schmitz_samples_meth_status, seqnames.field = "Chromosome", start.field = "Locus", end.field="Locus"), GBM_segments))

# Export longer gbM Schmitz segments
schmitz_relevant_gbm_segments = table(subjectHits(findOverlaps(CG_site_ranges, GBM_segments)))
schmitz_relevant_gbm_segments = schmitz_relevant_gbm_segments[schmitz_relevant_gbm_segments>3]

# Amy wants a gene ID for each segment, so do the overlapping:
relevant_genes = gff.genes[findOverlaps(GBM_segments, gene_ranges, select="first"),"gene_ID"]
# Which gene is best for each segment? We take the first one, because we expect that this is at least overlapping at 5', so in most cases, will hopefully be the longest segment

#write.table(as.data.frame(GBM_segments)[names(schmitz_relevant_gbm_segments)[schmitz_relevant_gbm_segments>3],c("seqnames","start","end")], file="Schmitz_gbM_segments_gt_3_CG_sites.tsv", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
write.table(cbind("gene_ID"=relevant_genes, as.data.frame(GBM_segments))[names(schmitz_relevant_gbm_segments)[schmitz_relevant_gbm_segments>3],c("gene_ID","seqnames","start","end")], file="Schmitz_gbM_segments_gt_3_CG_sites.tsv", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
# NB This produces duplicate entries for overlapping genes. As a fix, the first of each of such entries was chosen, in Excel, before the file was passed to Amy for modelling


accession_gbm_levels = rep(0, length(list.files(path=model_dir)))
this_model_no = 0
for (this_model_file in list.files(path=model_dir)) {
  this_model_no = this_model_no + 1
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  names(accession_gbm_levels)[this_model_no] = this_sample
  accession_gbm_levels[this_model_no] = sum(all_samples_meth_status[relevant_CG_sites,this_sample]=="M", na.rm=TRUE) / sum(all_samples_meth_status[relevant_CG_sites,this_sample] %in% c("M","U"), na.rm=TRUE)
}

hist(accession_gbm_levels, breaks = seq(0,1, 0.01))

mean(accession_gbm_levels)
#[1] 0.5772374
sd(accession_gbm_levels)
#[1] 0.06789811

mean(accession_gbm_levels)-sd(accession_gbm_levels)
#[1] 0.5093393
mean(accession_gbm_levels)+sd(accession_gbm_levels)
#[1] 0.6451355

# Col-0 accession levels:
accession_gbm_levels[c("SRX446038", "SRX446039", "SRX446040", "SRX248644")]
#SRX446038 SRX446039 SRX446040 SRX248644 
#0.5794306 0.6366895 0.6524200 0.7174731

sum(schmitz_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"]=="M", na.rm=TRUE) / sum(merged_samples_meth_status[schmitz_relevant_CG_sites,"Schmitz_parental"] %in% c("M","U"), na.rm=TRUE)
#[1] 0.6124375


#*** we pick up again here



# Let's export a spreadsheet of metadata combined with accession_gbm_levels
# read sample info in
raw_data_dir = "../0-raw_data/"
sample_list = read.table(paste0(raw_data_dir,"sample_info_synthesis.txt"), sep="\t", header=TRUE)
sample_list$Name = str_split_fixed(sample_list$Title, "\\(",2)[,1]
#sample_list1$Ecotype_id = as.numeric(str_split_fixed(str_split_fixed(sample_list1$Title, "\\(",2)[,2], "\\)", 2)[,1])
sample_list$Ecotype_id = as.numeric(str_split_fixed(sample_list$Sample, "_", 2)[,1])

accession_gbm_levels = as.data.frame(accession_gbm_levels)
accession_gbm_levels$SRA.Accession = row.names(accession_gbm_levels)
sample_list = merge(sample_list, accession_gbm_levels, by="SRA.Accession", all=TRUE)

# before we do this we need to make some space
# as well as removing these, I closed Chrome and Teams
rm(schmitz_all_samples_meth_status)
rm(schmitz_samples_meth_status)
gc()

# add in proportion of sites well-covered as a proxy for coverage
#accession_coverage_ratio = colSums((all_samples_meth_status[,4:ncol(all_samples_meth_status)]=="M") | (all_samples_meth_status[,4:ncol(all_samples_meth_status)]=="U"), na.rm=TRUE)
#accession_coverage_ratio = colSums((all_samples_meth_status[,4:1214]=="M") | (all_samples_meth_status[,4:1214]=="U"), na.rm=TRUE)
#accession_coverage_ratio = colSums((all_samples_meth_status[,4:1214] %in% c("M","U")), na.rm=TRUE)
#accession_coverage_ratio = colSums((all_samples_meth_status[,4:ncol(all_samples_meth_status)] %in% c("M","U")), na.rm=TRUE)
accession_coverage_ratio = rep(NA, ncol(all_samples_meth_status))
names(accession_coverage_ratio) = colnames(all_samples_meth_status)
for (this_col in 4:ncol(all_samples_meth_status)) {
  accession_coverage_ratio[this_col] = sum((all_samples_meth_status[,this_col]=="M") | (all_samples_meth_status[,this_col]=="U"), na.rm=TRUE)
}
accession_coverage_ratio = accession_coverage_ratio[4:ncol(all_samples_meth_status)]
accession_coverage_ratio = accession_coverage_ratio / nrow(all_samples_meth_status)
accession_coverage_ratio = as.data.frame(accession_coverage_ratio)
accession_coverage_ratio$SRA.Accession = row.names(accession_coverage_ratio)
sample_list = merge(sample_list, accession_coverage_ratio, by="SRA.Accession", all=TRUE)

write.table(sample_list, file="m1001_samples_mCG_density_coverage_2.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#How many accessions have methylation within 1SD of the mean?
#accession_gbm_levels is now a data frame, let's turn it back into a vector
accession_gbm_level_names = accession_gbm_levels$SRA.Accession
accession_gbm_levels=accession_gbm_levels$accession_gbm_levels
names(accession_gbm_levels)=accession_gbm_level_names
sum(accession_gbm_levels>=mean(accession_gbm_levels)-sd(accession_gbm_levels) & accession_gbm_levels<=mean(accession_gbm_levels)+sd(accession_gbm_levels))  
#sum(accession_gbm_levels$accession_gbm_levels>=mean(accession_gbm_levels$accession_gbm_levels)-sd(accession_gbm_levels$accession_gbm_levels) & accession_gbm_levels$accession_gbm_levels<=mean(accession_gbm_levels$accession_gbm_levels)+sd(accession_gbm_levels$accession_gbm_levels))  
#[1] 884
  
#These 884 accessions are the 'Col-0-like' ones, or at least, Col-0 and these all lie within 1SD of the mean of the methylation level of the Col-0 gbM segments
#So we are interested in genes which contain gbM segments in at least 20% of this group

relevant_samples = names(accession_gbm_levels[accession_gbm_levels>=mean(accession_gbm_levels)-sd(accession_gbm_levels) & accession_gbm_levels<=mean(accession_gbm_levels)+sd(accession_gbm_levels)])

gene_gbm_counts = rep(0, length(gene_ranges))
no_models = 0
for (this_model_file in list.files(path=model_dir)) {
  this_sample = substr(this_model_file, 1, str_length(this_model_file)-16)
  if (this_sample %in% relevant_samples) {
    no_models = no_models + 1
    cat(paste0(this_sample,"\n"))
    this_model = read.table(paste0(model_dir, "/", this_model_file), sep="\t", header=TRUE)
    this_model$Chromosome = paste0("Chr",this_model$Chromosome)
    this_model.gr = makeGRangesFromDataFrame(df=this_model, seqnames.field="Chromosome", start.field="Start", end.field="End")
    this_model_trimmed = this_model.gr[as.vector(table(queryHits(findOverlaps(this_model.gr, CG_site_ranges)))>=min_segment_sites)]
    model_gene_olaps = findOverlaps(this_model_trimmed, gene_ranges)
    gene_gbm_counts[subjectHits(model_gene_olaps)] = gene_gbm_counts[subjectHits(model_gene_olaps)] + 1
    
    #if (no_models == 1) {
    #  overlap_model = this_model_trimmed
    #} else {
    #  overlap_model = union(overlap_model, this_model_trimmed, ignore.strand=TRUE)
    #}
  }
}

relevant_gbM_loci = (cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts, gene_gbm_counts)[substr(gff.genes$gene_name,1,2)=="AT",])
write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts, gene_gbm_counts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,"_gbm_counts.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# limit it to genes with gbM segments in at least 20% of Col-0-like accessions
relevant_gbM_loci = relevant_gbM_loci[(!is.na(relevant_gbM_loci$gene_firsts)) & (relevant_gbM_loci$gene_gbm_counts>884/5),]
relevant_gbM_loci.gr = makeGRangesFromDataFrame(relevant_gbM_loci, seqnames.field = "seqnames", start.field = "gene_firsts", end.field = "gene_lasts")

levels(relevant_gbM_loci.gr@seqnames@values) = substr(as.character(levels(relevant_gbM_loci.gr@seqnames@values)),4,4)
relevant_gbM_loci.gr@seqinfo@seqnames = levels(relevant_gbM_loci.gr@seqnames@values)


#*** we don't bother with this


# The above method calculates relevant gbM_loci by taking first site to last site which are covered by gbM segments of at least 3 sites in at least 5 accessions.
# We now want to only include segments which meet that criteria, rather than potentially unmethylated segments between two methylated ones

# We traverse CG_site_gbM_counts one gene at a time (which already contains counts for each site of overlapping gbM segments longer than 3 sites) and identify contiguous segments within each gene which are covered by relevant gbM segments in at least 5 accessions

#for (min_seg_cover in 1:10) {
for (min_seg_cover in 5) {
  #gene_firsts = rep(0, length(gene_ranges))
  #gene_lasts = rep(0, length(gene_ranges))
  gene_gbM_segments = data.frame()
  for (this_gene in 1:length(gene_ranges)) {
    cat(paste0("cover: ", min_seg_cover, "  gene: ", this_gene, "\n"))
    doing_one = FALSE
    this_chr = as.data.frame(gene_ranges)[this_gene,"seqnames"]
    this_start = 0
    this_end = 0
    gene_site_olaps = findOverlaps(gene_ranges[this_gene], CG_site_ranges)
    gene_site_counts = CG_site_gbM_counts[subjectHits(gene_site_olaps)]
    no_sites = length(gene_site_olaps)
    if (no_sites>0) {
      prev_site = 0
      last_site = FALSE
      save_this = FALSE
      for (this_site in 1:no_sites) {
        if (this_site==no_sites) {}
        if (doing_one) {
          if (gene_site_counts[this_site]>=min_seg_cover) {
            if (this_site==no_sites) {
              this_end = this_site
              save_this = TRUE
            }
          } else {
            this_end = this_site -1
            doing_one = FALSE
            save_this = TRUE
          }
        } else {
          if (gene_site_counts[this_site]>=min_seg_cover) {
            doing_one=TRUE
            this_start = this_site
            if (this_site==no_sites) {
              this_end = this_site
              save_this = TRUE
            }
          }
        }
        if (save_this) {
          gene_gbM_segments = rbind(gene_gbM_segments, c("gene_ID"=gff.genes[this_gene,"gene_ID"], "chromosome"=this_chr, "gene_start" = gff.genes[this_gene,"V4"], "gene_end" = gff.genes[this_gene,"V5"], "start"=CG_site_ranges[subjectHits(gene_site_olaps)[this_start]]@ranges@start, "end"=CG_site_ranges[subjectHits(gene_site_olaps)[this_end]]@ranges@start+1))
          colnames(gene_gbM_segments) = c("gene_ID", "chromosome", "gene_start", "gene_end", "start", "end")
          save_this=FALSE
        }
      }
    } else {
      #gene_firsts[this_gene] = NA
      #gene_lasts[this_gene] = NA
    }
  }
  #write.table(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",], file=paste0("gene_gbm_extents_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(gene_gbM_segments, file=paste0("gene_gbm_segments_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  # These are all segments which meet the coverage criteria. If we want to limit it to genes which contain gbM in 20% of Col-0-like accessions:
  gene_gbM_segments.gr = makeGRangesFromDataFrame(gene_gbM_segments, seqnames.field="chromosome", start.field="start", end.field="end")
  write.table(gene_gbM_segments[queryHits(findOverlaps(gene_gbM_segments.gr, relevant_gbM_loci.gr)),], file=paste0("gene_gbm_segments_min_accesions_",min_seg_cover,"_min_CG_",min_segment_sites,"_20percent.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

View(cbind(gff.genes[,c("gene_ID")],as.data.frame(gene_ranges), gene_firsts, gene_lasts)[substr(gff.genes$gene_name,1,2)=="AT",])


#*** we pick up again here


# partition the sample list to find some for modelling
# no need to bother with anything below accession_coverage_ratio of 0.6 - that still allows us to include the most ad least methylated ones (apart from the three lowly methylated accessions with coverage<10%)
relevant_samples = sample_list[sample_list$accession_coverage_ratio>0.6,]

# Pick the 10 most and least methylated accessions
# 10 least methylated accessions
least_methylated = head(relevant_samples[order(relevant_samples$accession_gbm_levels),], 10)
# 10 most methylated
most_methylated = head(relevant_samples[order(-relevant_samples$accession_gbm_levels),], 10)

# remainder in the middle:
middle_methylated = relevant_samples[(relevant_samples$accession_gbm_levels>max(least_methylated$accession_gbm_levels)) & (relevant_samples$accession_gbm_levels<min(most_methylated$accession_gbm_levels)),]

min(least_methylated$accession_gbm_levels)
#[1] 0.1180222
max(least_methylated$accession_gbm_levels)
#[1] 0.2006623
min(most_methylated$accession_gbm_levels)
#[1] 0.4102528
max(most_methylated$accession_gbm_levels)
#[1] 0.5467242

min(middle_methylated$accession_gbm_levels)
#[1] 0.2027648
max(middle_methylated$accession_gbm_levels)
#[1] 0.4068121

#deciles for the in-between accessions:
modelling_deciles = quantile(middle_methylated$accession_gbm_levels, prob = seq(0, 1, length = 11), type = 5)

#       0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
#0.2027648 0.2522585 0.2658764 0.2769889 0.2855385 0.2953658 0.3039275 0.3134926 0.3269840 0.3479528 0.4068121 

#1 0.2027648 - 0.2522585
#2 0.2522585 - 0.2658764 
#3 0.2658764 - 0.2769889
#4 0.2769889 - 0.2855385
#5 0.2855385 - 0.2953658 
#6 0.2953658 - 0.3039275 
#7 0.3039275 - 0.3134926 
#8 0.3134926 - 0.3269840 
#9 0.3269840 - 0.3479528
#10 0.3479528 - 0.4068121 


# pick the accessions for export to Amy for modelling
# we are going to choose 10 sets of 10 accessions each, each set of 10 to represent a segment from the distribution of mCG density in Col-0 gbM segments

site_segment_olaps = findOverlaps(CG_site_ranges, relevant_gbM_loci.gr)
relevant_CG_sites = queryHits(site_segment_olaps)
relevant_gene_IDs = relevant_gbM_loci[subjectHits(site_segment_olaps),1]

middle_start = min(middle_methylated$accession_gbm_levels)
middle_end = max(middle_methylated$accession_gbm_levels)

#write.table(cbind("gene_ID"=relevant_gene_IDs, all_samples_meth_status[relevant_CG_sites, c("Chromosome", "Locus", least_methylated$SRA.Accession)]), file=paste0("m1001_samples_meth_status_mCG_density_low.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(cbind("gene_ID"=relevant_gene_IDs, all_samples_meth_status[relevant_CG_sites, c("Chromosome", "Locus", most_methylated$SRA.Accession)]), file=paste0("m1001_samples_meth_status_mCG_density_high.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

modelling_accessions=list()
for (this_section in 1:10) {
  #section_start = middle_start + (this_section-1)*(middle_end-middle_start)/10
  #section_end = middle_start + (this_section)*(middle_end-middle_start)/10
  section_start = modelling_deciles[this_section]
  section_end = modelling_deciles[this_section+1]
  modelling_accessions[[this_section]] = head(middle_methylated[order(-middle_methylated$accession_coverage_ratio),][(middle_methylated$accession_gbm_levels>=section_start) & (middle_methylated$accession_gbm_levels<section_end),"SRA.Accession"],10)
  cat(paste0(modelling_accessions[[this_section]]), "\n")  
}

modelling_accessions[[11]] = least_methylated$SRA.Accession
modelling_accessions[[12]] = most_methylated$SRA.Accession

site_gene_olaps = findOverlaps(CG_site_ranges, gene_ranges)
whole_gene_sites = queryHits(site_gene_olaps)
whole_gene_IDs = gff.genes$gene_ID[subjectHits(site_gene_olaps)]

for (this_section in 1:12) {
  #Need to make this gene, chromosome, start, end, status
  for (this_accession in modelling_accessions[[this_section]]) {
    #export segments only
    #write.table(cbind("gene_ID"=relevant_gene_IDs, "Chromosome"=all_samples_meth_status[relevant_CG_sites,"Chromosome"], "Start"=all_samples_meth_status[relevant_CG_sites,"Locus"], "End"=all_samples_meth_status[relevant_CG_sites,"Locus"]+1, "Status"=all_samples_meth_status[relevant_CG_sites,this_accession]), file=paste0("m1001_samples_meth_status_mCG_density_",this_section,"_",this_accession,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    # export whole genes
    write.table(cbind("gene_ID"=whole_gene_IDs, "Chromosome"=all_samples_meth_status[whole_gene_sites,"Chromosome"], "Start"=all_samples_meth_status[whole_gene_sites,"Locus"], "End"=all_samples_meth_status[whole_gene_sites,"Locus"]+1, "Status"=all_samples_meth_status[whole_gene_sites,this_accession]), file=paste0("m1001_samples_meth_status_mCG_density_",this_section,"_",this_accession,"_3.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}

# Initial extremities for modelling:

#Least methylated (both Ecker):
#SRX248646 Cvi_0
#SRX2190740 UKID116 (5822)-2

#Most methylated (Dubin) and very methylated (Ecker):
#SRX445897 5856_10C
#SRX1664724 Gron 12 (9386)


#Now we have exported accessions for modelling, we check the range of values found in Col-0 generally

wt_samples_meth_status = readRDS(file="../../Jay-SMP_variation/5-analysis/AncestralAnalysis_CG_all_samples_meth_status.rds")
wt_samples_meth_status = wt_samples_meth_status[wt_samples_meth_status$Chromosome %in% c("CHR1","CHR2","CHR3","CHR4","CHR5"),]
wt_samples_meth_status$Chromosome = as.integer(substr(wt_samples_meth_status$Chromosome,4,4))
wt_samples_meth_status = merge(wt_samples_meth_status, schmitz_samples_meth_status[schmitz_relevant_CG_sites,], by=c("Chromosome", "Locus"))

wt_samples_global_mCG = rep(0,ncol(wt_samples_meth_status))
for (this_sample in 4:ncol(wt_samples_meth_status)) {
  this_sample_name = colnames(wt_samples_meth_status)[this_sample]
  cat(paste(this_sample, this_sample_name, "\n"))
  wt_samples_global_mCG[this_sample] = sum(wt_samples_meth_status[,this_sample]=="M", na.rm=TRUE) / sum(wt_samples_meth_status[,this_sample] %in% c("M","U"), na.rm=TRUE)
}
names(wt_samples_global_mCG) = colnames(wt_samples_meth_status)
wt_samples_global_mCG = wt_samples_global_mCG[4:length(wt_samples_global_mCG)]
wt_samples_global_mCG = as.data.frame(wt_samples_global_mCG)
wt_samples_global_mCG$sample = row.names(wt_samples_global_mCG)

ggplot(as.data.frame(cbind("global_mCG"=wt_samples_global_mCG, "z"=rep(0, length(wt_samples_global_mCG))))) + geom_histogram(aes(x=global_mCG), binwidth=0.005)

mean(wt_samples_global_mCG)
# [1] 0.2671663
sd(wt_samples_global_mCG)
# [1] 0.03776547

# add in proportion of sites well-covered as a proxy for coverage
wt_accession_coverage_ratio = rep(NA, ncol(wt_samples_meth_status))
names(wt_accession_coverage_ratio) = colnames(wt_samples_meth_status)
for (this_col in 4:ncol(wt_samples_meth_status)) {
  wt_accession_coverage_ratio[this_col] = sum((wt_samples_meth_status[,this_col]=="M") | (wt_samples_meth_status[,this_col]=="U"), na.rm=TRUE)
}
wt_accession_coverage_ratio = wt_accession_coverage_ratio[4:ncol(wt_samples_meth_status)]
wt_accession_coverage_ratio = wt_accession_coverage_ratio / nrow(wt_samples_meth_status)
wt_accession_coverage_ratio = as.data.frame(wt_accession_coverage_ratio)
wt_accession_coverage_ratio$sample = row.names(wt_accession_coverage_ratio)

wt_samples_global_mCG = merge(wt_accession_coverage_ratio, wt_samples_global_mCG, by="sample")

ggplot(wt_samples_global_mCG[wt_samples_global_mCG$wt_accession_coverage_ratio>0.6,]) + geom_histogram(aes(x=wt_samples_global_mCG), binwidth=0.005) +xlab("global gbM")


mean(wt_samples_global_mCG[wt_samples_global_mCG$wt_accession_coverage_ratio>0.6,"wt_samples_global_mCG"])
# [1] 0.2727985
sd(wt_samples_global_mCG[wt_samples_global_mCG$wt_accession_coverage_ratio>0.6,"wt_samples_global_mCG"])
# [1] 0.01789724
#write.table(sample_list, file="m1001_samples_mCG_density_coverage_2.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
