### This part optimises the mCG cutoff between unmethylated and methylated CG sites by aiming to maximise agreement between rep siblingc across all sites

max_concordance_fraction = 0
optimal_meth_cutoff = 0
for (partial_meth_cutoff in seq(0.1, 0.9, by = 0.05)) {
	#Alternative approach: consider <50% m calls to be U-
	meth_data[,1] = ifelse(fisher_class=="I","I",ifelse(binom_class=="U","U",ifelse((coverage_data$Cov_C/(coverage_data$Cov_C+coverage_data$Cov_T))>=partial_meth_cutoff,"M","U")))

	# Choose the relevant p-value to retain, dependent on which status was chosen
	#### This needs more work
	#### Fisher's method (using a chisquared test) can be used to combine multiple p-values addressing the same general hypothesis.  Implemented in chiCombP function in the methylPipe package.
	#meth_data[,2] = ifelse((meth_data_fisher[,1]<fisher_p_value_cutoff),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),pmax(meth_data_fisher[,1],meth_data_fisher[,2]),meth_data_fisher[,1]),ifelse((meth_data_fisher[,2]<fisher_p_value_cutoff),meth_data_fisher[,2],pmin(meth_data_fisher[,1],meth_data_fisher[,2])))
	meth_data[,2] = 0
	
	# Discard status and replace with "I" if coverage is too high to consider trusting the call (5 sigma)
	# This assumes that coverage is roughly Gaussian distributed (although it typically is long-tailed)
	for(sample_name in rownames(sample_info)) {
		coverage_upper_cutoff = sd(coverage_data[Sample==sample_name,Cov_C+Cov_T],na.rm=TRUE)*5
		meth_data[(coverage_data$Sample==sample_name) & (!is.na(coverage_data$Cov_C)) & ((coverage_data$Cov_C+coverage_data$Cov_T)>coverage_upper_cutoff),1] = "I"
	}
	
	# Organise meth_data from the large matrix to a data frame
	meth_data=data.frame(meth_data)
	colnames(meth_data) = c("meth_status","p_value")
	#meth_data$p_value=as.numeric(levels(meth_data$p_value)[meth_data$p_value])
	
	meth_data=data.table(meth_data)
	
	# Create empty tables to accumulate each of the 'all samples' columns data
	all_samples_meth_status=NULL
	all_samples_p_value=NULL

	sample_no=0
	for(sample_name in rownames(sample_info)) {
		sample_no = sample_no+1
		#if(opt$verbose) {cat(paste0("Adding methylation calls to all-samples per-site table: ",sample_name,"\n"))}

		sample_meth_status=NULL
		sample_p_value=NULL
		
		# Make separate tables for the sample's meth_call, p_value and meth_call_basis, keyed on locus
		sample_meth_status=cbind(coverage_data,meth_data)[coverage_data$Sample==sample_name,c(2,3,4,7)]
		sample_p_value=cbind(coverage_data,meth_data)[coverage_data$Sample==sample_name,c(2,3,4,8)]
		colnames(sample_meth_status)[4]=sample_name
		colnames(sample_p_value)[4]=sample_name

		# Merge each of the sample's 3 tables into the relevant accumulation table
		if (sample_no==1) {
			all_samples_meth_status=sample_meth_status
			all_samples_p_value=sample_p_value
		} else {
			all_samples_meth_status=merge(all_samples_meth_status, sample_meth_status, by=c("Chromosome","Locus","Strand"), all=TRUE)
			all_samples_p_value=merge(all_samples_p_value, sample_p_value, by=c("Chromosome","Locus","Strand"), all=TRUE)
		}

		# Tidy up the temporary objects
		rm(sample_meth_status, sample_p_value, sample_meth_call_basis)
	} # end for each sample

	#Add a column to identify the number of samples with missing data
	#######CG_pair_meth_status$m_count <- apply(CG_pair_meth_status, 1, function(x) sum(is.na(x))/2)
	
	all_samples_meth_status = data.frame(all_samples_meth_status)
	all_samples_p_value = data.frame(all_samples_p_value)

	# Identify site with either an "M" or "U" call in all samples
	# Identify sites among these without identical calls in all samples

	# Some samples may need to be excluded from this analysis
	# For Becker and Schmitz data sets we want to exclude Line 69 as it is a methylation mutant
	excluded_samples = c()
	valid_samples=c()
	if (project_id=="SRA035939") {
		excluded_samples = c("SRR342380","SRR342390")
	} else if (project_id=="PRJEB2678") {
		excluded_samples = c("ERR046563","ERR046564")
	}
	valid_samples=rownames(sample_info)[!(rownames(sample_info) %in% excluded_samples)]
	valid_samples_meth_status=all_samples_meth_status[,!names(all_samples_meth_status) %in% excluded_samples]

	# Accumulate methylation calls which are concordant among replicates per line_generation

	line_generation_no = 0
	for(this_line in unique(sample_info$Line)) {
		for(this_generation in unique(sample_info[sample_info$Line==this_line,]$Generation)) {
			sample_count=0
			reps_concordant_meth_status=NULL
			for(sample_name in rownames(sample_info[(!rownames(sample_info) %in% excluded_samples) & sample_info$Line==this_line & sample_info$Generation==this_generation,])) {
				#cat(paste(this_line, this_generation, sample_name, "\n"))
				sample_count=sample_count+1
				if (sample_count == 1) {
					# This is the first replicate of replicate set
					line_generation_no = line_generation_no + 1
					# Cache the methylation status and move on
					reps_concordant_meth_status=all_samples_meth_status[,sample_name]
				} else {
					# We are doing a second or higher rep - check for concordance between this sample and the consensus of the previous reps_concordant_meth_status
					# If the current sample disagrees with the previous ones, record a status of "D" for disagreement/discordant for the locus, unless status is opposite (M/U) or "O" already found, in which case "O"
					#reps_concordant_meth_status=as.factor(ifelse(levels(reps_concordant_meth_status)[reps_concordant_meth_status]==all_samples_meth_status[,sample_name],levels(reps_concordant_meth_status)[reps_concordant_meth_status],ifelse(((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="M") & (all_samples_meth_status[,sample_name]=="U")) | ((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="U") & (all_samples_meth_status[,sample_name]=="M")) | ((levels(reps_concordant_meth_status)[reps_concordant_meth_status]=="O") & ((all_samples_meth_status[,sample_name]=="M") | (all_samples_meth_status[,sample_name]=="U"))),"O","D")))
					reps_concordant_meth_status=as.factor(ifelse(reps_concordant_meth_status==all_samples_meth_status[,sample_name],reps_concordant_meth_status,ifelse(((reps_concordant_meth_status=="M") & (all_samples_meth_status[,sample_name]=="U")) | ((reps_concordant_meth_status=="U") & (all_samples_meth_status[,sample_name]=="M")) | ((reps_concordant_meth_status=="O") & ((all_samples_meth_status[,sample_name]=="M") | (all_samples_meth_status[,sample_name]=="U"))),"O","D")))
					### Could do something here with the Fisher method to combine p-values from the concordant reps
				}
			}
			# All reps done for this line and generation, report and cache the results
			if (sample_count>0) {
				if (line_generation_no==1) {
					all_reps_meth_status=cbind(all_samples_meth_status[,c("Chromosome","Locus","Strand")],reps_concordant_meth_status)
					#all_reps_p_value=sample_p_value
				} else {
					all_reps_meth_status=cbind(all_reps_meth_status, reps_concordant_meth_status)
					#all_reps_p_value=cbind(all_samples_p_value, sample_p_value, by=c("Chromosome","Locus","Strand"), all=TRUE)
				}
				colnames(all_reps_meth_status)[line_generation_no+3]=paste0("Line",this_line,"Gen",this_generation)
				#cat(paste0(sum(all_reps_meth_status[,line_generation_no+3]=="O")," sites have opposite calls (M/U) among reps.\n"))
				# Tidy up the temporary objects
				rm(reps_concordant_meth_status, sample_count)
			}
		}
	}

	# Summarise the results of the rep concordance analysis
	line_gen_no = 0
	rep_concordance_summary=matrix(, nrow=ncol(all_reps_meth_status)-3, ncol=6)
	rownames(rep_concordance_summary)=colnames(all_reps_meth_status[,4:ncol(all_reps_meth_status)])
	for(line_gen in rownames(rep_concordance_summary)) {
		line_gen_no = line_gen_no + 1
		#cat(paste0("Summarising methylation call basis for ",line_gen,"\n"))
		status_counts = count(all_reps_meth_status[,line_gen])
		# Needs a hack for when "P" is zero (nrow should be 6 otherwise):
		if(nrow(status_counts)==5) {
			status_counts$x = as.character(status_counts$x)
			status_counts[6,]=status_counts[5,]
			status_counts[5,]=c("P",0)
			status_counts$x = as.factor(status_counts$x)
			status_counts$freq=as.numeric(status_counts$freq)
		}
		rep_concordance_summary[line_gen_no,]=status_counts$freq
		colnames(rep_concordance_summary)=status_counts$x
	}
	# output the summary of numbers of concordant and discordant sites
	#rep_concordance_summary
	# output the proportion of sites in each line with a 'clean' M/U call
	this_concordance_fraction = mean((rep_concordance_summary[,3]+rep_concordance_summary[,6])/sum(rep_concordance_summary[,1:6])*nrow(rep_concordance_summary))
	if (this_concordance_fraction > max_concordance_fraction) {
		max_concordance_fraction = this_concordance_fraction
		optimal_meth_cutoff = partial_meth_cutoff
	}
	cat(paste0("M/U cutoff: ",partial_meth_cutoff," Mean Rep-Concordant call proportion:",this_concordance_fraction,"\n"))
}
cat(paste0("Optimal mCG cutoff between 'U' and 'M' calls = ",optimal_meth_cutoff," giving ",max_concordance_fraction," concordant calls between reps.\n"))

#M/U cutoff: 0.1 Mean Rep-Concordant call proportion:0.801618844342882
#M/U cutoff: 0.15 Mean Rep-Concordant call proportion:0.80170362488033
#M/U cutoff: 0.2 Mean Rep-Concordant call proportion:0.802286382741762
#M/U cutoff: 0.25 Mean Rep-Concordant call proportion:0.802397928030959
#M/U cutoff: 0.3 Mean Rep-Concordant call proportion:0.802279296455108
#M/U cutoff: 0.35 Mean Rep-Concordant call proportion:0.802321559272633
#M/U cutoff: 0.4 Mean Rep-Concordant call proportion:0.802457575192001
#M/U cutoff: 0.45 Mean Rep-Concordant call proportion:0.802505038018438
#M/U cutoff: 0.5 Mean Rep-Concordant call proportion:0.802423979055995
#M/U cutoff: 0.55 Mean Rep-Concordant call proportion:0.801996762637589
#M/U cutoff: 0.6 Mean Rep-Concordant call proportion:0.801106082679836
#M/U cutoff: 0.65 Mean Rep-Concordant call proportion:0.799226534360767
#M/U cutoff: 0.7 Mean Rep-Concordant call proportion:0.79578927748986
#M/U cutoff: 0.75 Mean Rep-Concordant call proportion:0.790470025236867
#M/U cutoff: 0.8 Mean Rep-Concordant call proportion:0.779225668712224
#M/U cutoff: 0.85 Mean Rep-Concordant call proportion:0.757954726377315
#M/U cutoff: 0.9 Mean Rep-Concordant call proportion:0.729341371459527

#Optimal mCG cutoff between 'U' and 'M' calls = 0.45 giving 0.802505038018438 concordant calls between reps.
 
 
 
 
 
 ### This part optimises the number of CG sites per segment by aiming to maximise AUROC, using the annotation as 'truth'
 ### Before it's run again, needs to have potential_CG_segments.gr replaced by MG_segments.gr
 
 # First we make the comparison between Unmethylated + GBM genes and Transposons + Heterochromatic genes
	
	# Create data frame to store the results of the optimisation:
	# For each value of m and n parameters, how many variable sites are overlapping with UMR and LMR segments? What proportion of supposed Transposons and Heterochromatic gene space overlaps with UMR or LMR segments (fp)?  What proportion of supposed Unmenthylated or GBM methylated gene space overlaps with UMR or LMR segments (tp)? 
	m_n_optimisation=data.frame(m=numeric(), n=numeric(), UMR=integer(), LMR=integer(), fp=numeric(), tp=numeric())

	# Loop through values for m and n parameters 
	m.seq=seq(0.005,0.1, by=0.005)
	n.seq=seq(1,15, by=1)
	for (m.sel in m.seq) {
		for (n.sel in n.seq) {
			#n.sel=4
			#UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth_parents_CG.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)
			UMRLMRsegments.gr <- segmentUMRsLMRs(m=potential_CG_segments.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=Athaliana, seqLengths=sLengths)

			# Capitalise chromosome names in segmentation object
			levels(UMRLMRsegments.gr@seqnames)=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")
			UMRLMRsegments.gr@seqinfo@seqnames=c("CHR1","CHR2","CHR3","CHR4","CHR5","CHRM","CHRC")

			# Find overlaps between segments and positive and negative 'control' loci (annotated genes and transposons)
			fp_hits=findOverlaps(methylated_loci.gr, UMRLMRsegments.gr)
			fp_overlaps <- pintersect(methylated_loci.gr[queryHits(fp_hits)], UMRLMRsegments.gr[subjectHits(fp_hits)])
			fp_overlap_prop <- sum(width(fp_overlaps)) / sum(width(methylated_loci.gr))
			tp_hits=findOverlaps(unmethylated_loci.gr, UMRLMRsegments.gr)
			tp_overlaps <- pintersect(unmethylated_loci.gr[queryHits(tp_hits)], UMRLMRsegments.gr[subjectHits(tp_hits)])
			tp_overlap_prop <- sum(width(tp_overlaps)) / sum(width(unmethylated_loci.gr))
			
			# Find overlaps between segments and variable loci
			variable_sites_UMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="UMR"])
			variable_sites_LMR_olaps = findOverlaps(variant_call_ranges, UMRLMRsegments.gr[UMRLMRsegments.gr@elementMetadata$type=="LMR"])
			m_n_optimisation=rbind(m_n_optimisation, list(m=m.sel, n=n.sel, UMR=length(variable_sites_UMR_olaps@from), LMR=length(variable_sites_LMR_olaps@from), fp=fp_overlap_prop, tp=tp_overlap_prop))
			# How many variable sites were captured?
			cat(paste0("m.sel=",m.sel," n.sel=",n.sel," Variable sites in UMRs: ",length(variable_sites_UMR_olaps@from),"  LMRs: ",length(variable_sites_LMR_olaps@from),"\n"))
		}
	}
	m_n_optimisation$variable_sites_captured=m_n_optimisation$UMR+m_n_optimisation$LMR

	# Plot ROC curves for each value of n, and estimate AUC
	pdf(paste0(project_id,"_",meth_context,"_MethylSeekR_m_n_optimisation_parents_ROC.pdf"))
	print(ggplot(m_n_optimisation, aes(x=fp, y=tp)) + geom_point() +geom_text(aes(label=paste0(m)),hjust=0, vjust=0) +geom_line() +facet_wrap(~ n))
	dev.off()
	
	# Estimate AUC from trapezium approximation
	cat(paste0("Estimated AUC values from ROC curves generated by varying m parameter for each value of n in MethylSeekR:\n"))
	best_n = 0
	prev_best_auc = 0
	for (n.sel in n.seq) {
		m_n_optimisation_auc=0
		prev_tp=0
		prev_fp=0
		for (m.sel in m.seq) {
			m_n_optimisation_auc = m_n_optimisation_auc + (1-m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$fp)*(m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$tp-prev_tp) + ((m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$tp-prev_tp)*(m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$fp-prev_fp))/2
			prev_tp=m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$tp
			prev_fp=m_n_optimisation[(m_n_optimisation$m==m.sel) & (m_n_optimisation$n==n.sel),]$fp
		}
		if (m_n_optimisation_auc > prev_best_auc) {
			best_n = n.sel
			prev_best_auc = m_n_optimisation_auc
		}
		cat(paste0(" n=",n.sel," AUC=",m_n_optimisation_auc,"\n"))
	}
	cat(paste0("n=",best_n," maximises AUC (",prev_best_auc,").\n"))
	
	# Schmitz data:

	#n=1 AUC=0.927485882894315
	#n=2 AUC=0.928660724765878
	#n=3 AUC=0.93052386921611
	#n=4 AUC=0.933146271390467
	#n=5 AUC=0.934533888056953
	#n=6 AUC=0.935387796748723
	#n=7 AUC=0.935456904324158
	#n=8 AUC=0.93575810599037
	#n=9 AUC=0.935844509791035 ***
	#n=10 AUC=0.935423601423369
	#n=11 AUC=0.935206312103168
	#n=12 AUC=0.934498784067983
	#n=13 AUC=0.934116963065039
	#n=14 AUC=0.933816992571899
	#n=15 AUC=0.933719054896366
	#n=16 AUC=0.933275398107931
	#n=17 AUC=0.932807108263812
	#n=18 AUC=0.932438172207426
	#n=19 AUC=0.931341560763276
	#n=20 AUC=0.930869799202132
	
	# n=9 generates the largest AUC.  m=0.12 captures 3994 variable sites in UMR segments.  Previously we identified 2.5% as an appropriate threshold for unmethylated segments, so we will try this.
	
	#n=1 AUC=0.928405597073089
	# n=2 AUC=0.929596914114101
	# n=3 AUC=0.931435630612216
	# n=4 AUC=0.933776693439385
	# n=5 AUC=0.935008031513872
	# n=6 AUC=0.935827266918043
	# n=7 AUC=0.935908380315107
	# n=8 AUC=0.936235681827533
	# n=9 AUC=0.936370415538498
	# n=10 AUC=0.936088277156938
	# n=11 AUC=0.935879866004354
	# n=12 AUC=0.935241736511006
	# n=13 AUC=0.934889800859761
	# n=14 AUC=0.934620701902889
	# n=15 AUC=0.934469145258861
	# n=16 AUC=0.934005460774851
	# n=17 AUC=0.933478637857632
	# n=18 AUC=0.933066764714414
	# n=19 AUC=0.931984252023881
	# n=20 AUC=0.931537191673741

	#n=9 maximises AUC (0.936370415538498).
	