#!/usr/bin/env Rscript

# 0-ETL_1001_methylomes.R
#
# Jay Moore, John Innes Centre, Norwich
# jonathan.moore@jic.ac.uk
#
# 28/10/2019
#
# Version 1 - Staging script for 1001 methylomes data
#
# Change log
#

# Summary of functions
#
# Done:
# 
# Generate batch scripts to grab 1001 methylomes data and prepare it for analysis
#
# Batch up the samples into small groups (10 samples per batch)
# Divide the process up into steps.  
# Step 0 gets run on the software7 node so is isolated from other steps.
# For the remaining steps, sbatch step 1 on the cluster, and when a step completes for its batch of 10 samples it calls the next step in the chain
#
# 0_wget_batch: retrieve the fastq.gz files by ftp from EBI SRA
# 1_dtool_create: make a local dtool dataset for each sample
# 2_dtool_readme: prepare a readme for the dtool dataset
# 3_dtool_freeze: freeze the dtool dataset by making hash
# 4_dtool_copy: transfer the dataset to iRODS
# 5_dtool_tidy: delete the local dtool dataset
# 6_align: bring the dataset back from iRODS to the cluster node SSD and align each run to the reference. Transfer the resulting .cov file to storage
# 7_call_methylation: Call R script to load the .cov file for each run, consolidate runs for each sample, apply models to call methylation status, and save the reulting all_samples_meth_status for the sample
# 
#
# Underway:
#
#
# To do:
# 
#
#
# Results:
#
# 


# R script to assemble raw data from manifest files from SRA and ENA, and generate a script to download and align these data
path = "X:/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/"
setwd(path)
linux_path = "/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/"

# Global Arabidopsis accessions from Ecker study (927 of)
sample_list1 = read.table(paste0(path,"GSE43857.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=43857
fastq_list1 = read.table(paste0(path,"ENA_SRA065807.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/SRA065807

# Swedish Arabidopsis accessions from separate study (284 of)
sample_list2 = read.table(paste0(path,"GSE54292.txt"), sep="\t", header=TRUE)   # from https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=54292
fastq_list2 = read.table(paste0(path,"ENA_SRP035593.txt"), sep="\t", header=TRUE)   # from https://www.ebi.ac.uk/ena/data/view/PRJNA236110

# consolidate the two sets of data
sample_list = rbind.data.frame(sample_list1, sample_list2)
fastq_list = rbind.data.frame(fastq_list1, fastq_list2)

# clean up
rm(sample_list1, fastq_list1, sample_list2, fastq_list2) 

# associate fastq files with runs, and runs with samples
sample_files = data.frame()
for (this_accession in sample_list$SRA.Accession) {
	these_files = fastq_list[fastq_list$experiment_accession == this_accession,"fastq_ftp"]
	for (this_file in these_files) {
		some_files=strsplit(this_file,";")
		sample_files = rbind.data.frame(sample_files, c("Accession" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Accession"]),  "Title" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Title"]), "SRA_Accession" = this_accession, "FASTQ_file" = some_files[[1]][[1]]), stringsAsFactors=FALSE)
		if (length(some_files[[1]]) > 1) {
			sample_files = rbind.data.frame(sample_files, c("Accession" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Accession"]), "Title" = as.character(sample_list[sample_list$SRA.Accession==this_accession,"Title"]), "SRA_Accession" = this_accession, "FASTQ_file" = some_files[[1]][[2]]), stringsAsFactors=FALSE)
			# it looks like a lot of the paired-end files are represented as 'SINGLE' but should be 'PAIRED' - replace where necessary
			if (fastq_list[fastq_list$fastq_ftp == this_file,"library_layout"] == "SINGLE") {
				fastq_list[fastq_list$fastq_ftp == this_file,"library_layout"] = "PAIRED"
			}
		} else {
			#cat(fastq_list[fastq_list$experiment_accession == this_accession,"library_layout"])
		}
	}
}
colnames(sample_files) = c("Accession", "Title", "SRA_Accession", "FASTQ_file")


# break data set up into batches of X samples, and write a separate script to download each batch
batch_size = 10

# Capture the end of each batch run in an overall staging log
log_file = paste0(linux_path, "staging_log.txt")

# Ouput all the steps to a launch script for completeness
launch_script = paste0("launch_script.sh")
launch_script_contents = list()

# Generate the script files for the first step (0)
current_step = "0-wget_batch_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  for (this_file in sample_files[sample_files$SRA_Accession == sample_list[this_accession,"SRA.Accession"],"FASTQ_file"]) {
    script_contents = c(script_contents, paste0("wget ",this_file, "\n"))
  }
  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
    # at the start do stuff (not needed here)
    script_header = ""
    script_header = paste0("if [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\n")
	# at the end a flag file indicating the step is complete
    #script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt")
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (1)
previous_step =  "0-wget_batch_"
current_step = "1-dtool_create_"
next_step = "2-dtool_readme_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # make a new dtool_dataset for the accession
  script_contents = c(script_contents, paste0("dtool create ",sample_list[this_accession,"SRA.Accession"], "\n"))
  for (this_file in sample_files[sample_files$SRA_Accession == sample_list[this_accession,"SRA.Accession"],"FASTQ_file"]) {
    # add each of the files for the accession to the dtool_dataset
	this_filename = strsplit(this_file,"/")[[1]][[length(strsplit(this_file,"/")[[1]])]]
    script_contents = c(script_contents, paste0("mv ", linux_path, this_filename, " ", linux_path, "dtool_datasets/", sample_list[this_accession, "SRA.Accession"], "/data \n"))	
  }
  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (2)
previous_step =  "1-dtool_create_"
current_step = "2-dtool_readme_"
next_step = "3-dtool_freeze_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # make a new dtool readme for the accession based on the pre-prepared README.yml
  script_contents = c(script_contents, paste0("echo --- > ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("echo description: Bisulphite sequencing, SRA accession ", sample_list[this_accession,"SRA.Accession"], " >> ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("cat README.yml >> ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("dtool readme write ",sample_list[this_accession,"SRA.Accession"], " ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("rm ", batch_counter, "_readme.yml\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (3)
previous_step =  "2-dtool_readme_"
current_step = "3-dtool_freeze_"
next_step = "4-dtool_copy_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # freeze the dtool dataset
  script_contents = c(script_contents, paste0("/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh ",sample_list[this_accession,"SRA.Accession"], "\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

irods_uri = "irods:/jic_overflow/rg-daniel-zilberman/"

# Generate the script files for the next step (4)
previous_step =  "3-dtool_freeze_"
current_step = "4-dtool_copy_"
next_step = "5-dtool_tidy_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # commit the dtool dataset to iRODS
  #script_contents = c(script_contents, paste0("/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh ",sample_list[this_accession,"SRA.Accession"], " ", irods_uri, "\n"))
  # Commit the dataset to iRODS and capture the URI
  script_contents = c(script_contents, paste0("OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh ",sample_list[this_accession,"SRA.Accession"], " ", irods_uri, ")\n"))
  # Add the current accession and the iRODS URI to a URIs log file for this batch
  script_contents = c(script_contents, paste0("echo ", sample_list[this_accession,"SRA.Accession"], " $OUTPUT_URI >>", current_step, batch_number, "_uris.txt\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-long # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (5)
previous_step =  "4-dtool_copy_"
current_step = "5-dtool_tidy_"
next_step = "6-align_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # delete the local dtool dataset
  script_contents = c(script_contents, paste0("rm -r ", linux_path, "dtool_datasets/", sample_list[this_accession,"SRA.Accession"], "\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")

    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}


# This is the Col-0 reference
reference_file = "/jic/scratch/groups/Daniel-Zilberman/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta"
reference_id = "TAIR10"
# This is the base location of the 1001 pseudogenomes references
#reference_file = ""
#reference_id=""

alignment_location = paste0(linux_path, "../3-alignments/")

# Generate the script files for the next step (6)
previous_step =  "5-dtool_tidy_"
current_step = "6-align_"
next_step = "7-call_methylation_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
task_array = list()
for (this_accession in 1:nrow(sample_list)) {
  # add the accession to the task array
  #task_array = c(task_array, as.character(sample_list[this_accession,"SRA.Accession"]))
  # abandoned task arrays as they would involve having the whole batch of alignment files on the SSD at the same time (ran out of space for some batches).  Instead, decided to do them serially
  
  batch_counter = batch_counter + 1
  # go through the individual sets of reads (runs) for the accession, and check each one for SINGLE or PAIRED status
  for (this_run in fastq_list[as.character(fastq_list$experiment_accession) == as.character(sample_list[this_accession,"SRA.Accession"]),"run_accession"]) {
    this_library_layout = as.character(fastq_list[fastq_list$run_accession==this_run,"library_layout"])
	if (this_library_layout == "SINGLE") {
      # We are aligning SRRxxxxxxx.fastq.gz to the reference
	  script_contents = c(script_contents, "bsmap -a $BASE_DIR/", this_run, ".fastq.gz -d ", reference_file, " -o $BASE_DIR/", this_run, "_", reference_id, ".bsp\n")
	 #  bsmap -a $R1fq -b $R2fq -d "${rawroot}"/Mp_TGAC_V1.1_scaffolds.fa -o "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TGAC1.1.bsp

	} else if (this_library_layout == "PAIRED") {
      # we are aligning SRRxxxxxxx_1.fastq.gz and SRRxxxxxxx_2.fastq.gz to the reference
	  script_contents = c(script_contents, "bsmap -a $BASE_DIR/", this_run, "_1.fastq.gz -b $BASE_DIR/", this_run, "_2.fastq.gz -d ", reference_file, " -o $BASE_DIR/", this_run, "_", reference_id, ".bsp\n")
    } else {
      # there is no else because this can't happen
	}
	# The run is now aligned to the reference
	
    # once alignment is done, make the CG sites concatenated file of c/t counts (.methratio)
    script_contents = c(script_contents, "python ../Scripts/methratio.py -z -g -o $BASE_DIR/", this_run, "_", reference_id, ".methratio2 -d ", reference_file, " $BASE_DIR/", this_run, "_", reference_id, ".bsp\n")
    # make the per-cytosine file
    #python Scripts/methratio.py -z -o "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TGAC1.1.methratio -d "${rawroot}"/Mp_TGAC_V1.1_scaffolds.fa "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TGAC1.1.bsp

    # filter the .methratio file for just CG sites (takes file size down from 5G to 0.6G so easier to load in R)
    script_contents = c(script_contents, "head -n 1 $BASE_DIR/", this_run, "_", reference_id, ".methratio2 > $BASE_DIR/", this_run, "_", reference_id, ".cg.methratio2\n")
    script_contents = c(script_contents, "awk '{if (substr($4,3,2) == \"CG\") { print } }' $BASE_DIR/", this_run, "_", reference_id, ".methratio2 >> $BASE_DIR/", this_run, "_", reference_id, ".cg.methratio2\n")
    # move the finished gffs (.methratio and .methratio2) to storage
    script_contents = c(script_contents, "mv $BASE_DIR/", this_run, "_", reference_id, ".methratio2 ", alignment_location, "\n")
    script_contents = c(script_contents, "mv $BASE_DIR/", this_run, "_", reference_id, ".cg.methratio2 ", alignment_location, "\n")
    # clean up the temporary space to make space for the next alignment
	script_contents = c(script_contents, "rm $BASE_DIR/", this_run, "_", reference_id, ".*\n")
  }
  
  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch

	# make a slurm header with the right number of array entries for the job array
	# output the contents of task_array as the ARRAY variable in the script
	#slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-long # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=0-",batch_counter - 1,"\n#SBATCH --mem=64000\nARRAY=(", paste(task_array, collapse=" "), ")\n")
    # abandoned task arrays as they would involve having the whole batch of alignment files on the SSD at the same time (ran out of space for some batches).  Instead, decided to do them serially

    # version with no task array
	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-long # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=64000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n#SBATCH --localscratch=ssd:100    # request 100GB SSD\n")

    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\n")
    script_header = paste0(script_header, "BASE_DIR=$SLURM_LOCAL_SCRATCH\n")
	script_header = paste0(script_header, "export DTOOL_CACHE_DIRECTORY=$BASE_DIR\n")
    script_header = paste0(script_header, "#iinit\n")
	script_header = paste0(script_header, "source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")

	script_header = paste0(script_header, "source  bsmap-2.90.0\n")
	script_header = paste0(script_header, "source samtools-1.9\n")

    # add lines at start of script to go through the list of iRODS data sets associated with this batch and retrieve the data files
    script_header = paste0(script_header, "input=dtool_datasets/4-dtool_copy_",batch_number,"_uris.txt\n")
	script_header = paste0(script_header, "while IFS='\\n' read -r line\n")
	script_header = paste0(script_header, "do\n")
	script_header = paste0(script_header, "this_dataset=($line)\n")
	script_header = paste0(script_header, "#echo ${this_dataset[1]}\n")
	script_header = paste0(script_header, "dtool ls ${this_dataset[1]} | while read item; do\n") 
	script_header = paste0(script_header, "echo $item\n")
	script_header = paste0(script_header, "this_item=($item)\n")
	script_header = paste0(script_header, "#echo ${this_item[0]}\n")
	script_header = paste0(script_header, "#echo ${this_item[1]}\n")
	script_header = paste0(script_header, "READ1_ABSPATH=$(dtool item fetch ${this_dataset[1]} ${this_item[0]})\n")
	script_header = paste0(script_header, "#echo $READ1_ABSPATH\n")
	script_header = paste0(script_header, "mv $READ1_ABSPATH $BASE_DIR/${this_item[1]}\n")
	script_header = paste0(script_header, "done\n")	
	script_header = paste0(script_header, "done<\"$input\"\n")
    # we now have all the data files for the batch in SSD storage with their original filenames
	
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	  # start a new task_array
	  task_array = list()
	}
  }
}

# Generate the script files for the next step (7)
previous_step =  "6-align_"
current_step = "7-call_methylation_"
next_step = "8-non_CG_methylation_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # Call R and pass the sample ID to the call_methylation script
  script_contents = c(script_contents, paste0("R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=", sample_list[this_accession,"SRA.Accession"], "\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-medium # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=64000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")

    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\n#cd ",linux_path,"../5-analysis\nsource R_zg-1.0.2\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    # launch the next step in the process
    #launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (8)
previous_step =  "7-call_methylation_"
current_step = "8-non_CG_methylation_"
next_step = "9-dtool_create_mc_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1

  for (this_run in fastq_list[as.character(fastq_list$experiment_accession) == as.character(sample_list[this_accession,"SRA.Accession"]),"run_accession"]) {
    # make a file of the non-cg methylation
    #script_contents = c(script_contents, "head -n 1 $BASE_DIR/", this_run, "_", reference_id, ".methratio2 > $BASE_DIR/", this_run, "_", reference_id, ".non-cg.methratio2\n")
    # don't need to add header, as header does not contain "CG" in relevant position so gets output anyway. Also changed >> to > in line below
    script_contents = c(script_contents, "awk '{if (substr($4,3,2) != \"CG\") { print } }' $BASE_DIR/", this_run, "_", reference_id, ".methratio2 > $BASE_DIR/", this_run, "_", reference_id, ".non-cg.methratio2\n")
  }
  
  # Call R and pass the sample ID to the call_methylation script
  script_contents = c(script_contents, paste0("R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=", sample_list[this_accession,"SRA.Accession"], "\n"))

  # delete the non-cg methylation files as redundant
  for (this_run in fastq_list[as.character(fastq_list$experiment_accession) == as.character(sample_list[this_accession,"SRA.Accession"]),"run_accession"]) {
    script_contents = c(script_contents, "rm $BASE_DIR/", this_run, "_", reference_id, ".non-cg.methratio2\n")
  }
  
  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=64000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")

    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\n#cd ",linux_path,"../5-analysis\nsource R_zg-1.0.2\n")
	
    script_header = paste0(script_header, "BASE_DIR=../3-alignments\n")

	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    # launch the next step in the process
    #launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (9)
previous_step =  "8-non_CG_methylation_"
current_step = "9-dtool_create_mc_"
next_step = "10-dtool_readme_mc_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # make a new dtool_dataset for the accession's methylation data
  script_contents = c(script_contents, paste0("dtool create ",sample_list[this_accession,"SRA.Accession"], "_meth\n"))
  for (this_run in fastq_list[as.character(fastq_list$experiment_accession) == as.character(sample_list[this_accession,"SRA.Accession"]),"run_accession"]) {
    # add each of the files for the accession to the dtool_dataset
	#this_filename = strsplit(this_file,"/")[[1]][[length(strsplit(this_file,"/")[[1]])]]
    script_contents = c(script_contents, paste0("mv ", alignment_location, this_run, "_TAIR10.cg.methratio2 ", linux_path, "dtool_datasets/", sample_list[this_accession, "SRA.Accession"], "_meth/data \n"))	
    script_contents = c(script_contents, paste0("mv ", alignment_location, this_run, "_TAIR10.methratio2 ", linux_path, "dtool_datasets/", sample_list[this_accession, "SRA.Accession"], "_meth/data \n"))	
  }
  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (10)
previous_step =  "9-dtool_create_mc_"
current_step = "10-dtool_readme_mc_"
next_step = "11-dtool_freeze_mc_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # make a new dtool readme for the accession based on the pre-prepared README.yml
  script_contents = c(script_contents, paste0("echo --- > ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("echo description: Bisulphite sequencing methylation coverage data, SRA accession ", sample_list[this_accession,"SRA.Accession"], " >> ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("cat README.yml >> ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("dtool readme write ",sample_list[this_accession,"SRA.Accession"], "_meth ", batch_counter, "_readme.yml\n"))
  script_contents = c(script_contents, paste0("rm ", batch_counter, "_readme.yml\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (11)
previous_step =  "10-dtool_readme_mc_"
current_step = "11-dtool_freeze_mc_"
next_step = "12-dtool_copy_mc_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # freeze the dtool dataset
  script_contents = c(script_contents, paste0("/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh ",sample_list[this_accession,"SRA.Accession"], "_meth\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

irods_uri = "irods:/jic_overflow/rg-daniel-zilberman/"

# Generate the script files for the next step (12)
previous_step =  "11-dtool_freeze_mc_"
current_step = "12-dtool_copy_mc_"
next_step = "13-dtool_tidy_mc_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # commit the dtool dataset to iRODS
  #script_contents = c(script_contents, paste0("/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh ",sample_list[this_accession,"SRA.Accession"], " ", irods_uri, "\n"))
  # Commit the dataset to iRODS and capture the URI
  script_contents = c(script_contents, paste0("OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh ",sample_list[this_accession,"SRA.Accession"], "_meth ", irods_uri, ")\n"))
  # Add the current accession and the iRODS URI to a URIs log file for this batch
  script_contents = c(script_contents, paste0("echo ", sample_list[this_accession,"SRA.Accession"], "_meth $OUTPUT_URI >>", current_step, batch_number, "_uris.txt\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
  	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-long # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")
    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n#iinit\nsource /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Generate the script files for the next step (13)
previous_step =  "12-dtool_copy_mc_"
current_step = "13-dtool_tidy_mc_"
next_step = ""
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1
  # delete the local dtool dataset
  script_contents = c(script_contents, paste0("rm -r ", linux_path, "dtool_datasets/", sample_list[this_accession,"SRA.Accession"], "_meth\n"))

  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-short # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=16000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")

    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\ncd ",linux_path,"dtool_datasets\n")
	# at the end make a flag file indicating the step is complete
    script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\ncd ..\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}

# Step 6a is for retrieving alignments that have already been booked into iRODS
# Generate the script files for the next step (6a)
previous_step =  "6-align_"
current_step = "6a-retrieve_alignments_"
next_step = "8a-segmentation_"
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
task_array = list()
for (this_accession in 1:nrow(sample_list)) {
  # add the accession to the task array
  #task_array = c(task_array, as.character(sample_list[this_accession,"SRA.Accession"]))
  # abandoned task arrays as they would involve having the whole batch of alignment files on the SSD at the same time (ran out of space for some batches).  Instead, decided to do them serially
  
  batch_counter = batch_counter + 1
  # go through the individual sets of reads (runs) for the accession, and check each one for SINGLE or PAIRED status
  
  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch

	# make a slurm header with the right number of array entries for the job array
	# output the contents of task_array as the ARRAY variable in the script
	#slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-long # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=0-",batch_counter - 1,"\n#SBATCH --mem=64000\nARRAY=(", paste(task_array, collapse=" "), ")\n")
    # abandoned task arrays as they would involve having the whole batch of alignment files on the SSD at the same time (ran out of space for some batches).  Instead, decided to do them serially

    # version with no task array
	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-medium # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=64000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n#SBATCH --localscratch=ssd:100    # request 100GB SSD\n")

    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\n")
    script_header = paste0(script_header, "BASE_DIR=$SLURM_LOCAL_SCRATCH\n")
    #script_header = paste0(script_header, "BASE_DIR=../3-alignments\n")
	script_header = paste0(script_header, "export DTOOL_CACHE_DIRECTORY=$BASE_DIR\n")
    script_header = paste0(script_header, "#iinit\n")
	script_header = paste0(script_header, "source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh\n")

    # add lines at start of script to go through the list of iRODS data sets associated with this batch and retrieve the data files
    script_header = paste0(script_header, "input=dtool_datasets/12-dtool_copy_mc_",batch_number,"_uris.txt\n")
	script_header = paste0(script_header, "while IFS='\\n' read -r line\n")
	script_header = paste0(script_header, "do\n")
	script_header = paste0(script_header, "this_dataset=($line)\n")
	script_header = paste0(script_header, "#echo ${this_dataset[1]}\n")
	script_header = paste0(script_header, "dtool ls ${this_dataset[1]} | while read item; do\n") 
	script_header = paste0(script_header, "echo $item\n")
	script_header = paste0(script_header, "this_item=($item)\n")
	script_header = paste0(script_header, "#echo ${this_item[0]}\n")
	script_header = paste0(script_header, "#echo ${this_item[1]}\n")
	script_header = paste0(script_header, "READ1_ABSPATH=$(dtool item fetch ${this_dataset[1]} ${this_item[0]})\n")
	script_header = paste0(script_header, "#echo $READ1_ABSPATH\n")
	script_header = paste0(script_header, "mv $READ1_ABSPATH ../3-alignments/${this_item[1]}\n")
	script_header = paste0(script_header, "done\n")	
	script_header = paste0(script_header, "done<\"$input\"\n")
    # we now have all the data files for the batch in SSD storage with their original filenames
	
	# at the end make a flag file indicating the step is complete
    #script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
    # version which also calls script to retrieve next set in series:
	script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nsbatch ", next_step, batch_number, ".sh\nsbatch ", current_step, batch_number+1, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	  # start a new task_array
	  task_array = list()
	}
  }
}

# 8a is the step for segmenting the methylome
# Generate the script files for the next step (8a)
previous_step =  "6a-retrieve_alignments_"
current_step = "8a-segmentation_"
next_step = ""
batch_counter = 0
batch_number = 0
# start a new batch
script_contents = list()
for (this_accession in 1:nrow(sample_list)) {
  batch_counter = batch_counter + 1

  for (this_run in fastq_list[as.character(fastq_list$experiment_accession) == as.character(sample_list[this_accession,"SRA.Accession"]),"run_accession"]) {
    # make a file of the non-cg methylation
    #script_contents = c(script_contents, "head -n 1 $BASE_DIR/", this_run, "_", reference_id, ".methratio2 > $BASE_DIR/", this_run, "_", reference_id, ".non-cg.methratio2\n")
    # don't need to add header, as header does not contain "CG" in relevant position so gets output anyway. Also changed >> to > in line below
    script_contents = c(script_contents, "awk '{if (substr($4,3,2) != \"CG\") { print } }' $BASE_DIR/", this_run, "_", reference_id, ".methratio2 > $BASE_DIR/", this_run, "_", reference_id, ".non-cg.methratio2\n")
  }
  
  # Call R and pass the sample ID to the call_methylation script
  script_contents = c(script_contents, paste0("mkdir ../5-analysis/", sample_list[this_accession,"SRA.Accession"], "\n"))
  script_contents = c(script_contents, paste0("R_zg <../Scripts/8a-non_CG_methylation.R --no-save --args --sample=", sample_list[this_accession,"SRA.Accession"], "\n"))

  # delete the non-cg methylation files as redundant
  for (this_run in fastq_list[as.character(fastq_list$experiment_accession) == as.character(sample_list[this_accession,"SRA.Accession"]),"run_accession"]) {
    script_contents = c(script_contents, "rm $BASE_DIR/", this_run, "_", reference_id, ".non-cg.methratio2\n")
  }
  
  if ((batch_counter==batch_size) | (this_accession==nrow(sample_list))) {
    # finish batch
	slurm_header = paste0("#!/bin/bash -e\n#SBATCH -p jic-medium # partition (queue)\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address\n#SBATCH --array=\n#SBATCH --mem=64000\n#SBATCH -o ", current_step, batch_number, ".out\n#SBATCH -e ", current_step, batch_number, ".err\n")

    # at the start check to ensure flag file exists for previous step, and flag file doesn't exist for current step
    script_header = paste0(slurm_header, "until [ -f ", linux_path, previous_step, batch_number, ".txt ]; do\n  sleep 60\ndone\nif [ ! -f ", linux_path, current_step, batch_number, ".txt ]; then\n#cd ",linux_path,"../5-analysis\nsource R_zg-1.0.2\n")
	
    script_header = paste0(script_header, "BASE_DIR=../3-alignments\n")

	# at the end make a flag file indicating the step is complete
    #script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nsbatch ", next_step, batch_number, ".sh\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
    # version which does not call next step:
	script_footer = paste0("echo complete>", linux_path, current_step, batch_number, ".txt\necho ", current_step, batch_number, " complete>>", log_file, "\nelse\necho ", current_step, batch_number, " already_complete>>", log_file, "\nfi\n")
	script_contents = c(script_header, script_contents, script_footer)
	write.table(script_contents, file=paste0(path,current_step,batch_number,".sh"), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
    # launch the next step in the process
    #launch_script_contents = c(launch_script_contents, "sbatch ", paste0(path,current_step,batch_number,".sh\n"))
	if (this_accession<nrow(sample_list)) {
	  # start new batch
	  script_contents = list()
	  batch_number = batch_number + 1
	  batch_counter = 0
	}
  }
}


write.table(launch_script_contents, file=paste0(path,launch_script), sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)

# After scripts have been generated in windows using this script, they need to be converted for use on the cluster. Run these commands at an interactive prompt:
perl -pi -e 's/\r\n/\n/g' 0-wget_batch_*.sh
perl -pi -e 's/\r\n/\n/g' 1-dtool_create_*.sh
perl -pi -e 's/\r\n/\n/g' 2-dtool_readme_*.sh
perl -pi -e 's/\r\n/\n/g' 3-dtool_freeze_*.sh
perl -pi -e 's/\r\n/\n/g' 4-dtool_copy_*.sh
perl -pi -e 's/\r\n/\n/g' 5-dtool_tidy_*.sh
perl -pi -e 's/\r\n/\n/g' 6-align_*.sh
perl -pi -e 's/\r\n/\n/g' 7-call_methylation_*.sh
perl -pi -e 's/\r\n/\n/g' 8-non_CG_methylation_*.sh
perl -pi -e 's/\r\n/\n/g' 9-dtool_create_mc_*.sh
perl -pi -e 's/\r\n/\n/g' 10-dtool_readme_mc_*.sh
perl -pi -e 's/\r\n/\n/g' 11-dtool_freeze_mc_*.sh
perl -pi -e 's/\r\n/\n/g' 12-dtool_copy_mc_*.sh
perl -pi -e 's/\r\n/\n/g' 13-dtool_tidy_mc_*.sh
perl -pi -e 's/\r\n/\n/g' launch_script.sh
perl -pi -e 's/\r\n/\n/g' 6a-retrieve_alignments_*.sh
perl -pi -e 's/\r\n/\n/g' 8a-segmentation_*.sh

# Count SRX datasets in iRODS
dtool ls irods:/jic_overflow/rg-daniel-zilberman|grep SRX|wc -l

