# gbm_evo
**Gene body methylation mediates epigenetic inheritance of plant traits**: the scripts

Our paper shows evidence for the role of gene body methylation in shaping the evolution of plants' responses to their environment.  

The analysis scripts underlying the work are archived here for transparency and reproducibility.


**Primary analysis/data reduction:**

**0-ETL_1001_methylomes.R** is the ETL script that starts the process of analysing the 1001 methylomes data set.  It is an R script whose purpose is to create shell scripts to use a SLURM cluster to carry out a batch-wise methylation analysis of the 1001 methylomes data sets, downloaded from SRA.

This script reads the lists of samples and fastq files, downloaded from ENA and SRA (**ENA_SR*.txt** and **GSE*.txt**). It associates FASTQ files with samples, and samples with Arabidopsis accession identifiers from the 1001 genomes project.  It then works through the accessions, ten at a time, making SLURM shell scripts to carry out each subsequent stage of the primary analysis for the 10 accessions in question:

  . download the relevant fastq files from SRA
  . check the FASTQ files into an object store
  . retrieve them to local cluster node temporary SSD
  . align them to the reference
  . call methylation at each CG site
  . segment the methylomes
  . process the non-CG methylation
  . check the results back into the object store

The process is carried out batch-wise to limit primary disk space occupied by raw FASTQ and alignment files at any given time, and if run one batch at a time takes up to 100 days.

Example scripts of each of the stages are given (e.g. **0-wget_batch_0.sh**, **1-dtool_create_0.txt** etc) and the scripts are daisychained, so each stage calls the script to initiate the subsequent stage. 

Once the primary data reduction is complete, the remaining R scripts relate to the secondary analyses of methylation from the paper.

**Secondary analysis:**

Methylation calls, and methylation segmentations, are integrated to generate methylation status calls, and methylation level estimates, for each gene in each accession.

**parse_1001genome_VCF.pl** parses the 1001 genome SNP file to extract list of SNPs w.r.t. Col-0 in each accession.

**8b-analyse_segmentation.R** segments the CG and non-CG methylomes for a sample to generate a model of the genome partitioned into contiguous UM, gbM, teM and gbM-like segments.

**10-consolidate_methylation_calls.R** makes a single table with methylation calls at each site in all accessions, filtered for CG sites with SNPs, then estimates mean mCG per gene per sample.

**maps_code.R** plots locations of accessions on outline maps of Europe.

**pygwas.sh** runs the PyGWAS application on a SLURM cluster to estimate associations between loci and phenotypes.

**consolidate_pygwas_results.R** and **consolidate_pygwas_p-values.R** process PyGWAS results for consolidated reporting. We tried a number of alternative minimum coverage thresholds.

**merge_segmentation_models.R** merges segmentation models across accessions to identify the methylatable space in Arabidopsis. A range of parameter options are tried.

**parse_haplo_data.R** load up SNP data, create a haplogroups file for each gene, correlate gbM and teM retained eQTL and write out matrices of p-values etc.
