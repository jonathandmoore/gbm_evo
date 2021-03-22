# gbm_evo
**Gene body methylation mediates epigenetic inheritance of plant traits**: the scripts

Our paper shows evidence for the role of gene body methylation in shaping the evolution of plants' responses to their environment.  
https://www.biorxiv.org/content/10.1101/2021.03.15.435374v1

The analysis scripts underlying the work are shared here in the spirit of transparency and reproducibility.

Please note that these scripts are not intended as a ready-to-run software package. They are the original scripts used for our analysis, and are only likely to function in the specific software and data environment they were designed for, without re-engineering.

**Primary analysis/data reduction:**

**0-ETL_1001_methylomes.R** is the Extract/Translate/Load script that starts the process of analysing the 1001 methylomes data set.  It is an R script whose purpose is to create shell scripts to use a SLURM cluster to carry out a batch-wise methylation analysis of the 1001 methylomes data sets, downloaded from SRA/ENA.

This script reads the lists of samples and FASTQ files, downloaded from ENA and SRA (**ENA_SR*.txt** and **GSE*.txt** files define the set of relevant ENA/SRA accessions). It associates FASTQ files with samples, and samples with Arabidopsis accession identifiers from the 1001 genomes project.  It then works through the samples, ten at a time, making SLURM shell scripts to carry out each subsequent stage of the primary analysis for the 10 samples in question:

  - download the relevant FASTQ files of bisulphite sequence data from SRA (wget)
  - check the FASTQ files into an object store (dtool)
  - retrieve them to local cluster node temporary SSD (dtool)
  - align them to the genome reference and merge the results from multiple FASTQ files per sample (bsmap)
  - call methylation status at each CG site (7-call_methylation.R)
  - segment each methylome into contiguous UM, gbM, teM, or gbM-like segments (8b-analyse_segmentation.R)
  - process the non-CG methylation
  - check the results back into the object store

This workflow is carried out batch-wise to limit primary disk space occupied by raw FASTQ and alignment files at any given time, and if run one batch at a time takes up to 100 days elapsed.

Example scripts of each of the stages (0 to 13), for a single batch of 10 samples, are given here (e.g. **0-wget_batch_0.sh**, **1-dtool_create_0.sh**, **2-dtool_readme_0.sh**, **3-dtool_freeze_0.sh**, 4-dtool_copy_0.sh** etc). The ETL script causes the scripts to be daisychained, so as each stage completes it calls the script which initiates the subsequent stage. 

Once the primary data reduction is complete, the remaining R scripts relate to the secondary analyses of methylation, which forms the basis of our paper.

**Secondary analysis:**

**parse_1001genome_VCF.sh**, **parse_1001genome_VCF.pl** parse the 1001 genome SNP file to extract list of SNPs and short indels w.r.t. Col-0 in each accession.

**7-call_methylation.R** calls the methylation status of each CG site in a sample as M(ethylated), U(nmethylated), P(artially methylated), I(ndeterminate), X(excised due to presence of SNP/indel). Fisher's exact test is used to check for sufficient coverage, taking into account bisulphite conversion rates in the chloroplast of the sample in question, and for sites with suficient coverage binomial test is used to compare number of reads indicating M or U at the locus.

**8b-analyse_segmentation.R** segments the raw CG and non-CG methylomes for a sample to generate a model of the genome partitioned into contiguous UM, gbM, teM and gbM-like segments. Segments are subsequently adjusted so that only those supported by sufficient sites with reliable methylation calls are relied upon for further analysis.

Methylation calls and methylation segmentations are integrated, to generate methylation status calls and methylation level estimates for each gene in each accession:

**10-consolidate_methylation_calls.R** makes a single table with methylation calls at each site in all accessions, filtered for CG sites with SNPs, and estimates mean mCG per gene per sample based on calls.

**classify_genes_methylation_status.R** uses the segmentation and CG site methylation state in each accession to classify each gene's methylation status in each accession as UM, gbM, teM, or both, or indeterminate. A variety of sensitivity thresholds are implemented.

**maps_code.R** plots locations of accessions on outline maps of Europe.

**pygwas_3.sh** is an example of running the PyGWAS application on a SLURM cluster to estimate associations between loci and phenotypes. Various models (lm and amm), and filtering on p-value and MAF thresholds are implemented.

**consolidate_pygwas_results.R** and **consolidate_pygwas_p-values.R** process PyGWAS results for consolidated reporting. A number of alternative minimum coverage thresholds are implemented.

**run_tassel.sh** is an example of running tassel on SLURM cluster for estimating associations between ~200 environmental variables and gene methylation status of individual genes. A set of relevant example config files are provided: 881_mlm_config_template.xml, popstructure881acc_7pca_commandlineformat.txt, gbm_881_acc_k_matrix.txt, 881accessions_gbm_umr_calls_2_polymorphiclocifortassel.txt, 200_environments_for_pyGWAS.txt

**merge_segmentation_models.R** merges segmentation models across accessions to identify the broadly methylatable gene body space in Arabidopsis. A range of parameter and sensitivity options are implemented.

**parse_haplo_data.R** loads up SNP data, creates a haplogroups file for each gene, correlates gbM and teM retained eQTL by haplogroup, and writes out matrices of p-values etc.

**correl_genes_total_mcg.R** correlates mCG by gene with genome-wide mean mCG. 
