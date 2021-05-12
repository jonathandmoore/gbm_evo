# gbm_evo
**Gene body methylation mediates epigenetic inheritance of plant traits**: the scripts

Our paper shows evidence for the role of gene body methylation in shaping the evolution of plants' responses to their environment.  
https://www.biorxiv.org/content/10.1101/2021.03.15.435374v1

The main analysis scripts underlying the work are shared here in the spirit of transparency and reproducibility.

Please note that these scripts are not intended as a ready-to-run software package. They are the original scripts used for our analysis, and are only likely to function in the specific software and data environment they were designed for, without re-engineering.  Some parts of some scripts implement experimental analyses, which were not ultimately included in the published work, and some parts are not executable.

**Primary analysis/data reduction:**

**ENA_SR....txt** and **GSE....txt** files define the set of relevant ENA/SRA accessions.

**0-ETL_1001_methylomes.R** is the Extract/Translate/Load script that starts the process of analysing the 1001 methylomes data set.  It is an R script whose purpose is to create shell scripts to use a SLURM cluster to carry out a batch-wise methylation analysis of the 1001 methylomes data sets, downloaded from SRA/ENA.

This script reads the lists of samples and FASTQ files, downloaded from ENA and SRA. It associates FASTQ files with samples, and samples with Arabidopsis accession identifiers from the 1001 genomes project.  It then works through the samples, ten at a time, making SLURM shell scripts to carry out each subsequent stage of the primary analysis for the 10 samples in question:

  - download the relevant FASTQ files of bisulphite sequence data from SRA (wget)
  - check the FASTQ files into an object store (dtool)
  - retrieve them to local cluster node temporary SSD (dtool)
  - align them to the genome reference and merge the results from multiple FASTQ files per sample (bsmap)
  - call methylation status at each CG site (7-call_methylation.R)
  - segment each methylome into contiguous UM, gbM, teM, or gbM-like segments (8b-analyse_segmentation.R)
  - process non-CG methylation (8-non_CG_methylation.R - no analyses implemented)
  - check the results back into the object store (dtool)

This workflow is carried out batch-wise to limit primary disk space occupied by raw FASTQ and alignment files at any given time, and if run one batch at a time takes up to 100 days elapsed.

Example scripts of each of the stages (0 to 13), for a single batch of 10 samples, are given here (e.g. **0-wget_batch_0.sh**, **1-dtool_create_0.sh**, **2-dtool_readme_0.sh**, **3-dtool_freeze_0.sh**, 4-dtool_copy_0.sh** etc). The ETL script causes the scripts to be daisychained, so as each stage completes it calls the script which initiates the subsequent stage. 

Once the primary data reduction is complete, the remaining R scripts relate to the secondary analyses of methylation, the results of which form the basis of our paper.

**Secondary analysis:**

**parse_1001genome_VCF.sh**, **parse_1001genome_VCF.pl** parse the 1001 genome SNP file to extract lists of SNPs and short indels w.r.t. Col-0 at cytosines and CG sites in each accession.

**7-call_methylation.R** calls the methylation status at each CG site in a sample as M(ethylated), U(nmethylated), P(artially methylated), I(ndeterminate) or X(excised due to presence of SNP/indel). Fisher's exact test is used to check for sufficient coverage, taking into account bisulphite conversion rates in the chloroplast of the sample in question, and for sites with sufficient coverage binomial test is used to compare numbers of reads indicating M or U at the locus.

**8b-analyse_segmentation.R** declares and executes functions derived from a modified version of the MethylSeqR package, adapted to analyse plant methylomes. The script segments the raw CG and non-CG methylomes for a sample to generate a model of the genome partitioned into contiguous UM, gbM, teM and gbM-like segments. Segments are subsequently adjusted to reflect the genome annotation (for gbM segments), and missing data, and segments are trimmed so that only those supported by sufficient sites with reliable methylation calls are relied upon for further analysis.

Methylation calls and methylation segmentation models are integrated, to generate methylation status calls and methylation level estimates for each gene in each accession:

**10-consolidate_methylation_calls.R** makes a single table with methylation calls at each site in all accessions, filtered for CG sites with SNPs, and estimates mean mCG per gene per sample based on these calls.

**classify_genes_methylation_status.R** uses the methylome segmentation model and CG site methylation state in each accession to classify each gene's methylation status in each accession as UM, gbM, teM, or both, or indeterminate. Methylation calls at individual sites are used as confirmatory evidence of results from segmentation of raw methylomes. A variety of sensitivity thresholds are implemented, so that we could test how performance of epiGWA varied depending on how genes' methylation status was defined. Variant 2 was used ultimately in the paper: gene is classed as gbM in an accession if it is overlapped by a gbM segment in the segmentation model of raw methylation data, but not substantially by a teM segment, and the gbM segment is supported by at least one site with a methylation call of 'M' by statistical analysis of read coverage at the site.

**correl_genes_total_mcg.R** correlates mCG by gene with genome-wide mean mCG. 

**maps_code.R** plots locations of accessions on outline maps of Europe.

**pygwas_3.sh** is an example of running the PyGWAS application on SLURM cluster to estimate associations between loci and phenotypes. Various models (lm and amm), and filtering on p-value and MAF thresholds are implemented.

**consolidate_pygwas_results.R** and **consolidate_pygwas_p-values.R** process PyGWAS results for consolidated reporting. A number of alternative minimum coverage thresholds are implemented.

**run_tassel.sh** is an example of running tassel on SLURM cluster for estimating associations between ~200 environmental variables and gene methylation status of individual genes. A set of relevant example config files are provided: 881_mlm_config_template.xml, popstructure881acc_7pca_commandlineformat.txt, gbm_881_acc_k_matrix.txt, 881accessions_gbm_umr_calls_2_polymorphiclocifortassel.txt, 200_environments_for_pyGWAS.txt

**parse_haplo_data.R** loads up SNP data, creates a haplogroups file for each gene, correlates gbM and teM at retained eQTL by haplogroup, and writes out matrices of p-values etc.

**merge_segmentation_models.R** merges segmentation models across accessions to identify the broadly methylatable gene body space in Arabidopsis. A wide range of parameter and sensitivity options are implemented.

**met1_vs_col0.R** generates comparisons of bisulphite-seq and H3K9me2 chipseq data for met1 vs. Col-0.

**boxcox_eQTLs.R** applies boxcox transform to gbM and teM eQTLs for association analysis of trans-acting factors.

**extract_gbm_example_states.R** replaces a section from the middle of merge_segmentation_models.R. The purpose of the section is to extract example accession CG methylomes from each decile of global genome mean mCG, for modelling and analysis.

**ends_analysis.pl** analyses two gff files, one for 'genes' or regions of interest, and one for scores at loci. It generates a metaplot histogram spanning the ends of the 'genes' showing how the mean score averaged over all genes varies at bp resolution around the borders of 'genes'.  It is copied from DZLab_tools from zilbermanlab.net, with one fix applied at line 214. 

**parse_h2az_data.R** loads Chip-Seq or CHIP-Chip data relating to H2A.Z from published experiments, parses and smooths, and uses the H2A.Z abundance at the boundaries of annotated gbM segments to trim those segments before modelling. Also filters segments based on their remaining mean H2A.Z, and allows to remove segments from genes with teM calls > 0.05 MAF. This script also generates .GFF data files and perl commands to run ends_analysis.pl to produce annotation metaplots of filtered segments, with respect to H2A.Z, CG site density and methylation density, in order to assess model performance.

**parse_h2az_associations.R** loads lists of genes associated by PyGWAS with variation around HTA9 gene (H2A.Z), and compares these genes with those not associated, in terms of their H2A.Z.


**Scripts specific to the analysis of the Schmitz and Becker MA line data sets**
  
  - Schmitz et al, 2011 - NCBI SRA accession SRA035939
  - Becker et al, 2011 - EMBL ENA accession PRJEB2678

  - assign_gene_methylation_status_2020-08-21.R
  - 5_estimate_de_novo_rate_2019-09-12.R
  - 5-estimate_GBM_rates_2020-08-07.R
  - 5-estimate_GBM_boundary_size_2019-04-25.R
  - 5-methylation_calls_2018-06-26.R
  - 5-methylation_calls_optimise_non_CG_segmentationf_2018-03-14.R
  - 5-methylation_calls_optimise_mCG_cutoff_2018-03-14.R
  - 5-methylation_calls.sh
  - 4-find_CGs.pl
  - 4-find_CGs.sh
  - 4-generate_covs.pl
  - 4-annotate_gffs.pl
  - 4-annotate_gffs.sh
  - 3-align_bs-sequel_PRJEB2678.sh
  - 3-align_bs-sequel_SRA035939.sh
  - 3-align_SRA035939.sh
  - 3-align-PRJEB2678.sh
  - 2-trim_SRA035939.sh
  - 2-trim_PRJEB2678.sh
  - 1-qc.sh
  - 0-stage.sh


**Phylogenetic analyses**
  - 5-RAxML_2018-07-30.bat
  - 5-phylogenetic_analysis_2019-08-19.R
