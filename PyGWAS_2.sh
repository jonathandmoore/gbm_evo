#!/bin/bash -e
#SBATCH -p jic-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jonathan.moore@jic.ac.uk # send-to address
#SBATCH --array=2-9
#SBATCH --mem=32000

#source /jic/software/staging/CISSUPPORT-11453/stagingloader
#source pygwas-1.7.4
source package aceb215e-0972-41bf-b1bd-47da4fa9aa94

base_dir=/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes
ref_dir=${base_dir}/0-reference/PyGWAS
pheno_dir=${base_dir}/pygwas/phenotype

#for threshold in 8 10 12 14
#for threshold in 8
#do

#results_dir=${base_dir}/pygwas/results
#results_dir=${base_dir}/pygwas/results/8
#results_dir=${base_dir}/pygwas/results/10
#results_dir=${base_dir}/pygwas/results/12
#results_dir=${base_dir}/pygwas/results/14
#results_dir=${base_dir}/pygwas/results/${threshold}
#results_dir=${base_dir}/pygwas/results/412
#results_dir=${base_dir}/pygwas/results/412_mCG
#results_dir=${base_dir}/pygwas/results/625
#results_dir=${base_dir}/pygwas/results/625_mCG
#results_dir=${base_dir}/pygwas/results/Spain_Germany
#results_dir=${base_dir}/pygwas/results/Spain_Germany2
#results_dir=${base_dir}/pygwas/results/mcg_k-means_9_cluster_centres
#results_dir=${base_dir}/pygwas/results/flowering_related_traits
results_dir=${base_dir}/pygwas/results/boxcox_transformed_rfitness

#pheno_file=expression_fitness_associated_genes_fitness_accessions_panel_412.csv
#pheno_file=mCG_fitness_associated_genes_fitness_accessions_panel_412.csv
#pheno_file=expression_fitness_associated_genes_625_accessions_panel.csv
#pheno_file=mCG_fitness_associated_genes_625_accessions_panel.csv
#pheno_file=Bonafide_gbm_948_0.60_accession.csv
#pheno_file=Bonafide_gbm_948_all_accession.csv
#pheno_file=20_genes_mCG_gene_body_trial_for_pyGWAS.csv
#pheno_file=rfitness_arabidopsis_pehnotypes_Spain_Germany.csv
#pheno_file=mcg_k-means_9_cluster_centres.csv

#pheno_file=rfiteness412accessoinsspainandgermany.csv
#pheno_file=FLowering_related_traits_phenotypes_for_pyGWAS.csv
pheno_file=boxcox_transformed_rfitness_data_for_pyGWAS.csv

#pheno_file=gbm_${threshold}_fold_coverage_Chr01_top500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_501_to_2500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_2501_to_4500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_4501_to_6500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_6501_to_8500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_8501_to_10500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_10501_to_12500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_12501_to_14500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_14501_to_16500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_16501_to_18500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_18501_to_20500_genes.csv
#pheno_file=bona_fide_gbm_${threshold}x_20501_to_22349_genes.csv



# set filter conditions for output files

max_pval=0.0001
min_macs=15

cd ${results_dir}

#for i in {2..21}
#do
  # make file containing the relevant gene phenotype data
  cut -d ',' -f 1,${SLURM_ARRAY_TASK_ID} ${pheno_dir}/${pheno_file} > temp_pheno_${SLURM_ARRAY_TASK_ID}.csv
  gene=`head temp_pheno_${SLURM_ARRAY_TASK_ID}.csv -n 1|cut -d ',' -f 2`
  mv temp_pheno_${SLURM_ARRAY_TASK_ID}.csv ${gene}.csv
  #rm temp_pheno_${SLURM_ARRAY_TASK_ID}.csv
 
  # run pygwas for the gene
  srun pygwas run -t none -a amm -g ${ref_dir}/1001_SNP_MATRIX -k ${ref_dir}/1001_SNP_MATRIX/kinship_ibs_mac5.hdf5 -o ${gene}_results_amm.csv ${gene}.csv
  srun pygwas run -t none -a lm -g ${ref_dir}/1001_SNP_MATRIX -k ${ref_dir}/1001_SNP_MATRIX/kinship_ibs_mac5.hdf5 -o ${gene}_results_lm.csv ${gene}.csv
  
#  head -n 1 ${gene}_results_amm.csv > ${gene}_results_amm_filtered.csv
#  head -n 1 ${gene}_results_lm.csv > ${gene}_results_lm_filtered.csv
#  awk -F "," -v max_pval=${max_pval} '($3 + 0) <= max_pval' ${gene}_results_amm.csv | awk -F "," -v min_macs=${min_macs} '($5 + 0) >= min_macs' >> ${gene}_results_amm_filtered.csv
#  awk -F "," -v max_pval=${max_pval} '($3 + 0) <= max_pval' ${gene}_results_lm.csv | awk -F "," -v min_macs=${min_macs} '($5 + 0) >= min_macs' >> ${gene}_results_lm_filtered.csv

  # This version filters only on min_macs but returns all p-values
  head -n 1 ${gene}_results_amm.csv > ${gene}_results_amm_filtered_macs_${min_macs}.csv
  head -n 1 ${gene}_results_lm.csv > ${gene}_results_lm_filtered_macs_${min_macs}.csv
  awk -F "," -v min_macs=${min_macs} '($5 + 0) >= min_macs' ${gene}_results_amm.csv >> ${gene}_results_amm_filtered_macs_${min_macs}.csv
  awk -F "," -v min_macs=${min_macs} '($5 + 0) >= min_macs' ${gene}_results_lm.csv >> ${gene}_results_lm_filtered_macs_${min_macs}.csv
  
  #  rm ${gene}_results.csv
#done

#done

#srun tar cvzf AT1G11XXXX.tar.gz AT1G11*_results.csv
#mv AT1G11XXXX.tar.gz ~
#rm AT1G11*_results.csv
  
