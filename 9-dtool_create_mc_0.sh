#!/bin/bash -e
#SBATCH -p jic-short # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=16000
#SBATCH -o 9-dtool_create_mc_0.out
#SBATCH -e 9-dtool_create_mc_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/8-non_CG_methylation_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/9-dtool_create_mc_0.txt ]; then
cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
dtool create SRX2190716_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295279_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190716_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295279_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190716_meth/data 
dtool create SRX2190717_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295280_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190717_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295280_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190717_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295281_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190717_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295281_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190717_meth/data 
dtool create SRX2190718_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295282_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190718_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295282_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190718_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295283_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190718_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295283_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190718_meth/data 
dtool create SRX2190719_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295284_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190719_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295284_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190719_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295285_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190719_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295285_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190719_meth/data 
dtool create SRX2190720_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295286_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190720_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295286_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190720_meth/data 
dtool create SRX2190721_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295287_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190721_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295287_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190721_meth/data 
dtool create SRX2190722_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295288_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190722_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295288_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190722_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295289_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190722_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295289_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190722_meth/data 
dtool create SRX2190723_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295290_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190723_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295290_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190723_meth/data 
dtool create SRX2190724_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295291_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190724_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295291_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190724_meth/data 
dtool create SRX2190725_meth
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295292_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190725_meth/data 
mv /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/SRR4295292_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190725_meth/data 
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/9-dtool_create_mc_0.txt
echo 9-dtool_create_mc_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
cd ..
sbatch 10-dtool_readme_mc_0.sh
else
echo 9-dtool_create_mc_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

