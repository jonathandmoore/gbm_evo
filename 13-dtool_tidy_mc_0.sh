#!/bin/bash -e
#SBATCH -p jic-short # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=16000
#SBATCH -o 13-dtool_tidy_mc_0.out
#SBATCH -e 13-dtool_tidy_mc_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/12-dtool_copy_mc_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/13-dtool_tidy_mc_0.txt ]; then
cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190716_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190717_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190718_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190719_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190720_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190721_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190722_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190723_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190724_meth
rm -r /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets/SRX2190725_meth
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/13-dtool_tidy_mc_0.txt
echo 13-dtool_tidy_mc_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
cd ..
sbatch 0.sh
else
echo 13-dtool_tidy_mc_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

