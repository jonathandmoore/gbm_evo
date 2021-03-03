#!/bin/bash -e
#SBATCH -p jic-short # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=16000
#SBATCH -o 11-dtool_freeze_mc_0.out
#SBATCH -e 11-dtool_freeze_mc_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/10-dtool_readme_mc_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/11-dtool_freeze_mc_0.txt ]; then
cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190716_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190717_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190718_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190719_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190720_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190721_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190722_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190723_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190724_meth
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190725_meth
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/11-dtool_freeze_mc_0.txt
echo 11-dtool_freeze_mc_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
cd ..
sbatch 12-dtool_copy_mc_0.sh
else
echo 11-dtool_freeze_mc_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

