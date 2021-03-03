#!/bin/bash -e
#SBATCH -p jic-short # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=16000
#SBATCH -o 3-dtool_freeze_0.out
#SBATCH -e 3-dtool_freeze_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/2-dtool_readme_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/3-dtool_freeze_0.txt ]; then
cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190716
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190717
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190718
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190719
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190720
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190721
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190722
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190723
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190724
/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_freeze.slurm.sh SRX2190725
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/3-dtool_freeze_0.txt
echo 3-dtool_freeze_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
cd ..
sbatch 4-dtool_copy_0.sh
else
echo 3-dtool_freeze_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

