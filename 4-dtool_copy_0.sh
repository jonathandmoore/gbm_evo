#!/bin/bash -e
#SBATCH -p jic-long # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=16000
#SBATCH -o 4-dtool_copy_0.out
#SBATCH -e 4-dtool_copy_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/3-dtool_freeze_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/4-dtool_copy_0.txt ]; then
cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190716 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190716 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190717 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190717 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190718 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190718 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190719 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190719 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190720 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190720 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190721 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190721 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190722 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190722 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190723 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190723 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190724 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190724 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190725 irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190725 $OUTPUT_URI >>4-dtool_copy_0_uris.txt
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/4-dtool_copy_0.txt
echo 4-dtool_copy_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
cd ..
sbatch 5-dtool_tidy_0.sh
else
echo 4-dtool_copy_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

