#!/bin/bash -e
#SBATCH -p jic-long # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=16000
#SBATCH -o 12-dtool_copy_mc_0.out
#SBATCH -e 12-dtool_copy_mc_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/11-dtool_freeze_mc_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/12-dtool_copy_mc_0.txt ]; then
cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190716_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190716_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190717_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190717_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190718_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190718_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190719_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190719_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190720_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190720_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190721_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190721_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190722_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190722_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190723_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190723_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190724_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190724_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
OUTPUT_URI=$(/jic/scratch/groups/Daniel-Zilberman/dtool-slurm-scripts/dtool_copy.slurm.sh SRX2190725_meth irods:/jic_overflow/rg-daniel-zilberman/)
echo SRX2190725_meth $OUTPUT_URI >>12-dtool_copy_mc_0_uris.txt
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/12-dtool_copy_mc_0.txt
echo 12-dtool_copy_mc_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
cd ..
sbatch 13-dtool_tidy_mc_0.sh
else
echo 12-dtool_copy_mc_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

