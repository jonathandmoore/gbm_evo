#!/bin/bash -e
#SBATCH -p jic-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=64000
#SBATCH -o 7-call_methylation_0.out
#SBATCH -e 7-call_methylation_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/6-align_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/7-call_methylation_0.txt ]; then
#cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../5-analysis
source R_zg-1.0.1
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190716
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190717
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190718
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190719
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190720
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190721
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190722
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190723
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190724
R_zg <../Scripts/7-call_methylation.R --no-save --args --sample=SRX2190725
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/7-call_methylation_0.txt
echo 7-call_methylation_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
sbatch 8-non_CG_methylation_0.sh
else
echo 7-call_methylation_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

