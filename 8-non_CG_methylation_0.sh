#!/bin/bash -e
#SBATCH -p jic-short # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=64000
#SBATCH -o 8-non_CG_methylation_0.out
#SBATCH -e 8-non_CG_methylation_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/7-call_methylation_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/8-non_CG_methylation_0.txt ]; then
#cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../5-analysis
source R_zg-1.0.1
BASE_DIR=../3-alignments
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295279_TAIR10.methratio2 > $BASE_DIR/SRR4295279_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190716
rm $BASE_DIR/SRR4295279_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295280_TAIR10.methratio2 > $BASE_DIR/SRR4295280_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295281_TAIR10.methratio2 > $BASE_DIR/SRR4295281_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190717
rm $BASE_DIR/SRR4295280_TAIR10.non-cg.methratio2
rm $BASE_DIR/SRR4295281_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295282_TAIR10.methratio2 > $BASE_DIR/SRR4295282_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295283_TAIR10.methratio2 > $BASE_DIR/SRR4295283_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190718
rm $BASE_DIR/SRR4295282_TAIR10.non-cg.methratio2
rm $BASE_DIR/SRR4295283_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295284_TAIR10.methratio2 > $BASE_DIR/SRR4295284_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295285_TAIR10.methratio2 > $BASE_DIR/SRR4295285_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190719
rm $BASE_DIR/SRR4295284_TAIR10.non-cg.methratio2
rm $BASE_DIR/SRR4295285_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295286_TAIR10.methratio2 > $BASE_DIR/SRR4295286_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190720
rm $BASE_DIR/SRR4295286_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295287_TAIR10.methratio2 > $BASE_DIR/SRR4295287_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190721
rm $BASE_DIR/SRR4295287_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295288_TAIR10.methratio2 > $BASE_DIR/SRR4295288_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295289_TAIR10.methratio2 > $BASE_DIR/SRR4295289_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190722
rm $BASE_DIR/SRR4295288_TAIR10.non-cg.methratio2
rm $BASE_DIR/SRR4295289_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295290_TAIR10.methratio2 > $BASE_DIR/SRR4295290_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190723
rm $BASE_DIR/SRR4295290_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295291_TAIR10.methratio2 > $BASE_DIR/SRR4295291_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190724
rm $BASE_DIR/SRR4295291_TAIR10.non-cg.methratio2
awk '{if (substr($4,3,2) != "CG") { print } }' $BASE_DIR/SRR4295292_TAIR10.methratio2 > $BASE_DIR/SRR4295292_TAIR10.non-cg.methratio2
R_zg <../Scripts/8-non_CG_methylation.R --no-save --args --sample=SRX2190725
rm $BASE_DIR/SRR4295292_TAIR10.non-cg.methratio2
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/8-non_CG_methylation_0.txt
echo 8-non_CG_methylation_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
sbatch 9-dtool_create_mc_0.sh
else
echo 8-non_CG_methylation_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

