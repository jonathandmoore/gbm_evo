#!/bin/bash -e
#SBATCH -p jic-short # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=16000
#SBATCH -o 2-dtool_readme_0.out
#SBATCH -e 2-dtool_readme_0.err
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/1-dtool_create_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/2-dtool_readme_0.txt ]; then
cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/dtool_datasets
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
echo --- > 1_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190716 >> 1_readme.yml
cat README.yml >> 1_readme.yml
dtool readme write SRX2190716 1_readme.yml
rm 1_readme.yml
echo --- > 2_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190717 >> 2_readme.yml
cat README.yml >> 2_readme.yml
dtool readme write SRX2190717 2_readme.yml
rm 2_readme.yml
echo --- > 3_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190718 >> 3_readme.yml
cat README.yml >> 3_readme.yml
dtool readme write SRX2190718 3_readme.yml
rm 3_readme.yml
echo --- > 4_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190719 >> 4_readme.yml
cat README.yml >> 4_readme.yml
dtool readme write SRX2190719 4_readme.yml
rm 4_readme.yml
echo --- > 5_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190720 >> 5_readme.yml
cat README.yml >> 5_readme.yml
dtool readme write SRX2190720 5_readme.yml
rm 5_readme.yml
echo --- > 6_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190721 >> 6_readme.yml
cat README.yml >> 6_readme.yml
dtool readme write SRX2190721 6_readme.yml
rm 6_readme.yml
echo --- > 7_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190722 >> 7_readme.yml
cat README.yml >> 7_readme.yml
dtool readme write SRX2190722 7_readme.yml
rm 7_readme.yml
echo --- > 8_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190723 >> 8_readme.yml
cat README.yml >> 8_readme.yml
dtool readme write SRX2190723 8_readme.yml
rm 8_readme.yml
echo --- > 9_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190724 >> 9_readme.yml
cat README.yml >> 9_readme.yml
dtool readme write SRX2190724 9_readme.yml
rm 9_readme.yml
echo --- > 10_readme.yml
echo description: Bisulphite sequencing, SRA accession SRX2190725 >> 10_readme.yml
cat README.yml >> 10_readme.yml
dtool readme write SRX2190725 10_readme.yml
rm 10_readme.yml
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/2-dtool_readme_0.txt
echo 2-dtool_readme_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
cd ..
sbatch 3-dtool_freeze_0.sh
else
echo 2-dtool_readme_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

