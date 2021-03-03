#!/bin/bash -e
#SBATCH -p jic-long # partition (queue)
#SBATCH -t 30-00:00
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=64000

cd /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-reference

# get the data (takes ages - it's 132GB)
#wget  https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf

# parse the data into single files of modified cytosones per sample. Also takes ages
perl ../Scripts/parse_1001genome_VCF.pl <1001genomes_snp-short-indel_only_ACGTN.vcf >parse_1001genome_VCF.log

