#!/bin/bash -e
#SBATCH -p jic-long # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=64000
#SBATCH -o 6-align_0_fix.out
#SBATCH -e 6-align_0_fix.err
#SBATCH --localscratch=ssd:100    # request 100GB SSD
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/5-dtool_tidy_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/6-align_0.txt ]; then
BASE_DIR=$SLURM_LOCAL_SCRATCH
export DTOOL_CACHE_DIRECTORY=$BASE_DIR
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
source  bsmap-2.90.0
source samtools-1.9
input=dtool_datasets/4-dtool_copy_0_uris.txt
while IFS='\n' read -r line
do
this_dataset=($line)
#echo ${this_dataset[1]}
dtool ls ${this_dataset[1]} | while read item; do
echo $item
this_item=($item)
#echo ${this_item[0]}
#echo ${this_item[1]}
READ1_ABSPATH=$(dtool item fetch ${this_dataset[1]} ${this_item[0]})
#echo $READ1_ABSPATH
mv $READ1_ABSPATH $BASE_DIR/${this_item[1]}
done
done<"$input"
bsmap -a $BASE_DIR/SRR4295284.fastq.gz -d /jic/scratch/groups/Daniel-Zilberman/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -o $BASE_DIR/SRR4295284_TAIR10.bsp
python ../Scripts/methratio.py -z -g -o $BASE_DIR/SRR4295284_TAIR10.methratio2 -d /jic/scratch/groups/Daniel-Zilberman/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta $BASE_DIR/SRR4295284_TAIR10.bsp
head -n 1 $BASE_DIR/SRR4295284_TAIR10.methratio2 > $BASE_DIR/SRR4295284_TAIR10.cg.methratio2
awk '{if (substr($4,3,2) == "CG") { print } }' $BASE_DIR/SRR4295284_TAIR10.methratio2 >> $BASE_DIR/SRR4295284_TAIR10.cg.methratio2
mv $BASE_DIR/SRR4295284_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/
mv $BASE_DIR/SRR4295284_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/
rm $BASE_DIR/SRR4295284_TAIR10.*
bsmap -a $BASE_DIR/SRR4295285.fastq.gz -d /jic/scratch/groups/Daniel-Zilberman/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -o $BASE_DIR/SRR4295285_TAIR10.bsp
python ../Scripts/methratio.py -z -g -o $BASE_DIR/SRR4295285_TAIR10.methratio2 -d /jic/scratch/groups/Daniel-Zilberman/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta $BASE_DIR/SRR4295285_TAIR10.bsp
head -n 1 $BASE_DIR/SRR4295285_TAIR10.methratio2 > $BASE_DIR/SRR4295285_TAIR10.cg.methratio2
awk '{if (substr($4,3,2) == "CG") { print } }' $BASE_DIR/SRR4295285_TAIR10.methratio2 >> $BASE_DIR/SRR4295285_TAIR10.cg.methratio2
mv $BASE_DIR/SRR4295285_TAIR10.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/
mv $BASE_DIR/SRR4295285_TAIR10.cg.methratio2 /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/../3-alignments/
rm $BASE_DIR/SRR4295285_TAIR10.*

echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/6-align_0.txt
echo 6-align_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
#sbatch 7-call_methylation_0.sh
else
echo 6-align_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

