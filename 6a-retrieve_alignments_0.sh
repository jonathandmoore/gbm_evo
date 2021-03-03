#!/bin/bash -e
#SBATCH -p jic-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=
#SBATCH --mem=64000
#SBATCH -o 6a-retrieve_alignments_0.out
#SBATCH -e 6a-retrieve_alignments_0.err
#SBATCH --localscratch=ssd:100    # request 100GB SSD
until [ -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/6-align_0.txt ]; do
  sleep 60
done
if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/6a-retrieve_alignments_0.txt ]; then
BASE_DIR=$SLURM_LOCAL_SCRATCH
export DTOOL_CACHE_DIRECTORY=$BASE_DIR
#iinit
source /jic/scratch/groups/Daniel-Zilberman/load_dtool.sh
input=dtool_datasets/12-dtool_copy_mc_0_uris.txt
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
mv $READ1_ABSPATH ../3-alignments/${this_item[1]}
done
done<"$input"
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/6a-retrieve_alignments_0.txt
echo 6a-retrieve_alignments_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
sbatch 8a-segmentation_0.sh
sbatch 6a-retrieve_alignments_1.sh
else
echo 6a-retrieve_alignments_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

