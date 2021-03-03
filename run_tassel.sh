#!/bin/bash -e
#SBATCH -p jic-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jonathan.moore@jic.ac.uk # send-to address
#SBATCH --array=2-199
#SBATCH --mem=32000

source tassel-5.0

base_dir=/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/tassel5-standalone/881

mkdir ${base_dir}/${SLURM_ARRAY_TASK_ID}
cd ${base_dir}/${SLURM_ARRAY_TASK_ID}
cp ../200_envs_for_tassel_${SLURM_ARRAY_TASK_ID}.txt ./climate_variable_for_tassel.txt

perl ../../run_pipeline.pl  -configFile ../881_mlm_config_template.xml

