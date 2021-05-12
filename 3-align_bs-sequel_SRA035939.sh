#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=0-12
#SBATCH --mem=16000

source bssequel-0.0.1
source bowtie-1.0.1

trimroot=2-trim/SRA035939
alignroot=3-alignments/SRA035939

# bs-sequel.pl makes use of /tmp for sorting operations. This will rapidly exceed the available space on the execution node(s), so it makes sense to relocate temporary files
mkdir "${alignroot}"/tmp
TMPDIR= "${alignroot}"/tmp

ARRAY=(SRR342347 SRR342348 SRR342349 SRR342379 SRR342380 SRR342381 SRR342382 SRR342383 SRR342384 SRR342385 SRR342389 SRR342390 SRR342391 SRR342353 SRR342354 SRR342378)

##Raw reads are 100bp long so using -ls 1 50 -rs 51 10
srun bs-sequel.pl -f ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -l "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}.fastq -b ${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -d "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -k 0 -n 2 -t Arabidopsis -2 0 -rnd 1 -mh 10 --parallel 0 --new-cm

