#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=0-29
#SBATCH --mem=16000

source trimmomatic-0.33

rawroot=0-raw_data/PRJEB2678
trimroot=2-trim/PRJEB2678

ARRAY=(ERR046546 ERR046547 ERR046548 ERR046549 ERR046550 ERR046551 ERR046552 ERR046553 ERR046554 ERR046555 ERR046556 ERR046557 ERR046558 ERR046559 ERR046560 ERR046561 ERR046562 ERR046563 ERR046564 ERR046565 ERR046566 ERR046567 ERR046568 ERR046569 ERR046570 ERR046571 ERR046572 ERR046573 ERR046574 ERR046575)

##Raw reads are 100bp long so using -ls 1 50 -rs 51 10
#srun bs-sequel.pl -f ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -l "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_1.fastq -r "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_2.fastq -b ${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -d "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -k 0 -n 2 -t Arabidopsis -2 0 -rnd 1 -mh 10 --parallel 0 --new-cm


srun java -jar ../../Software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 "${rawroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz "${rawroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_2.fastq.gz "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_1_unpaired.fastq.gz "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_2.fastq.gz "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_2_unpaired.fastq.gz ILLUMINACLIP:../../Software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

