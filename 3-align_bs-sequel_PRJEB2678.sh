#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jmoore@jic.ac.uk # send-to address
#SBATCH --array=0-29
#SBATCH --mem=16000

source bssequel-0.0.1
source bowtie-1.0.1

#trimroot=2-trim/PRJEB2678
trimroot=0-raw_data/PRJEB2678
alignroot=3-alignments/PRJEB2678

# bs-sequel.pl makes use of /tmp for sorting operations. This will rapidly exceed the available space on the execution node(s), so it makes sense to relocate temporary files
mkdir "${alignroot}"/tmp
TMPDIR= "${alignroot}"/tmp

#ARRAY=(ERR046572)
#ARRAY=(ERR046546 ERR046547 ERR046548)
ARRAY=(ERR046546 ERR046547 ERR046548 ERR046549 ERR046550 ERR046551 ERR046552 ERR046553 ERR046554 ERR046555 ERR046556 ERR046557 ERR046558 ERR046559 ERR046560 ERR046561 ERR046562 ERR046563 ERR046564 ERR046565 ERR046566 ERR046567 ERR046568 ERR046569 ERR046570 ERR046571 ERR046572 ERR046573 ERR046574 ERR046575)

#Paired-end version using trimmed data
#srun bs-sequel.pl -f ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -l "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_1.fastq -r "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_2.fastq -b ${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -d "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -k 300 -n 2 -t Arabidopsis -2 1 -rnd 1 -mh 10 --parallel 0 --new-cm
# This failed and complained that the ends are different lengths

#Paired-end version using untrimmed data - Raw reads are 100bp long so using -ls 1 50 -rs 51 10
#srun bs-sequel.pl -f ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -l "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_1.fastq -r "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_2.fastq -b ${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -d "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -k 300 -n 2 -t Arabidopsis -ls 1 50 -rs 51 100 -2 1 -rnd 1 -mh 10 --parallel 0 --new-cm
# This aligned successfully but failed to stitch the parts together with a handfull of reads at a partcular point on CHR4

#Single-end version using trimmed data
#srun bs-sequel.pl -f ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -l "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}.fq -b ${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -d "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -k 0 -n 2 -t Arabidopsis -2 0 -rnd 1 -mh 10 --parallel 0 --new-cm
# This ran successfully but only aligned around 50% reads

#Single-end version using untrimmed data - Raw reads are 100bp long so using -ls 1 50 -rs 51 10
srun bs-sequel.pl -f ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -l "${trimroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}.fq -b ${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -d "${alignroot}"/${ARRAY[$SLURM_ARRAY_TASK_ID]}_TAIR10 -k 0 -n 2 -t Arabidopsis -ls 1 50 -rs 51 100 -2 0 -rnd 1 -mh 10 --parallel 0 --new-cm
