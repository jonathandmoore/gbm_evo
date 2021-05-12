#!/bin/bash

sample=$1

trimroot=2-trim/SRA035939/"${sample}"
alignroot=3-alignments/SRA035939/"${sample}"

rm -r "${alignroot}"
mkdir "${alignroot}"

#(single-end reads)

perl ../../Software/YAMA-master/yama.pl --keep-all-temp -x ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -Q "${trimroot}".fastq -o "${alignroot}" -n "${sample}" --bowtie2="-p 16 -k 10"

