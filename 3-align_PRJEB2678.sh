#!/bin/bash

sample=$1

trimroot=2-trim/PRJEB2678/"${sample}"
alignroot=3-alignments/PRJEB2678/"${sample}"

rm -r "${alignroot}"
mkdir "${alignroot}"

#(single-end reads)

perl ../../Software/YAMA-master/yama.pl -x ../../Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.fasta -Q "${trimroot}".fastq.gz -o "${alignroot}" -n "${sample}" --bowtie2="-p 16 -k 50"

