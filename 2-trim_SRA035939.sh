#!/bin/bash

sample=$1

rawroot=0-raw_data/SRA035939/"${sample}"
trimroot=2-trim/SRA035939/"${sample}"

java -jar ../../Software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 "${rawroot}".fastq.gz "${trimroot}".fastq.gz ILLUMINACLIP:../../Software/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

