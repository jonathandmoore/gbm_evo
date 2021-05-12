#!/bin/bash

#mkdir 0-raw_data
#mkdir 0-raw_data/PRJEB2678
#mkdir 0-raw_data/SRA035939

mkdir 1-qc

mkdir 2-trim
mkdir 2-trim/PRJEB2678
mkdir 2-trim/SRA035939

mkdir 3-alignments
mkdir 3-alignments/PRJEB2678
mkdir 3-alignments/SRA035939

mkdir 4-summary
mkdir 4-summary/PRJEB2678
mkdir 4-summary/SRA035939

mkdir 5-analysis
mkdir 5-analysis/PRJEB2678
mkdir 5-analysis/SRA035939

#download PRJEB2678 fastq files from ENA
cd 0-raw_data/PRJEB2678
. PRJEB2678_download_fastqs_OSX.sh >PRJEB2678_download_fastqs_OSX.out 2>PRJEB2678_download_fastqs_OSX.err
cd ../..

#download SRA035939 fastq files from ENA
cd 0-raw_data/SRA035939
. SRA035939_download_fastqs_OSX.sh >SRA035939_download_fastqs_OSX.out 2>SRA035939_download_fastqs_OSX.err
cd ../..

