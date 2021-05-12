#!/bin/bash

../../Software/FastQC.app/Contents/MacOS/fastqc -o 1-qc 0-raw_data/PRJEB2678/ERR*.fastq.gz
../../Software/FastQC.app/Contents/MacOS/fastqc -o 1-qc 0-raw_data/SRA035939/SRR*.fastq.gz
