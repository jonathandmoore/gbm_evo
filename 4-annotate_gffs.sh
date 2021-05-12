#!/bin/bash

# Call perl script to go through all the bs_sequel output GFF files in single-c and windows.
# For each file replace the 3rd column with informative identifiers for the sample and context

sbatch --mem=12000 --wrap 'perl Scripts/4-annotate_gffs.pl SRA035939'
sbatch --mem=12000 --wrap 'perl Scripts/4-annotate_gffs.pl PRJEB2678'

