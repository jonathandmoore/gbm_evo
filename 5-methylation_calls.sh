#!/bin/bash

source R-3.4.0-20170516

# NB. System R did not contain the optparse package, and user has no permissions to install packages on cluster, nor a link to CRAN
# It was installed locally from the CRAN package by something like this (but wget was done on PC first, then the file copied to the user's R directory):

#cd ~/R/x86_64-pc-linux-gnu-library/3.4/
#wget https://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz
#wget https://cran.r-project.org/src/contrib/optparse_1.4.4.tar.gz
#R CMD INSTALL getopt_1.20.0.tar.gz
#R CMD INSTALL optparse_1.4.4.tar.gz

# Call Rscript to run the R script for the test project, CG context
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a load'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a plot'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a call'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a merge'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a analyse'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a all'
#sbatch --mem=64000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CHG -a all'
sbatch --mem=64000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CHG -a plot'
sbatch --mem=64000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CHG -a call'

sbatch --mem=64000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CHH -a all'

# Call Rscript to run the R script for the Schmitz project, CHG context
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p SRA035939 -c CG -a load'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p SRA035939 -c CHG -a plot'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p SRA035939 -c CHG -a call'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p SRA035939 -c CHG -a merge'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p SRA035939 -c CHG -a analyse'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p SRA035939 -c CHG -a all'

# Call Rscript to run the R script for the test project, CG context
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a load'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a plot'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a call'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a merge'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a analyse'
#sbatch --mem=24000 --wrap 'Rscript Scripts/5-methylation_calls.R -p test -c CG -a all'

