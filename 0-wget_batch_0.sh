if [ ! -f /jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/0-wget_batch_0.txt ]; then
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/009/SRR4295279/SRR4295279_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/009/SRR4295279/SRR4295279_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/000/SRR4295280/SRR4295280.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/001/SRR4295281/SRR4295281.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/002/SRR4295282/SRR4295282.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/003/SRR4295283/SRR4295283.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/004/SRR4295284/SRR4295284.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/005/SRR4295285/SRR4295285.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/006/SRR4295286/SRR4295286_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/006/SRR4295286/SRR4295286_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/007/SRR4295287/SRR4295287_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/007/SRR4295287/SRR4295287_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/008/SRR4295288/SRR4295288.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/009/SRR4295289/SRR4295289.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/000/SRR4295290/SRR4295290_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/000/SRR4295290/SRR4295290_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/001/SRR4295291/SRR4295291_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/001/SRR4295291/SRR4295291_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/002/SRR4295292/SRR4295292_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR429/002/SRR4295292/SRR4295292_2.fastq.gz
echo complete>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/0-wget_batch_0.txt
echo 0-wget_batch_0 complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
else
echo 0-wget_batch_0 already_complete>>/jic/scratch/groups/Daniel-Zilberman/Projects/Jay-1001_methylomes/0-raw_data/staging_log.txt
fi

