#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=download
#SBATCH --time=02:00:00 # HH/MM/SS
#SBATCH --mem=8G

cd /athena/angsd/scratch/cnm4001/hg38
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz' -O hg38.fa.gz
wget 'https://genome.ucsc.edu/cgi-bin/hgTables' -O hg38.sgd.gtf.gz


