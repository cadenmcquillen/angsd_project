#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=download
#SBATCH --time=02:00:00 # HH/MM/SS
#SBATCH --mem=8G

cd /athena/angsd/scratch/cnm4001/hg38
wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz' -O hg38.gencodev43.gtf.gz


