#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=download
#SBATCH --time=02:00:00 # HH/MM/SS
#SBATCH --mem=8G

cd /athena/angsd/scratch/cnm4001/hg38
wget 'https://www.gencodegenes.org/human/' -O hg38.gencodev43.gtf.gz


