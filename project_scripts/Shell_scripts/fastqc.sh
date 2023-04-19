#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=test
#SBATCH --time=01:00:00 # HH/MM/SS
#SBATCH --mem=6G


cd /athena/angsd/scratch/cnm4001/NPC_fastq_data/
mamba activate angsd

for file in *; do
  fastqc $file --extract;
done
