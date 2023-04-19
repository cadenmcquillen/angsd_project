#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=test
#SBATCH --time=06:0:00 # HH/MM/SS
#SBATCH --mem=6G



cd /athena/angsd/scratch/cnm4001/NPC_fastq_data


while read ID; do
  wget "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=${ID}" -O ${ID}.fastq.gz 


done < SRR_ACC_List.txt
