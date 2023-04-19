#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=index
#SBATCH --time=00:20:00 # HH/MM/SS
#SBATCH --mem=8G

cd /athena/angsd/scratch/cnm4001/
mamba activate angsd


for index in 597 598 599 600 601 602; do

  samtools index /athena/angsd/scratch/cnm4001/NPC_bams/SRR11849${index}.Aligned.sortedByCoord.out.bam

done
