#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=index
#SBATCH --time=24:00:00 # HH/MM/SS
#SBATCH --mem=64G


cd /athena/angsd/scratch/cnm4001/hg38
mamba activate angsd

STAR --runMode genomeGenerate \
  --runThreadN 1 \
  --genomeDir hg38_STARindex \
  --genomeFastaFiles /athena/angsd/scratch/cnm4001/hg38/hg38.fa \
  --sjdbGTFfile /athena/angsd/scratch/cnm4001/hg38/hg38.gencodev41.gtf \
  --sjdbOverhang 49
