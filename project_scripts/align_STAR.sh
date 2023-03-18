#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=align
#SBATCH --time=4:00:00 # HH/MM/SS
#SBATCH --mem=32G

cd /athena/angsd/scratch/cnm4001/
mamba activate angsd


for index in 597 598 599 600 601 602; do


  STAR --runMode alignReads \
    --runThreadN 1 \
    --genomeDir  hg38/hg38_STARindex \
    --readFilesIn  NPC_fastq_data/SRR11849${index}.fastq.gz\
    --readFilesCommand zcat \
    --outFileNamePrefix NPC_bams/SRR11849${index}. \
    --outSAMtype BAM SortedByCoordinate ;

  samtools index /athena/angsd/scratch/cnm4001/NPC_bams/SRR11849${index}.Aligned.sortedByCoord.out.bam
done
