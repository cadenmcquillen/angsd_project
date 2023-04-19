#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=align_qc
#SBATCH --time=00:30:00 # HH/MM/SS
#SBATCH --mem=8G

mamba activate angsd

featureCounts -a /athena/angsd/scratch/cnm4001/hg38/hg38.gencodev41.gtf \
-o /athena/angsd/scratch/cnm4001/NPC_counts/NPC_combined_featureCounts.txt \
/athena/angsd/scratch/cnm4001/NPC_bams/SRR11849597.Aligned.sortedByCoord.out.bam /athena/angsd/scratch/cnm4001/NPC_bams/SRR11849598.Aligned.sortedByCoord.out.bam
/athena/angsd/scratch/cnm4001/NPC_bams/SRR11849599.Aligned.sortedByCoord.out.bam
/athena/angsd/scratch/cnm4001/NPC_bams/SRR11849600.Aligned.sortedByCoord.out.bam
/athena/angsd/scratch/cnm4001/NPC_bams/SRR11849601.Aligned.sortedByCoord.out.bam
/athena/angsd/scratch/cnm4001/NPC_bams/SRR11849602.Aligned.sortedByCoord.out.bam
