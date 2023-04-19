#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=align_qc
#SBATCH --time=03:00:00 # HH/MM/SS
#SBATCH --mem=20G


cd /athena/angsd/scratch/cnm4001/NPC_alignmentQC
mkdir ReRun_Qorts
cd ReRun_Qorts
for index in 597 598 599 600 601 602; do
  mkdir reRun_SRR11849${index}_alignmentQC
done

cd /athena/angsd/scratch/cnm4001/NPC_bams
mamba activate qorts

for index in 597 598 599 600 601 602; do
  qorts -Xmx18000M QC \
  --generatePlots \
  --singleEnded \
   SRR11849${index}.Aligned.sortedByCoord.out.bam \
  /athena/angsd/scratch/cnm4001/hg38/hg38.gencodev43.gtf  \
  /athena/angsd/scratch/cnm4001/NPC_alignmentQC/ReRun_Qorts/reRun_SRR11849${index}_alignmentQC
done


