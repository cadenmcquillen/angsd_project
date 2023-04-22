### Analysis of Next-Gen Sequencing Data Final Project Spring 2023

#### Scripts

#### Shell Scripts
Pipeline from fastq file to count matrix. Scripts located at ` ./angsd_project/project_scripts/Shell_scripts `
1. ##### Data Download
  1. ` download_fastq.sh `
  2. ` download_genome.sh `
  3. ` download_gtf.sh `
2. ##### Fastq QC
  4. ` fastqc.sh `
3. ##### Genome indexing and alginment
  5. ` index_star.sh `
  6. ` align_STAR.sh `
  7. ` index_alignments.sh `
4. ##### Alignment QC
  8. ` alignmentQC.sh `
  9. (optional) ` reRun_alignmentQC.sh `
5. ##### Feature Counts
  10. ` featureCounts.sh `
  11. (optional) ` reRun_featureCounts.sh `

#### Rmds
Rmarkdown files located at ` ./angsd_project/project_scripts/Rmd `
1. ##### Shell script documentation
  1. ` Download_Begin_Data_Process.Rmd `
  2. ` Visualizing_Alignments.Rmd `
2. ##### R based downstream analysis
  3. ` Feature_Counts.Rmd `
  4.` Project_Feddback.Rmd `
  5. ` Differential_expression.Rmd `
3. ##### Cumulative report
  6. ` Final_Report.Rmd `

