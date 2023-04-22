### Analysis of Next-Gen Sequencing Data Final Project Spring 2023

#### Scripts

#### Shell Scripts
Pipeline from fastq file to count matrix. Scripts located at ` ./angsd_project/project_scripts/Shell_scripts `
- ##### Data Download
* ` download_fastq.sh `
* ` download_genome.sh `
* ` download_gtf.sh `
- ##### Fastq QC
* ` fastqc.sh `
- ##### Genome indexing and alginment
* ` index_star.sh `
* ` align_STAR.sh `
* ` index_alignments.sh `
- ##### Alignment QC
* ` alignmentQC.sh `
* (optional) ` reRun_alignmentQC.sh `
- ##### Feature Counts
* ` featureCounts.sh `
* (optional) ` reRun_featureCounts.sh `

#### Rmds
Rmarkdown files located at ` ./angsd_project/project_scripts/Rmd `
- ##### Shell script documentation
* ` Download_Begin_Data_Process.Rmd `
* ` Visualizing_Alignments.Rmd `
- ##### R based downstream analysis
* ` Feature_Counts.Rmd `
* ` Project_Feddback.Rmd `
* ` Differential_expression.Rmd `
- ##### Cumulative report
* ` Final_Report.Rmd `

