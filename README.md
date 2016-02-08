# RefCNV

<b>RefCNV</b> is a R script with novel approach for the identification of copy number variants from whole exome sequencing data. The algorithm is empirically based, and allows for local adjustment for differences in sequence read coverage across genomic regions. The method is based on establishing a set of normal reference controls.

# RefCNV
<b>Title</b> RefCNV: Identification of gene-based copy number variants using whole exome sequencing

<b>Date</b> 2016-1-13

<b>Author</b> Lun-Ching Chang

<b>Maintainer</b> Lun-Ching Chang lun-ching.chang@nih.gov; lunching@gmail.com

<b>Description</b> <b>RefCNV</b> is a R toolkit for copy number variants using whole exome sequencing data. Current pipeline includes two major compartments: (1) run regression model for all provided replicated references and perform LOOCV (leave-one-out cross validation) for CNV threshold; (2) summarize exom-based results into gene-based and predict CNVs (A: amplification, N: normal and D: deletion).

<b>Depends</b> R (>=3.2.2)

# Main function of RefCNV
Before running our CNV algorithm, all users have to run "bedtools" (http://bedtools.readthedocs.org/en/latest/) to generate the coverage data. Sample code for running coverage based on the capture regions in ".bed" file. (please set the output file name end with ".all.coverage")

bedtools coverage -abam file_name.bam -b file_name.bed > file_name.all.coverage

## [GB_CNVs_EL.R](https://github.com/lunching/RefCNV.git) 

<b>Usage</b>
* GB_CNVs_EL(case, ref, CNV_q = .05, avg_cov_rm = 30, TMR, TMR_new, gene_map)

<b>Arguments</b>
* case 
  * a list of cases generated by function <b>Cov_list</b> (see sample code [Cov_list](https://github.com/lunching/RefCNV.git) here).
* ref 
  * a list of references generated by function xxx.r.
* CNV_q
  * quantiles of the median standardized residuals across all references from whole genome (default = .05).
* avg_cov_rm
  * the minimum coverage (default = 30) that the capture region will retain for futher analysis.
* TMR
  * numeric vector of total number of mapped read for each reference (the length of the vector need to match the number of references). 
* TMR_new
  * numeric vector of total number of mapped read for each case (the length of the vector need to match the number of cases). 
* gene_map
  * a matrix of gene map file need to provided by user which includes gene symbol, gene ID, chromosome, start position and end position (see sample map file [Sample_gene.map](https://github.com/lunching/RefCNV.git)).

<b>Output</b>

Outputs will be automatically generated by function "GB_CNVs_EL" into excel file with name "CNVs_summary.xls" which contains gene symbol, median of standardized residuals of all cases followed by colnames sMSR_S with constitutive numbers and CNV prediction (D: deletion, N: normal, and A: amplification) of all cases followed by colnames CNV_S with constitutive numbers (see sample output [CNVs_summary.xls](https://github.com/lunching/RefCNV.git)). 
