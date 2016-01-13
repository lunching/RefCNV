# RefCNV
<b>RefCNV</b> is a novel approach for the identification of copy number variants from whole exome sequencing data. The algorithm is empirically based, and allows for local adjustment for differences in sequence read coverage across genomic regions. The method is based on establishing a set of normal reference controls.

# RefCNV
<b>Title</b> RefCNV: Identification of gene-based copy number variants using whole exome sequencing

<b>Date</b> 2016-1-13

<b>Description</b> <b>RefCNV</b> is a R toolkit for copy number variants using whole exome sequencing data. Current pipeline includes two major compartments: (1) run regression model for all provided replicated references and perform LOOCV (leave-one-out cross validation) for CNV threshold; (2) summarize exom-based results into gene-based and predict CNVs (A: amplification, N: normal and D: deletion).

<b>Depends</b> R (>=3.2.2)

# Main function of RefCNV
Before running our CNV algorithm, all users have to run "bedtools" (http://bedtools.readthedocs.org/en/latest/) to generate the coverage data. Sample code for running coverage based on the capture regions in ".bed" file. (please set the output file name end with ".all.coverage")

bedtools coverage -abam file_name.bam -b file_name.bed > file_name.all.coverage

## GB_CNVs_EL 

<b>Usage</b>
* GB_CNVs_EL(case, ref, CNV_q = .05, avg_cov_rm = 30, TMR, TMR_new)

** Arguments **
