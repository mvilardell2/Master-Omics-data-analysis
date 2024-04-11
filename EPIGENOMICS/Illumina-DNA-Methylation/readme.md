# Illumina DNA methylation data analysis

In this folder you'll find all the information regarding Illumina Infinium 450K and EPIC DNA methylation analysis.

Beta-values-->how methylated is a specific CpG island (not working in single-cell experiment).

## Content
- Script_EPIC.R: Contains the steps to analyse DNA methylation data coming from EPIC arrays.
- Script_450K.R: Steps to analyse DNA methylation data coming from 450K arrays.
- Practical_Vilardell.Rmd: Rmarkdown that contains the final practice for my master's assignment.

## Workflow:
- Read raw methylation data: from IDATS and Sample sheet (details of the samples).
- Quality control: Identify outliers samples.
  Include signal detection, quality check, density plots of beta values, detection p-values, and filtering unwanted probes and samples.
- Normalization with minfi (default NOOB normalization). Density plot of beta-values after normalization.
- Exploratory data analysis: Check for correlation variables in the data.
- Heatmaps of top variable CpG probes based on norm-beta-values.
- Create the model
- Differentially methylated positions--> Limma. Extract DMPs per contrast. Annotate Genomic context and genic regions. 
- Differentially methylated regions--> DMRcate. Minimum of 3-5 cpG. Annotation to the nearest gene.
- Functional analysis of DMP and DMRs with Gprofiler2
 





