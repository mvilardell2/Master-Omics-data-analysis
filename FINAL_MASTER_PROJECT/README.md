# Final Master Project
You can find here all the necessary information that aligns with the written report. This includes the supplementary material and the code scripts developed to perform the analysis. Due to data protection normative, I cannot upload the data used for the project.

## Title: Transcriptomic analysis of BRCA1 and BRCA2 carriers.
For my master's final project I performed a complex transcriptomic analysis from blood samples of individuals carrying BRCA1 and BRCA2 mutations. The identification of BRCA1 and BRCA2 carriers is essential, as it indicates which treatments and therapeutic approaches are more effective against cancer, based on the specific mutation of the individual. 
Gene expression profiling data was obtained using microarrays (Affymetrix GeneChhip Human).

Steps performed: 

Quality control of array images, annotation, deconvolution, filtering, variable selection, differential expression and gene set enrichment analysis. 

## Content: 
- Preprocessing steps report. This R markdown file contains: Loading data and exploration, quality control, data normalization, sample aggregation and PCA, and annotation (assign probe ID to gene).
- Deconvolution_final. Apply different methods of deconvolution, statistical tests performed to asses its variability in different conditions, correlation analysis between the 3 methods used for deconvolution.
- LASSO_DEG_GESA. Script that contains the filtering, variable selection with LASSO, differential expression gene analysis with Limma, and functional analysis.
- Final_excelfile_code. R script to generate a GTF file that contains detailed annotation information for each gene. It also generates an Excel file with the results of differential expression analysis for each of the comparisons that have been studied. 
