# Microarray Based Tumor Classification

This is the repository for Project 1 in BF528, Spring 2021.
This analysis will focus on reproducing the results from the comparison of C3 and C4 tumor subtypes, as demonstrated in the following paper:

Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391

## Contributors

- Daisy Wenyan Han daisyhan@bu.edu
- Divya Sundaresan divyas3@bu.edu
- Alec Jacobsen aggjacob@bu.edu
- Emmanuel Saake esaake@bu.edu

## Repository Contents

* __data_preprocessing_and_QC.R__ -- R script to normalize all of the microarrays together using the Robust Multiarray Averaging (RMA) algorithm. Standard quality control metrics are computed on the normalized data, and the distribution of samples is visualized using Principal Component Analysis (PCA).
* __analyst_p1.R__ --Analyst script to perform noise and data dimensionality reduction in three phases, a logarithmic test involving the selection of genes on the criteria that 20% of their gene-expression value be > log2(15)); a low variance filter with a threshold of p<0.01 (Involving the computation of critical region using chi-squared); and coefficient of variation (> 0.186) test. Further analysis including hierarchical clustering, Welch t-test evaluation, and then  p_value adjustment were undertaken. 
* __biologist_analysis.R__ -- The authors in Marisa et al. sought to understand the biological significance of the different gene expression profiles for each tumor subtype using gene set enrichment analysis. Specifically, KEGG, GO, and cancer hallmark genesets were compared with the top 1000 up- and down-regulated genes of each subtype compared with all the others using Fisher’s Exact test. This analysis was reproduced using KEGG gene sets and the differential expression results generated above.
