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
* __analyst_p1.R__ -- Marisa et al. selected genes based on three well defined metrics. These metrics are computed for the above normalized data. Further filtering using suggested cutoffs resulted in a reduced set of features as described in the supplemental methods of the paper. 
* __biologist_analysis.R__ -- The authors in Marisa et al. sought to understand the biological significance of the different gene expression profiles for each tumor subtype using gene set enrichment analysis. Specifically, KEGG, GO, and cancer hallmark genesets were compared with the top 1000 up- and down-regulated genes of each subtype compared with all the others using Fisher’s Exact test. This analysis was reproduced using KEGG gene sets and the differential expression results generated above.
