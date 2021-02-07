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

* biologist_analysis.R -- The authors in Marisa et al. sought to understand the biological significance of the different gene expression profiles for each tumor subtype using gene set enrichment analysis. Specifically, KEGG, GO, and cancer hallmark genesets were compared with the top 1000 up- and down-regulated genes of each subtype compared with all the others using Fisher’s Exact test. We will try to reproduce this analysis using KEGG gene sets and the differential expression results generated above.
