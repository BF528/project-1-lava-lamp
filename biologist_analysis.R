# Load Required Packages
library(tidyverse)
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(GSEABase)

### 6.1 - Map Probe ID's ###

# Read in data
# data <- read.csv("/project/bf528/project_1/data/differential_expression_results.csv") #Sample data
# data <- read.csv("/projectnb2/bf528/users/lava_lamp/project_1/differential_expression_s5.csv", row.names = 1, header= TRUE) #Passed all filters
data <- read.csv("/projectnb2/bf528/users/lava_lamp/project_1/diffexp_s5_6.csv", row.names = 1, header= TRUE) #Unfiltered
data <- data %>% 
  filter(adjusted_p < 0.05) %>% #Filter for FDR < 0.05
  arrange(desc(t_statistics))
matches <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(row.names(data)), columns = ("SYMBOL"))

# Remove duplicates
dedup_matches <- matches[!duplicated(matches[1]), ]

# Combine symbols with differential expression results
data <- cbind(dedup_matches, data)

# Multiple probes may align to the same gene symbol. Remove duplicates, keeping those with most significant adjusted p-value
data <- data %>%
  group_by(SYMBOL) %>%
  filter(adjusted_p == min(adjusted_p)) %>%
  ungroup(SYMBOL)

### 6.2 and 6.4 - Gene Sets ###

# Load in Gene Sets
hallmarks <- getGmt("genesets/h.all.v7.2.symbols.gmt")
kegg <- getGmt("genesets/c2.cp.kegg.v7.2.symbols.gmt")
go <- getGmt("genesets/c5.go.v7.2.symbols.gmt")

# Gene Set Lengths
cat("Hallmarks Length: ", length(names(hallmarks))) #50
cat("GO Length: ", length(names(go))) #10271
cat("KEGG Length: ", length(names(kegg))) #186

## 6.3 - Top Differentially Expressed Genes
top10_up <- head(data, 10)
top1000_up <- head(data, 1000)
top10_down <- head(data %>% arrange(t_statistics), 10)
top1000_down <- head(data %>% arrange(t_statistics), 1000)
# write.csv(top10_up, "top10_upregulated_genes.csv")
# write.csv(top10_down, "top10_downregulated_genes.csv")

### 6.5 - Fisher Test ###

# Fisher Test Functions
fisher <- function(gmt, dataset){
  pvalues <- c()
  df <- list()
  for(geneset in gmt){
    setname <- setName(geneset)
    geneids <- geneIds(geneset)
    differentially_expressed <- length(dataset$SYMBOL)
    in_set <- length(geneids)
    in_set_differential <- sum(geneids %in% dataset$SYMBOL)
    in_set_not_differential <- in_set - in_set_differential
    not_in_set_differential <- differentially_expressed - in_set_differential
    not_in_set_not_differential <- 0
    fishervals <- fisher.test(matrix(c(in_set_differential, in_set_not_differential, not_in_set_differential, not_in_set_not_differential), nrow = 2))
    pval <- fishervals$p.value
    est <- fishervals$estimate
    padj <- p.adjust(pval, method="fdr")
    df[[setname]] <- data.frame(geneset = setname, statistic = est, pvalue = pval, p.adj = padj)
    pvalues[setname] <- pval
  }
  rets <- list(pvalues, df)
  return(rets)
}

# KEGG Pathways
pvalues_kegg_up <- fisher(kegg, top1000_up)[[1]]
df_kegg_up <- fisher(kegg, top1000_up)[[2]]
pvalues_kegg_down <- fisher(kegg, top1000_down)[[1]]
df_kegg_down <- fisher(kegg, top1000_down)[[2]]

# GO Pathways
pvalues_go_up <- fisher(go, top1000_up)[[1]]
df_go_up <- fisher(go, top1000_up)[[2]]
pvalues_go_down <- fisher(go, top1000_down)[[1]]
df_go_down <- fisher(go, top1000_down)[[2]]

# Hallmark Pathways
pvalues_h_up <- fisher(hallmarks, top1000_up)[[1]]
df_h_up <- fisher(hallmarks, top1000_up)[[2]]
pvalues_h_down <- fisher(hallmarks, top1000_down)[[1]]
df_h_down <- fisher(hallmarks, top1000_down)[[2]]

# Top 3 Enriched Gene Sets
# KEGG Pathways
top3_kegg <- names(head(sort(pvalues_kegg_up), 3))
keggs_up <- rbind(df_kegg_up[[top3_kegg[[1]]]], df_kegg_up[[top3_kegg[[2]]]], df_kegg_up[[top3_kegg[[3]]]])
write.csv(keggs_up, "top_kegg_up.csv")
top3_kegg <- names(head(sort(pvalues_kegg_down), 3))
keggs_down <- rbind(df_kegg_down[[top3_kegg[[1]]]], df_kegg_down[[top3_kegg[[2]]]], df_kegg_down[[top3_kegg[[3]]]])
write.csv(keggs_down, "top_kegg_down.csv")

# GO Pathways
top3_go <- names(head(sort(pvalues_go_up), 3))
gos_up <- rbind(df_go_up[[top3_go[[1]]]], df_go_up[[top3_go[[2]]]], df_go_up[[top3_go[[3]]]])
write.csv(gos_up, "top_go_up.csv")
top3_go <- names(head(sort(pvalues_go_down), 3))
gos_down <- rbind(df_go_down[[top3_go[[1]]]], df_go_down[[top3_go[[2]]]], df_go_down[[top3_go[[3]]]])
write.csv(gos_down, "top_go_down.csv")

# Hallmarks 
top3_h <- names(head(sort(pvalues_h_up), 3))
hs_up <- rbind(df_h_up[[top3_h[[1]]]], df_h_up[[top3_h[[2]]]], df_h_up[[top3_h[[3]]]])
write.csv(hs_up, "top_h_up.csv")
top3_h <- names(head(sort(pvalues_h_down), 3))
hs_down <- rbind(df_h_down[[top3_h[[1]]]], df_h_down[[top3_h[[2]]]], df_h_down[[top3_h[[3]]]])
write.csv(hs_down, "top_h_down.csv")

# Number of significantly enriched gene sets at p<0.05
n_enriched_kegg_up <- sum(pvalues_kegg_up < 0.05)
n_enriched_kegg_down <- sum(pvalues_kegg_down < 0.05)
n_enriched_go_up <- sum(pvalues_go_up < 0.05)
n_enriched_go_down <- sum(pvalues_go_down < 0.05)
n_enriched_h_up <- sum(pvalues_h_up < 0.05)
n_enriched_h_down <- sum(pvalues_h_down < 0.05)

cat("Number of Upregulated Enriched KEGG gene sets:", n_enriched_kegg_up) #186
cat("Number of Downregulated Enriched KEGG gene sets:", n_enriched_kegg_down) #186
cat("Number of Upregulated Enriched GO gene sets:", n_enriched_go_up) #10271
cat("Number of Downregulated Enriched GO gene sets:", n_enriched_go_down) #10271
cat("Number of Upregulated Enriched Hallmark gene sets:", n_enriched_h_up) #50
cat("Number of Downregulated Enriched Hallmark gene sets:", n_enriched_h_down) #50
