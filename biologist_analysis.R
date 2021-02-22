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
# data <- read.csv("/project/bf528/project_1/data/differential_expression_results.csv")
data <- read.csv("/projectnb2/bf528/users/lava_lamp/project_1/differential_expression_s5.csv", row.names = 1, header= TRUE)
data <- data %>% arrange(desc(t_statistics))
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
top10_down <- tail(data, 10)
top1000_down <- tail(data, 1000)
write.csv(top10_up, "top10_upregulated_genes.csv")
write.csv(top10_down, "top10_downregulated_genes.csv")

### 6.5 - Fisher Test ###

# Fisher Test Functions
fisher <- function(gmt){
  pvalues <- c()
  df <- list()
  for(geneset in gmt){
    setname <- setName(geneset)
    geneids <- geneIds(geneset)
    differentially_expressed <- length(data$SYMBOL)
    in_set <- length(geneids)
    in_set_differential <- sum(geneids %in% data$SYMBOL)
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
pvalues_kegg <- fisher(kegg)[[1]]
df_kegg <- fisher(kegg)[[2]]

# GO Pathways
pvalues_go <- fisher(go)[[1]]
df_go <- fisher(go)[[2]]

# Hallmark Pathways
pvalues_h <- fisher(hallmarks)[[1]]
df_h <- fisher(hallmarks)[[2]]

# Top 3 Enriched Gene Sets
top3_kegg <- names(head(sort(pvalues_kegg), 3))
keggs <- rbind(df_kegg[[top3_kegg[[1]]]], df_kegg[[top3_kegg[[2]]]], df_kegg[[top3_kegg[[3]]]])
# write.csv(keggs, "top_kegg.csv")
top3_go <- names(head(sort(pvalues_go), 3))
gos <- rbind(df_go[[top3_go[[1]]]], df_go[[top3_go[[2]]]], df_go[[top3_go[[3]]]])
# write.csv(gos, "top_go.csv")
top3_h <- names(head(sort(pvalues_h), 3))
hs <- rbind(df_h[[top3_h[[1]]]], df_h[[top3_h[[2]]]], df_h[[top3_h[[3]]]])
# write.csv(hs, "top_h.csv")

# Number of significantly enriched gene sets at p<0.05
n_enriched_kegg <- sum(pvalues_kegg < 0.05)
n_enriched_go <- sum(pvalues_go < 0.05)
n_enriched_h <- sum(pvalues_h < 0.05)

cat("Number of Enriched KEGG gene sets:", n_enriched_kegg) #186
cat("Number of Enriched GO gene sets:", n_enriched_go) #10271
cat("Number of Enriched Hallmark gene sets:", n_enriched_h) #50
