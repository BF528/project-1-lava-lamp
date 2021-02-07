rm(list = ls())
setwd('/projectnb/bf528/users/lava_lamp/project_1')

# installing and loading required packages 
BiocManager::install(c('affy','affyPLM','sva','AnnotationDbi','hgu133plus2.db'))
install.packages("ggplot2")
install.packages("ggpubr")
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(ggarrange)
library(ggpubr)

# reading in CEL files and normalizing together 
data <- ReadAffy(celfile.path = '/projectnb/bf528/users/lava_lamp/project_1/samples')
eset <- rma(data)

# computing RLE and NUSE scores and plotting hist of means
Pset <- fitPLM(data, normalize = TRUE, background = TRUE)
medianRLE <- apply(RLE(Pset,type = 'values'),2,median)
medianNUSE <- apply(NUSE(Pset,type = 'values'),2,median)

RLEplot <- ggplot(as.data.frame(medianRLE), aes(x = medianRLE)) + geom_histogram(fill = 'cornflowerblue', color = 'black')+
              theme_light() + ylab('Count') + xlab('RLE') + ggtitle('Median RLE for Samples')
NUSEplot <- ggplot(as.data.frame(medianNUSE), aes(x = medianNUSE)) + geom_histogram(fill = 'darkorange', color = 'black')+
              theme_light() + ylab('Count') + xlab('NUSE') + ggtitle('Median NUSE for Samples')
pdf('Median_RLE_NUSE_hists.pdf')
ggarrange(RLEplot,NUSEplot, nrow = 2)
dev.off()

# correcting for batch effects and outputting CSV
metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
mod <- model.matrix(~as.factor(normalizationcombatmod), data = metadata)
combat_data <- ComBat(dat = exprs(eset), batch = metadata$normalizationcombatbatch, mod = mod)
write.csv(combat_data, file = 'expression_data.csv')

# scaling data and PCA
scaled_data <- t(scale(t(combat_data)))
pca <- prcomp(scaled_data, scale = F, center = F)

summary(pca) #to see % variance explained
pc12 <- as.data.frame(pca$rotation[,1:2]) #extracting first two PCs for plotting

pdf('Expression_PCA.pdf')
ggplot(data = pc12, aes(PC1,PC2)) + geom_point(color = 'darkred')+
  theme_light()+xlab('PC1 (11.47%)')+ylab('PC2 (8.41%)') + ggtitle('PCA of Gene Expression for 134 Samples')
dev.off()
